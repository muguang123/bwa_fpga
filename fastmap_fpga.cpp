//
// Created by jie on 3/9/19.
//

#include <zlib.h>
#include "fastmap_fpga.h"
#include "ConcurrentMap_Grampa.h"
#include "kseq.h"
#include "threadpool.h"
#include "concurrentqueue.h"
#include "fpga_sw.h"
#include "ksw.h"
#include "kvec.h"
#include "utils.h"

KSEQ_DECLARE(gzFile)
using namespace std;

extern "C"{     
    void smem_aux_destroy(smem_aux_t *a);
    extern smem_aux_t *smem_aux_init();
    int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, int n, mem_alnreg_t *a);
    int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);
    void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m);
    int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2]);

}

typedef struct {
    kseq_t *ks, *ks2;
    mem_opt_t *opt;
    mem_pestat_t *pes0;
    int64_t n_processed;
    int copy_comment, actual_chunk_size;
    bwaidx_t *idx;
} ktp_aux_t;

typedef struct {
    ktp_aux_t *aux;
    int n_seqs;
    bseq1_t *seqs;
} ktp_data_t;

//============global data===============

typedef moodycamel::BlockingConcurrentQueue<ktp_data_t*> BlockingQueueStep0;
static BlockingQueueStep0 blockingQueueStep0(BLOCKQUEUE_STEP0_SIZE); //queue0 ret


typedef moodycamel::BlockingConcurrentQueue<fpga_sw_output*> BlockingQueueStep2;
static BlockingQueueStep2 blockingQueueStep2(BLOCKQUEUE_STEP2_SIZE); //queue0

extern BlockingQueueStep1 blockingQueueStep1;

typedef moodycamel::BlockingConcurrentQueue<TaskReadInputType*> BlockingQueueStep3;
static BlockingQueueStep3 blockingQueueStep3(BLOCKQUEUE_STEP3_SIZE); //queue0


typedef junction::ConcurrentMap_Grampa<int,TaskReadInputType*> ConcurrentMap;
static ConcurrentMap concurrentMap(32);//hashmap

static ktp_aux_t aux_glb;//全局变量，用于process0, process1  SOFT_SW, process3

std::atomic_int mapTaskNum(1);
std::atomic_int mapTaskCount(1);
std::mutex mtx_readInputFile;
std::mutex mtx_writeToFpga;
std::mutex mtx_readFromFpga;
std::mutex mtx_erase_task;
std::mutex mtx_n_processed;
int file_finish_flag = 0;

int Fpga_flag = 0; //0 stand for no using FPGA
//===============================================

//===================First Step================

static void *process0(void *procees0_args)
{
    while(true){

        ktp_aux_t *aux = (ktp_aux_t*)procees0_args;
        int i;
        // process0 read the input file
        ktp_data_t *ret;
        int64_t size = 0;
        ret = (ktp_data_t*) calloc(1, sizeof(ktp_data_t)); 
        mtx_readInputFile.lock();//lock the bseq_read() function
        ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, aux->ks2);
        mtx_readInputFile.unlock();
        if(ret->seqs == 0)
            file_finish_flag = 1;
        if (ret->seqs == 0) {
            free(ret);
            return 0;
        }
        if (!aux->copy_comment)
            for (i = 0; i < ret->n_seqs; ++i) {
                free(ret->seqs[i].comment);
                ret->seqs[i].comment = 0;
            }
        for (i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
        if (bwa_verbose >= 3)
            fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, ret->n_seqs, (long)size);
        blockingQueueStep0.enqueue(ret); //

    }
}

//=====================Second Step=====================

static void process1(void *data)
{// 第二步：调用mem_process_seqs对seqs进行处理，直到写入SAM文件
    ktp_data_t* point;
    while(true){

        blockingQueueStep0.wait_dequeue(point);
        const mem_opt_t *opt = aux_glb.opt;// 用户配置的opt
        const bwaidx_t *idx = aux_glb.idx;// 读取的idx
        int read_num, i =0;

        read_num = point->n_seqs;
        TaskReadInputType* task_input;
        task_input = (TaskReadInputType*)malloc(sizeof(TaskReadInputType));//
        task_input->one_read =(OneReadInputType**)malloc(read_num*sizeof(OneReadInputType*));
        int taskId = mapTaskNum++;  //todo better to use %
        smem_aux_t * buf;
        buf = smem_aux_init();
        int num_seed = 0;// the number of seed per task
        for (i=0; i<read_num; i++){
            char *seq = point->seqs[i].seq;
            int l_seq = point->seqs[i].l_seq;
            mem_align1_core_pre(opt, idx, l_seq, seq, task_input, i,  &num_seed, buf);//会入队并且有hashmap
        }
        smem_aux_destroy(buf);
        task_input->task_id = taskId;
        task_input->num_seed = num_seed;
        task_input->count = 0;
        task_input->read_num = read_num;
        task_input->seqs = point->seqs;
        task_input->n_processed = aux_glb.n_processed + read_num;

        mtx_n_processed.lock();
        aux_glb.n_processed += read_num; //aux maybe change ,need to be locked ro automatic operation
        mtx_n_processed.unlock();

        concurrentMap.assign(taskId, task_input);
        cerr << "task_input = " << taskId << endl;

        TaskReadInputType* task = concurrentMap.get(taskId);
        for (i=0; i<read_num; i++){
            char *seq = task->one_read[i]->query;
            int l_seq = point->seqs[i].l_seq;
            int read_id = point->seqs[i].id; //todo maybe need to use %
            mem_chain_v *chn = task->one_read[i]->v;
            int64_t **rmax = task->one_read[i]->ramx_data;
            TasktoFpga(opt, idx, seq, chn, read_id, taskId, l_seq, rmax);
        }

        free(point);// free the ret, but the ret->seq have no free
    }
}

void sendDataToFpga(void *args){

    ByteBuffer SendBuffer(PACKAGE_SIZE);
    fpga_data_input* Data;
    Byte* temp = (Byte*)malloc(PACKAGE_SIZE*sizeof(Byte)); //todo need to free
    if (temp == NULL){
        cerr << "error to malloc space";
    }
    memset(temp,0, PACKAGE_SIZE);
    int cal_num = 0;
    int flag = 1;
    int pkg_count = 0;
    while(true){
        bool isGet = blockingQueueStep1.wait_dequeue_timed(Data,SEND_DELAY);//超时等待并且出队
        if(isGet){
            int length = Data->length;
            int num_zero = ((Data->length) % 16) ? 16 - (Data->length) % 16 : 0;
            if(SendBuffer.size()+length +num_zero <PACKAGE_SIZE){
                cal_num ++ ;
                if (flag == 1) {
                    SendBuffer.putLong(0);//first 128bit save the num of cal unit per 128k
                    SendBuffer.putLong(0);
                    flag = 0;
                }
                DataTobuffer(Data, SendBuffer);
                //todo free rs qs

            }else{
                pkg_count ++ ;
                SendBuffer.getBytes(temp,SendBuffer.size());//将SendBuffer中的数据拷贝到temp
                for (int i=0 ; i<16; i++){ //save the number of cal unit to the SendBuffer
                    if (i<12) temp[i] = 0;
                    else{
                        temp[i] = cal_num >> (8*(15-i));
                    }
                }

                cal_num = 1;
                // lock the SendBuffer is also the problem
                mtx_writeToFpga.lock();
                writeDataToFPGA(temp);//push the data to FPGA
                mtx_writeToFpga.unlock();
                SendBuffer.clear();//清空SendBuffer以用来下一份数据集的发送
                SendBuffer.putLong(0);//first 128bit save the num of cal unit per 128k
                SendBuffer.putLong(0);

                DataTobuffer(Data, SendBuffer);
            }
        }else{//超时时间里面没有获取到数据,接下来判断SendBuffer中是否存在有效数据,如果有则发送,如果没有则继续等待获取数据
            if(SendBuffer.size()>0){
                SendBuffer.getBytes(temp,SendBuffer.size());//将SendBuffer中的数据拷贝到temp
                for (int i=12 ; i<16; i++){ //save the number of cal unit to the SendBuffer
                    temp[i] = cal_num >> (8*(15-i));
                }
                mtx_writeToFpga.lock();
                writeDataToFPGA(temp);
                mtx_writeToFpga.unlock();
                pkg_count ++;
                SendBuffer.clear();
                flag = 1;
                cal_num = 0;

            }else{
                continue;
            }
        }
    }
}
void writeDataToFPGA(Byte* buffer){
    while (true) {
        if (IsDsWriteable(1)) {
            WriteDsPkg(buffer, 1);
            break;
        }
    }
}
void GetDataFromFpga(void* args){

    while (true){

        Byte* GetBuffer = (Byte*)malloc(FPGA_OUT_PACKAGE_SIZE * sizeof(char));
        memset(GetBuffer, 0, FPGA_OUT_PACKAGE_SIZE);
        mtx_readFromFpga.lock();
        while (true) {//
            //sleep(1);
            if (IsUsReadable(1)) {
                ReadUsPkg(GetBuffer, 1);
                break;
            }
        }
        mtx_readFromFpga.unlock();
        // directly push the output to the Query
        int num = FPGA_OUT_PACKAGE_SIZE/32;
        for (int i=0; i< num; i++){
            fpga_sw_output* fpgaoutput = (fpga_sw_output*)malloc(sizeof(fpga_sw_output));
            memset(fpgaoutput, 0, num);
            fpgaoutput->task_id = GetBuffer[i * 32 + 31] * 256 + GetBuffer[i * 32 + 30]; //todo maybe some problem
            fpgaoutput->read_id = GetBuffer[i * 32 + 29] * 256 * 256 * 256 + GetBuffer[i * 32 + 28] * 256 * 256 + GetBuffer[i * 32 + 27] * 256 + GetBuffer[i * 32 + 26];
            fpgaoutput->chain_id = GetBuffer[i * 32 + 25] * 256 + GetBuffer[i * 32 + 24];
            fpgaoutput->seed_id = GetBuffer[i * 32 + 23] * 256 + GetBuffer[i * 32 + 22];
            fpgaoutput->flag[0] = GetBuffer[i * 32 + 21] & 0x03; //stand for whether extend
            fpgaoutput->flag[1] = GetBuffer[i * 32 + 21] & 0x0c; //stand for whether coordinate need to plus 1
            //GetBuffer[11] is the zero
            fpgaoutput->result[0].max = GetBuffer[i * 32 + 19] * 256 + GetBuffer[i * 32 + 18];//left_result
            fpgaoutput->result[0].max_j = GetBuffer[i * 32 + 17] * 256 + GetBuffer[i * 32 + 16];
            fpgaoutput->result[0].max_i = GetBuffer[i * 32 + 15];
            fpgaoutput->result[0].g_score = GetBuffer[i * 32 + 14] * 256 + GetBuffer[i * 32 + 13];
            fpgaoutput->result[0].g_score_j = GetBuffer[i * 32 + 12] * 256 + GetBuffer[i * 32 + 11]; //todo maybe the coordinate need to add 1
            fpgaoutput->result[0].g_score_i = GetBuffer[i * 32 + 10];

            fpgaoutput->result[1].max = GetBuffer[i * 32 + 9] * 256 + GetBuffer[i * 32 + 8];//left_result
            fpgaoutput->result[1].max_j = GetBuffer[i * 32 + 7] * 256 + GetBuffer[i * 32 + 6];
            fpgaoutput->result[1].max_i = GetBuffer[i * 32 + 5];
            fpgaoutput->result[1].g_score = GetBuffer[i * 32 + 4] * 256 + GetBuffer[i * 32 + 3];
            fpgaoutput->result[1].g_score_j = GetBuffer[i * 32 + 2] * 256 + GetBuffer[i * 32 + 1];
            fpgaoutput->result[1].g_score_i = GetBuffer[i * 32 + 0];

            if ((GetBuffer[i*32+31] == 0xff) && (GetBuffer[i*32+29] == 0xff) && GetBuffer[i*32+21] == 0xff) {
                //fpga_output(fpgaoutput[0]);
                free(fpgaoutput);
                continue ;
            } //ddr上会有缓存数据没有清理干净
            
            blockingQueueStep2.enqueue(fpgaoutput); //push data to queue
        }
        free(GetBuffer);
    }
}

#define LIKELY(x) __builtin_expect((x),1)
void Soft_SW(void *args){
    while(true){

        fpga_data_input* point_fpga_input;  //the input of fpga,
        blockingQueueStep1.wait_dequeue(point_fpga_input);
        int qle, tle, gtle, gscore;
        int o_del = aux_glb.opt->o_del;
        int e_del = aux_glb.opt->e_del;
        int o_ins = aux_glb.opt->o_ins;
        int e_ins = aux_glb.opt->e_ins;
        const int8_t  *mat = aux_glb.opt->mat;
        int aw[2], max_off[2];
        aw[0] = aw[1] = aux_glb.opt->w;
        int pen_clip5 = aux_glb.opt->pen_clip5;

        int tlen = point_fpga_input->sw_cal[0].r_len;
        int qlen = point_fpga_input->sw_cal[0].q_len;
        int zdrop = aux_glb.opt->zdrop;
        uint8_t *query = point_fpga_input->sw_cal[0].qs;
        uint8_t *target = point_fpga_input->sw_cal[0].rs;
        int score_init = point_fpga_input->sw_cal[0].score_init;
        fpga_sw_output * sw_output;
        sw_output = (fpga_sw_output*)malloc(sizeof(fpga_sw_output)); //free
        sw_output->flag[1] = 0x00; //no plus 1
        if (qlen){
            sw_output->flag[0] = 0x2;
            sw_output->result[0].max = ksw_extend2(qlen, query, tlen, target, 5, mat, o_del, e_del, o_ins, e_ins, aw[0],
                                                   pen_clip5, zdrop, score_init, &qle, &tle, &gtle, &gscore, &max_off[0]);
            sw_output->result[0].max_i = qle;
            sw_output->result[0].max_j = tle;
            sw_output->result[0].g_score = gscore;
            sw_output->result[0].g_score_i = qlen; // todo maybe need to substrat 1
            sw_output->result[0].g_score_j = gtle;
        }else {
            sw_output->flag[0] = 0x0; // zero stand for no extend
            sw_output->result[0].max = score_init;
            sw_output->result[0].max_i = 0;
            sw_output->result[0].max_j = 0;
            sw_output->result[0].g_score = 0;
            sw_output->result[0].g_score_i = 0;
            sw_output->result[0].g_score_j = 0;
        }

        tlen = point_fpga_input->sw_cal[1].r_len;  //right extension
        qlen = point_fpga_input->sw_cal[1].q_len;
        query = point_fpga_input->sw_cal[1].qs;
        target = point_fpga_input->sw_cal[1].rs;
        score_init = sw_output->result[0].max;
        aw[0] = score_init;
        if (qlen){
            sw_output->flag[0] = sw_output->flag[0] | 0x1;
            sw_output->result[1].max = ksw_extend2(qlen, query, tlen, target, 5, mat, o_del, e_del, o_ins, e_ins, aw[0],
                                                   pen_clip5, zdrop, score_init, &qle, &tle, &gtle, &gscore, &max_off[0]);
            sw_output->result[1].max_i = qle;
            sw_output->result[1].max_j = tle;
            sw_output->result[1].g_score = gscore;
            sw_output->result[1].g_score_i = qlen;
            sw_output->result[1].g_score_j = gtle;

        }else{
            sw_output->flag[0] = sw_output->flag[0] | 0x0;
            sw_output->flag[1] = 0;
            sw_output->result[1].max = sw_output->result[0].max;
            sw_output->result[1].max_i = 0;
            sw_output->result[1].max_j = 0;
            sw_output->result[1].g_score = 0;
            sw_output->result[1].g_score_i = 0;
            sw_output->result[1].g_score_j = 0;
        }

        sw_output->task_id = point_fpga_input->task_id;
        sw_output->read_id = point_fpga_input->id.read_id;
        sw_output->chain_id = point_fpga_input->id.chain_id;
        sw_output->seed_id = point_fpga_input->id.seed_id;
        //==========test======
        // fpga_output(sw_output[0]);
        //===================
        blockingQueueStep2.enqueue(sw_output); //free in the process2

        // free the sw fpga input queue
        free(point_fpga_input->sw_cal[0].qs);
        free(point_fpga_input->sw_cal[0].rs);
        free(point_fpga_input->sw_cal[1].qs);
        free(point_fpga_input->sw_cal[1].rs);
        free(point_fpga_input);
    }
}

static void *process2(void *data)
{// 第二步：调用mem_process_seqs对seqs进行处理，直到写入SAM文件
    fpga_sw_output* point;

    while(true){

        blockingQueueStep2.wait_dequeue(point);
        int task_id = point->task_id;
        int read_id = point->read_id;
        int chain_id = point->chain_id;
        int seed_id = point->seed_id;
        TaskReadInputType* task = concurrentMap.get(task_id);

        task->one_read[read_id]->output_buf[chain_id][seed_id].result[0] = point->result[0]; //put the result in the buf
        task->one_read[read_id]->output_buf[chain_id][seed_id].result[1] = point->result[1];
        task->one_read[read_id]->output_buf[chain_id][seed_id].flag[0] = point->flag[0];
        task->one_read[read_id]->output_buf[chain_id][seed_id].flag[1] = point->flag[1];

        mtx_erase_task.lock();
        task->count ++ ;//todo maybe a atom operation
        if (task->count == task->num_seed){
            concurrentMap.erase(task_id);
            blockingQueueStep3.enqueue(task);
        }
        mtx_erase_task.unlock();

        free(point);
    }
}

static void process3(void *data)
{// 第二步：调用mem_process_seqs对seqs进行处理，直到写入SAM文件
    TaskReadInputType* point;
    int task_count_test = 0;

    while(true){

        blockingQueueStep3.wait_dequeue(point);
        //==============================
        task_count_test ++;  //the value is 15
        //cerr << task_count_test << endl ;
        int read_num = point->read_num;

        mem_pestat_t pes[4]; // todo maybe some problem
        const mem_pestat_t *pes0 = aux_glb.pes0;
        const mem_opt_t *opt = aux_glb.opt;
        const bntseq_t *bns = aux_glb.idx->bns;
        const uint8_t *pac = aux_glb.idx->pac;

        mem_alnreg_v *regs = (mem_alnreg_v*)malloc(read_num*sizeof(mem_alnreg_v));
        for(int i=0; i<read_num; i++){

            kv_init(regs[i]);
            mem_chain_v *chn = point->one_read[i]->v;
            fpga_sw_output **fpga_out = point->one_read[i]->output_buf;
            int64_t **rmax_data = point->one_read[i]->ramx_data;
            char* seq = point->one_read[i]->query;
            int l_seq = point->one_read[i]->l_seq;
            mem_align1_core_suf(aux_glb.opt, chn, l_seq, rmax_data, fpga_out, &regs[i]); //free the chan and seeds, rmax_data, fpga_output
            //===============================TEST===================
            //fpga_regs_output(regs[i]);
            //================================END===================
            regs[i].n = mem_sort_dedup_patch(opt, bns, pac, (uint8_t*)seq, regs[i].n, regs[i].a);


            unsigned int j = 0;
            if (bwa_verbose >= 4) {
                err_printf("* %ld chains remain after removing duplicated chains\n", regs[i].n);
                for (j = 0; j < regs[i].n; ++j) {
                    mem_alnreg_t *p = &regs[i].a[i];
                    printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
                }
            }
            for (j = 0; j < regs[i].n; ++j) {
                mem_alnreg_t *p = &regs[i].a[j];
                if (p->rid >= 0 && aux_glb.idx->bns->anns[p->rid].is_alt)// 如果对应的染色体的is_alt值为1，则将lan的is_alt也变成1
                    p->is_alt = 1;
            }

            free(point->one_read[i]->v);
            free(point->one_read[i]);

        }
        free(point->one_read);

        //todo
        if (aux_glb.opt->flag&MEM_F_PE) { // infer insert sizes if not provided
            if (pes0) memcpy(pes, pes0, 4 * sizeof(mem_pestat_t)); // if pes0 != NULL, set the insert-size distribution as pes0
            else mem_pestat(opt, bns->l_pac, read_num, regs, pes); // otherwise, infer the insert size distribution from data，从这一批reads的数据中得出每个direction下的infer_size的情况，保存在pes中
        }

        int64_t n_processed = point->n_processed - point->read_num;
        bseq1_t *seqs = point->seqs;//todo in the process2 the task is free ,the seqs may be free too
        read_num = (opt->flag&MEM_F_PE)? read_num>>1 : read_num;
        for (int i=0; i<read_num; i++){ //todo maybe it can can use multi thread to deal with// i=57628, taskCounttest = 29 报错
            if (!(opt->flag&MEM_F_PE)) {// 跳过，直接执行else MEM_F_PE = Ox2,如果输入只有一个文件，那么就执行if
                if (bwa_verbose >= 4) printf("=====> Finalizing read '%s' <=====\n", seqs[i].name);
                mem_mark_primary_se(opt, regs[i].n, regs[i].a, aux_glb.n_processed + i);
                mem_reg2sam(opt, bns, pac, &seqs[i], &regs[i], 0, 0);
                free(regs[i].a);
            } else {
                if (bwa_verbose >= 4) printf("=====> Finalizing read pair '%s' <=====\n", seqs[i<<1|0].name);
                mem_sam_pe(opt, bns, pac, pes, (n_processed>>1) + i, &seqs[i<<1], &regs[i<<1]);//移位的目的是为了使指针指向当前位置
                free(regs[i<<1|0].a); free(regs[i<<1|1].a);//todo have problem
            }

        }

        free(regs);
        //todo step3 relase the memory
        for (int i = 0; i < point->read_num; ++i) {
            if (seqs[i].sam) err_fputs(seqs[i].sam, stdout);
            free(seqs[i].name); free(seqs[i].comment);
            free(seqs[i].seq); free(seqs[i].qual); free(seqs[i].sam);
        }
        free(seqs);
        free(point);
        cerr << "task_count_test:" << task_count_test << endl ;
        mapTaskCount++;
    }
}

void threadpoolInit(ThreadPoolParameter threadPoolParameter){
    ThreadPool step0Thread(threadPoolParameter.Step0ThreadNum);
    ktp_aux_t* process0Args = &aux_glb;
    unsigned int i = 0;
    for (i = 0; i < threadPoolParameter.Step0ThreadNum; ++i) {
        step0Thread.enqueue(threadPoolParameter.Step0Func,(void*)process0Args);
    }

    ThreadPool step1Thread(threadPoolParameter.Step1ThreadNum);
    for (i = 0; i < threadPoolParameter.Step1ThreadNum; ++i) {
        step1Thread.enqueue(threadPoolParameter.Step1Func,(void*)NULL);//no arguments
    }

    ThreadPool sendToFpga(threadPoolParameter.SendToFpgaNum);
    for (i = 0; i < threadPoolParameter.SendToFpgaNum; ++i) {
        sendToFpga.enqueue(threadPoolParameter.SendToFpgaFunc, (void*)NULL);
    }

    ThreadPool GetFromFpga(threadPoolParameter.GetFromFpgaNum);
    for (i = 0; i < threadPoolParameter.GetFromFpgaNum; ++i) {
        GetFromFpga.enqueue(threadPoolParameter.GetFromFpgaFunc, (void*)NULL);
    }

    ThreadPool step2Thread(threadPoolParameter.Step2ThreadNum);
    for (i = 0; i < threadPoolParameter.Step2ThreadNum; ++i) {
        step2Thread.enqueue(threadPoolParameter.Step2Func,(void*)NULL);
    }

    ThreadPool step3Thread(threadPoolParameter.Step3ThreadNum);
    for (i = 0; i < threadPoolParameter.Step3ThreadNum; ++i) {
        step3Thread.enqueue(threadPoolParameter.Step3Func,(void*)NULL);
    }
}


void bwa_mem_fpga(void){
    
    junction::QSBR::Context  context = junction::DefaultQSBR.createContext();
    junction::DefaultQSBR.update(context);
    junction::DefaultQSBR.destroyContext(context);//todo have to test
    ThreadPoolParameter threadPoolParameter;
    if (Fpga_flag == 0){
        threadPoolParameter.Step0ThreadNum = 2;
        threadPoolParameter.Step1ThreadNum = 2;
        threadPoolParameter.SendToFpgaNum = 0;
        threadPoolParameter.GetFromFpgaNum = 2;
        threadPoolParameter.Step2ThreadNum = 2;
        threadPoolParameter.Step3ThreadNum = 2;
        threadPoolParameter.Step0Func = process0;
        threadPoolParameter.Step1Func = process1;
        threadPoolParameter.SendToFpgaFunc = 0;
        threadPoolParameter.GetFromFpgaFunc = Soft_SW;
        threadPoolParameter.Step2Func = process2;
        threadPoolParameter.Step3Func = process3;
        threadpoolInit(threadPoolParameter);

    }
    else{
        //================FPGA_TEST=================
        pthread_t tid;
        pthread_create(&tid, NULL, printDebugInfo, NULL);
        //====================================================
        int open_statue = OpenQueueCtrl();
        if (!open_statue){
            cerr<<"\nOpen xdma device failed!\n";
            exit(-1);
        }
        cerr<<"\nOpen xdma device success\n";
        threadPoolParameter.Step0ThreadNum = 4;
        threadPoolParameter.Step1ThreadNum = 1;
        threadPoolParameter.SendToFpgaNum = 1;
        threadPoolParameter.GetFromFpgaNum = 1;
        threadPoolParameter.Step2ThreadNum = 1;
        threadPoolParameter.Step3ThreadNum = 1;
        threadPoolParameter.Step0Func = process0;
        threadPoolParameter.Step1Func = process1;
        threadPoolParameter.SendToFpgaFunc = sendDataToFpga;
        threadPoolParameter.GetFromFpgaFunc = GetDataFromFpga;
        threadPoolParameter.Step2Func = process2;
        threadPoolParameter.Step3Func = process3;
        threadpoolInit(threadPoolParameter);
    }
    CloseQueueCtrl();
    
}