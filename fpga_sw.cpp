//
// Created by jie on 10/23/18.
//

#include "fpga_sw.h"
#include "utils.h"
#include "kvec.h"


extern "C"{
    mem_chain_v mem_chain(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, int len, const uint8_t *seq, void *buf);
    extern int mem_chain_flt(const mem_opt_t *opt, int n_chn, mem_chain_t *a);
    extern void mem_print_chain(const bntseq_t *bns, mem_chain_v *chn);
}

BlockingQueueStep1 blockingQueueStep1(BLOCKQUEUE_STEP1_SIZE); //

void Fpga_input(fpga_data_input* input){
    FILE *file;
    file = fopen("/home/jie/Documents/wes_input/fpga_sw_input.txt", "wb");
    if (!file) perror("/home/jie/Documents/wes_input/fpga_sw_input.txt"), exit(1);

    putc(input->id.chain_id, file);
    putc(input->id.read_id, file);
    putc(input->id.seed_id, file);
    putc(input->task_id, file);

    {
        putc(input->sw_cal[0].q_len,file);
        for (int i=0; i<input->sw_cal[0].q_len; i++) putc(input->sw_cal[0].qs[i], file);
        putc(input->sw_cal[0].r_len, file);
        for (int i=0; i<input->sw_cal[0].r_len; i++) putc(input->sw_cal[0].rs[i], file);
        putc(input->sw_cal[0].score_init, file);
    }

    {
        putc(input->sw_cal[1].q_len,file);
        for (int i=0; i<input->sw_cal[1].q_len; i++) putc(input->sw_cal[1].qs[i], file);
        putc(input->sw_cal[1].r_len, file);
        for (int i=0; i<input->sw_cal[1].r_len; i++) putc(input->sw_cal[1].rs[i], file);
        putc(input->sw_cal[1].score_init, file);
    }

    fclose(file);

}
extern void fpga_chain_output(mem_chain_v *chain);
void mem_align1_core_pre(const mem_opt_t *opt, const bwaidx_t *idx, int l_seq, char* seq, TaskReadInputType* task_input, int temp, int *num_seed, void* buf)
{

    int i, j;
    mem_chain_v chn;
    const bwt_t *bwt = idx->bwt;
    const bntseq_t *bns = idx->bns;
    for (i=0; i<l_seq; i++)
        seq[i] = seq[i] < 4 ? seq[i] : nst_nt4_table[(int)seq[i]];

    chn = mem_chain(opt, bwt, bns, l_seq, (uint8_t*)seq, buf);
    //========================TEST==========//
    //fpga_chain_output(&chn);
    //=======================END===========//

    chn.n = mem_chain_flt(opt, chn.n, chn.a);
    if (bwa_verbose >= 4) mem_print_chain(bns, &chn);
    int64_t **rmax;
    for (j = 0; j < chn.n; ++j){
        mem_chain_t *ctemp ;
        ctemp = &(chn.a[j]);
        int k = 0;
        for (k = ctemp->n-1; k>=0; --k){
            (*num_seed) ++ ;
        }
    }
    rmax = (int64_t**)malloc((chn.n)*sizeof(int64_t*));

    for (j = 0; j < chn.n; ++j) {
        mem_chain_t *c = &chn.a[j];
        if (bwa_verbose >= 4) err_printf("* ---> Processing chain(%d) <---\n", i);
        int64_t l_pac = bns->l_pac, max = 0;
        rmax[j] = (int64_t*)malloc(2*sizeof(int64_t));

        // get the max possible span
        rmax[j][0] = l_pac<<1; rmax[j][1] = 0;
        for (i = 0; i < c->n; ++i) {
            int64_t b, e;
            const mem_seed_t *t = &c->seeds[i];
            b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
            e = t->rbeg + t->len + ((l_seq - t->qbeg - t->len) + cal_max_gap(opt, l_seq - t->qbeg - t->len));
            rmax[j][0] = rmax[j][0] < b? rmax[j][0] : b;// reference the forward position
            rmax[j][1] = rmax[j][1] > e? rmax[j][1] : e;// reference the backward position
            if (t->len > max) max = t->len;
        }
        rmax[j][0] = rmax[j][0] > 0? rmax[j][0] : 0;
        rmax[j][1] = rmax[j][1] < l_pac<<1? rmax[j][1] : l_pac<<1;
        if (rmax[j][0] < l_pac && l_pac < rmax[j][1]) { // crossing the forward-reverse boundary; then choose one side
            if (c->seeds[0].rbeg < l_pac) rmax[j][1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
            else rmax[j][0] = l_pac;
        }
    }

    fpga_sw_output** output_buf = (fpga_sw_output**)malloc(chn.n * sizeof(fpga_sw_output*)) ;// chain space
    for (int cj=0; cj<chn.n; cj++){
        output_buf[cj] = (fpga_sw_output*)malloc(chn.a[cj].n * sizeof(fpga_sw_output)); //seed space
    }

    OneReadInputType* Data ;
    Data = (OneReadInputType*)malloc(sizeof(OneReadInputType));
    Data->l_seq = l_seq;
    Data->v = (mem_chain_v*)malloc(sizeof(mem_chain_v)); 
    *(Data->v) = chn;
    Data->ramx_data = rmax;
    Data->query = seq;
    Data->output_buf = output_buf;
    task_input->one_read[temp] = Data;
}

extern int Fpga_flag;
void mem_align1_core_suf(const mem_opt_t *opt, const mem_chain_v *chn, int l_query, int64_t **rmax_data, fpga_sw_output** fpga_output, mem_alnreg_v *av)
{
    mem_chain_t *c;
    int64_t rmax[2];
    int aw[2], max_off[2];

    for (int cj=0; cj<chn->n; cj++)
    {
        c = &(chn->a[cj]);
        rmax[0] = rmax_data[cj][0];
        rmax[1] = rmax_data[cj][1];

        uint64_t *srt;
        srt = (uint64_t*)malloc(c->n * 8);
        for (int j = 0; j < c->n; ++j)
            srt[j] = (uint64_t)c->seeds[j].score<<32 | j;
        ks_introsort_64(c->n, srt);
        int i, k = 0;
        const mem_seed_t *s;
        for (k = (c->n)-1; k>=0; --k){
            mem_alnreg_t *a ;
            s = &c->seeds[(uint32_t)srt[k]];
            for (i = 0; i < av->n; ++i) { // test whether extension has been made befor
                mem_alnreg_t *p = &av->a[i];//!这里的av是针对虽有chain的一个aln结果，即针对这个read的所有aln结果
                int64_t rd;
                int qd, w, max_gap;
                if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb || s->qbeg + s->len > p->qe) continue; // not fully contained，首先检查在reference中的覆盖情况，再检查在read中的匹配情况，如果没有被全部覆盖，则要继续处理
                if (s->len - p->seedlen0 > .1 * l_query) continue; // this seed may give a better alignment，如果这个seed的长度更长，那么要继续处理
                // qd: distance ahead of the seed on query; rd: on reference，判断seed在允许的gap范围内能否被已有的aln覆盖
                qd = s->qbeg - p->qb;/*当前seed在read中开始位置与已有aln在read中开始位置的距离*/ rd = s->rbeg - p->rb;/*当前seed在reference中开始位置与已有aln在reference中开始位置的距离*/
                max_gap = cal_max_gap(opt, qd < rd? qd : rd); // the maximal gap allowed in regions ahead of the seed
                w = max_gap < p->w? max_gap : p->w; // bounded by the band width，p->w是之前aln用的band，以p->w为上界，取最小值，表示最seed之前能延伸的最大长度
                if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
                // similar to the previous four lines, but this time we look at the region behind
                qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
                max_gap = cal_max_gap(opt, qd < rd? qd : rd);
                w = max_gap < p->w? max_gap : p->w;
                if (qd - rd < w && rd - qd < w) break;// break表示这个seed已经在之前某个aln的附近
            }// 进一步判断，这个for并不是从上一个for循环跳出的位置开始遍历，而是又重新遍历
            if (i < av->n) { // the seed is (almost) contained in an existing alignment; further testing is needed to confirm it is not leading to a different aln
                if (bwa_verbose >= 4)
                    printf("** Seed(%d) [%ld;%ld,%ld] is almost contained in an existing alignment [%d,%d) <=> [%ld,%ld)\n",
                           k, (long)s->len, (long)s->qbeg, (long)s->rbeg, av->a[i].qb, av->a[i].qe, (long)av->a[i].rb, (long)av->a[i].re);
                for (i = k + 1; i < c->n; ++i) { // check overlapping seeds in the same chain，遍历这个seed之后的所有seed
                    const mem_seed_t *t;
                    if (srt[i] == 0) continue;// str[k]==0表示这个seed没有进行extension，即不需要处理的seed，跳过
                    t = &c->seeds[(uint32_t)srt[i]];// s需要比较的seed，t是余下的所有seed
                    if (t->len < s->len * .95) continue; // only check overlapping if t is long enough;只有t的长度与待检测seed长度相比足够长时才比较 TODO: more efficient by early stopping
                    if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;// 1.s的起始位置位于t之前，2.s与t重叠部分的长度超过s长度的一半，3.s和t在read中的起始位置差和在reference中的起始位置差不同
                    if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;// 1.s的结束位置位于t之后，2.s与t重叠部分的长度超过s长度的一半，3.s和t在read中的起始位置差和在reference中的起始位置差不同
                }// break出来就代表有overlap的seed，有可能导致一个新的aln，需要继续往下面的extension执行
                if (i == c->n) { // no overlapping seeds; then skip extension，如果没有overlap的seeds，那么跳过这个seed
                    srt[k] = 0; // mark that seed extension has not been performed
                    continue;// 这个continue是针对最外层的for循环
                }
                if (bwa_verbose >= 4)
                    printf("** Seed(%d) might lead to a different alignment even though it is contained. Extension will be performed.\n", k);
            }

            a = kv_pushp(mem_alnreg_t, *av);// 在av中保存当前align的结果
            memset(a, 0, sizeof(mem_alnreg_t));
            a->w = aw[0] = aw[1] = opt->w;
            a->score = a->truesc = -1;
            a->rid = c->rid;
            // extend left
            int gscore, qle, tle, gtle ;
            gscore = fpga_output[cj][k].result[0].g_score;
            if ((fpga_output[cj][k].flag[1] & 0x08) == 0 ){
                qle = fpga_output[cj][k].result[0].max_i;
                tle = fpga_output[cj][k].result[0].max_j;
            }
            else{
                qle = fpga_output[cj][k].result[0].max_i + 1;
                tle = fpga_output[cj][k].result[0].max_j + 1;
            }
            if (Fpga_flag == 0)
                gtle = fpga_output[cj][k].result[0].g_score_j;
            else
                gtle =fpga_output[cj][k].result[0].g_score_j + 1;
            
            max_off[0] = abs(qle - tle);
            if ((fpga_output[cj][k].flag[0] & 0x2) == 0){ //flga[0] stand for whether extend
                a->score = s->len;
                a->truesc = s->len;
                a->qb = 0;
                a->rb = s->rbeg;
            }else{
                a->score = fpga_output[cj][k].result->max;
                aw[0] = opt->w << 0;
                if (max_off[0] >= (aw[0]>>1) + (aw[0]>>2)) aw[0] = opt->w << 1;
                if (gscore <= 0 || gscore <= a->score - opt->pen_clip5){
                    a->qb = s->qbeg - qle, a->rb = s->rbeg - tle;
                    a->truesc = a->score;
                }else{
                    a->qb = 0, a->rb = s->rbeg - gtle;
                    a->truesc = gscore;
                }
            }

            //extern right
            int qe = s->qbeg + s->len;
            int re = s->rbeg + s->len - rmax[0];
            if ((fpga_output[cj][k].flag[1] & 0x04) == 0){
                qle = fpga_output[cj][k].result[1].max_i;
                tle = fpga_output[cj][k].result[1].max_j;
            }
            else{
                qle = fpga_output[cj][k].result[1].max_i + 1;
                tle = fpga_output[cj][k].result[1].max_j + 1;
            }
            if (Fpga_flag == 0)
                gtle = fpga_output[cj][k].result[1].g_score_j;
            else
                gtle =fpga_output[cj][k].result[1].g_score_j + 1;

            if ((fpga_output[cj][k].flag[0] & 0x01) == 0){ //如果不做拓展
                a->qe = l_query;
                a->re = s->rbeg + s->len;
            }else{
                int sc0 = a->score;
                a->score = fpga_output[cj][k].result[1].max;
                aw[1] = opt->w << 0;

                gscore = fpga_output[cj][k].result[1].g_score;

                max_off[1] = abs(qle - tle);
                if (max_off[1] >= (aw[1]>>1) + (aw[1]>>2)) aw[1] = opt->w << 1;

                if (gscore <= 0 || gscore <= a->score - opt->pen_clip3){
                    a->qe = qe + qle, a->re = rmax[0] + re + tle;
                    a->truesc += a->score - sc0;
                }else{
                    a->qe = l_query, a->re = rmax[0] + re + gtle;
                    a->truesc += gscore - sc0;
                }
            }

            int m = 0;
            for (m=0, a->seedcov=0; m<c->n; m++){
                const mem_seed_t *t = &(c->seeds[m]);
                if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
                    a->seedcov += t->len;
            }
            a->w = aw[0] > aw[1]? aw[0] : aw[1];
            a->seedlen0 = s->len;
            a->frac_rep = c->frac_rep;
        }
        free(srt);
        free(rmax_data[cj]);
        free(chn->a[cj].seeds);
        free(fpga_output[cj]);

    }
    free(rmax_data);
    free(chn->a);
    free(fpga_output);
}

void DataTobuffer(fpga_data_input* Data, ByteBuffer &buffer){
    int i = 0;
    buffer.putLong(0);//second 128bit save the length of next cal unit, its unit is 128bit but the length is 8bit
    int num_zero = ((Data->length) % 16) ? 16 - (Data->length) % 16 : 0;
    buffer.putInt(0);
    buffer.putInt_Little((Data->length -1)/16 + 1);

    buffer.putShort_Little(Data->task_id);//use 2byte to save the id, taskID need to be recycled
    buffer.putInt_Little(Data->id.read_id);// use 4byte
    buffer.putShort_Little(Data->id.chain_id); // use 2 byte
    buffer.putShort_Little(Data->id.seed_id); // use 2 byte
    buffer.putShort_Little(0); // total is 12byte fill the zero is to fit the FPGA

    buffer.putShort_Little(Data->sw_cal[0].q_len);//use 2+2Byte sve read_len r_len
    buffer.putShort_Little(Data->sw_cal[0].r_len);

    buffer.putInt_Little(Data->sw_cal[0].score_init); //use 4byte save socre_init

    uint8_t *rs = Data->sw_cal[0].rs;
    uint8_t *qs = Data->sw_cal[0].qs;
    int r_len = Data->sw_cal[0].r_len;
    int q_len = Data->sw_cal[0].q_len;
    int len_left = 0;

    if (q_len != 0) {
        for (i = 0; i < ((q_len - 1) / 8) + 1; i++) {//4bit save one bp,save qs
            len_left = q_len - 8 * i;
            int q_query = 0;
            int j = 0;
            if (len_left <= 8) {
                for (j = 0; j < len_left; j++) {
                    int shift = 28 - 4 * j;
                    q_query = q_query | (qs[i*8+j] << shift);
                }
                buffer.putInt_Little(q_query);
            } else {
                for (j = 0; j < 8; j++) {
                    int shift = 28 - 4 * j;
                    q_query = q_query | (qs[i*8+j] << shift);
                }
                buffer.putInt_Little(q_query);
            }
        }
        for (i = 0; i < ((r_len - 1) / 8) + 1; i++) {//4bit save one bp save rs
            len_left = r_len - 8 * i;
            int q_query = 0;
            int j = 0;
            if (len_left <= 8) {
                for (j = 0; j < len_left; j++) {
                    int shift = 28 - 4 * j;
                    q_query = q_query | (rs[i*8+j] << shift);
                }
                buffer.putInt_Little(q_query);
            } else {
                for (j = 0; j < 8; j++) {
                    int shift = 28 - 4 * j;
                    q_query = q_query | (rs[i*8+j] << shift);
                }
                buffer.putInt_Little(q_query); //todo maybe wrong
            }
        }
    }

    buffer.putShort_Little(Data->sw_cal[1].q_len);//use 2+2Byte sve q_len r_len
    buffer.putShort_Little(Data->sw_cal[1].r_len);

    buffer.putInt_Little(Data->sw_cal[1].score_init); //use 4byte save socre_init

    rs = Data->sw_cal[1].rs;
    qs = Data->sw_cal[1].qs;
    r_len = Data->sw_cal[1].r_len;
    q_len = Data->sw_cal[1].q_len;
    if (q_len != 0){
        for(i=0; i< ((q_len-1)/8) +1; i++)  {//4bit save one bp,save qs
            len_left = q_len - 8*i;
            int q_query = 0;
            int j=0;
            if (len_left <= 8){
                for (j=0; j<len_left; j++){
                    int shift = 28 - 4*j;
                    q_query = q_query | (qs[i*8+j] << shift);
                }
                buffer.putInt_Little(q_query);
            }
            else{
                for (j=0; j<8; j++){
                    int shift = 28 - 4*j;
                    q_query = q_query | (qs[i*8+j] << shift);
                }
                buffer.putInt_Little(q_query);
            }
        }

        for(i=0; i< ((r_len-1)/8) +1; i++)  {//4bit save one bp save rs
            len_left = r_len - 8*i;
            int q_query = 0;
            int j=0;
            if (len_left <= 8){
                for (j=0; j<len_left; j++){
                    int shift = 28 - 4*j;
                    q_query = q_query | (rs[i*8+j] << shift);
                }
                buffer.putInt_Little(q_query);
            }
            else{
                for (j=0; j<8; j++){
                    int shift = 28 - 4*j;
                    q_query = q_query | (rs[i*8+j] << shift);
                }
                buffer.putInt_Little(q_query);
            }
        }

    }
    for (i=0; i<num_zero; i++) buffer.putChar(0);

    free(Data->sw_cal[0].qs);
    free(Data->sw_cal[0].rs);
    free(Data->sw_cal[1].qs);
    free(Data->sw_cal[1].rs);
    free(Data);

}

/*void writeDataToFPGA(Byte* buffer){
    while (true) {
        if (IsDsWriteable(1)) {
            WriteDsPkg(buffer, 1);
            break;
        }
    }
}*/

void TasktoFpga(const mem_opt_t *opt, const bwaidx_t *idx, char *seq, mem_chain_v *chn, int read_id, int task_id, int l_seq, int64_t **rmax ){

    bntseq_t *bns = idx->bns;
    uint8_t *pac = idx->pac;
    unsigned int j;
    for (j = 0; j < chn->n; ++j) {// 到这里为止，chain是按照weight大小从大到小排序的
        mem_chain_t *c = &(chn->a[j]);//对chain每个chain遍历处理，在调用函数中又对chain中的每一个seed处理
        if (bwa_verbose >= 4) err_printf("* ---> Processing chain(%d) <---\n", j);

        int i, k, rid; // aw: actual bandwidth used in extension
        int64_t  tmp;
        mem_seed_t *s;
        uint8_t *rseq = 0;
        uint64_t *srt;
        rseq = bns_fetch_seq(bns, pac, &rmax[j][0], c->seeds[0].rbeg, &rmax[j][1], &rid);
        assert(c->rid == rid);

        srt = (uint64_t*)malloc(c->n * 8);
        for (i = 0; i < c->n; ++i)
            srt[i] = (uint64_t)c->seeds[i].score<<32 | i;// score的初始值等于seed的长度，srt记录了seed的score和编号
        ks_introsort_64(c->n, srt);// 将seed按照长度从小到大排序。之前chain已经按照长度排序过了。这里的排序是对srt排序，并没有对chain中的seeds真正排序，后面以srt为index遍历seed

        for (k = c->n-1; k>=0; --k){
            s = &c->seeds[(uint32_t)srt[k]];

            id_struct id;
            id.read_id = read_id;
            id.chain_id = j;
            id.seed_id = k;
            //向左拓展

            fpga_data_input* fpga_input;
            fpga_input = (fpga_data_input*)malloc(sizeof(fpga_data_input));
            fpga_input->id = id ;
            fpga_input->sw_cal[0].q_len = s->qbeg;
            tmp = s->rbeg - rmax[j][0];// seed在reference中的起始位置与chain可能起始位置的差值
            fpga_input->sw_cal[0].r_len = (fpga_input->sw_cal[0].q_len) ? tmp : 0;
            fpga_input->sw_cal[0].q_len = (fpga_input->sw_cal[0].r_len) ? s->qbeg : 0;

            fpga_input->sw_cal[0].qs = (uint8_t*)malloc(fpga_input->sw_cal[0].q_len);//todo free after the sw output calculate
            for (i = 0; i < fpga_input->sw_cal[0].q_len; ++i) fpga_input->sw_cal[0].qs[i] = seq[s->qbeg - 1 - i];
            fpga_input->sw_cal[0].rs = (uint8_t*)malloc(tmp); //todo free after the sw output calculate
            for (i = 0; i < fpga_input->sw_cal[0].r_len; ++i) fpga_input->sw_cal[0].rs[i] = rseq[tmp - 1 - i];// rs存储相对于reference左边的部分，并且是从后往前读取的
            fpga_input->sw_cal[0].score_init = s->len * opt->a;

            //向右拓展
            int qe, re;
            qe= s->qbeg + s->len;
            fpga_input->sw_cal[1].q_len = l_seq - qe;
            re = s->rbeg + s->len -rmax[j][0];
            fpga_input->sw_cal[1].r_len = (fpga_input->sw_cal[1].q_len) ? rmax[j][1] - rmax[j][0] -re : 0;
            fpga_input->sw_cal[1].q_len = (fpga_input->sw_cal[1].r_len) ? fpga_input->sw_cal[1].q_len : 0;

            fpga_input->sw_cal[1].qs = (uint8_t *)malloc(fpga_input->sw_cal[1].q_len);
            for (i=0; i<fpga_input->sw_cal[1].q_len; i++) fpga_input->sw_cal[1].qs[i] = seq[qe+i];//
            fpga_input->sw_cal[1].rs = (uint8_t *)malloc(fpga_input->sw_cal[1].r_len );
            for (i=re; i<fpga_input->sw_cal[1].r_len+re; i++) {
                fpga_input->sw_cal[1].rs[i-re] = rseq[i];
            }

            fpga_input->sw_cal[1].score_init = 0; //这个值是没有用的，应该根据FPGA计算得到
            fpga_input->task_id = task_id;
            //lenthe is equal to how many byte
            if (fpga_input->sw_cal[0].q_len == 0){
                fpga_input->length = 12+8; //head+read_len+hap_len+socre_init task_id 2+ read_id 4+ chain_id 2+ seed_id 2+ 16'b0
            }else{
                fpga_input->length = 12+8+((fpga_input->sw_cal[0].q_len -1)/8 + 1)*4+((fpga_input->sw_cal[0].r_len -1)/8 + 1)*4;
            }
            if (fpga_input->sw_cal[1].q_len == 0){
                fpga_input->length = fpga_input->length + 8;
            }else{
                fpga_input->length = fpga_input->length + 8+((fpga_input->sw_cal[1].q_len -1)/8 +1)*4+((fpga_input->sw_cal[1].r_len -1)/8 + 1)*4;
            }

            blockingQueueStep1.enqueue(fpga_input);

        }
        if(srt) free(srt);
        free(rseq);

    }
}



