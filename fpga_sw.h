//
// Created by jie on 3/9/19.
//

#ifndef BWA_FPGA_FPGA_SW_H
#define BWA_FPGA_FPGA_SW_H

#include <stdint.h>
#include <stddef.h>
#include "bwt.h"
#include "bwa.h"
#include "blockingconcurrentqueue.h"
#include "ByteBuffer.hpp"
#include "bwamem.h"
#include "queue_ctrl.h"
#define BLOCKQUEUE_STEP1_SIZE 1000

typedef struct {
    int64_t rbeg;// beg = begin
    int32_t qbeg, len;
    int score;
}mem_seed_t; // unaligned memory

typedef struct {
    int n, m, first, rid;
    uint32_t w:29, kept:2, is_alt:1;
    float frac_rep;
    int64_t pos;
    mem_seed_t *seeds;
}mem_chain_t;

typedef struct {
    size_t n, m;
    mem_chain_t *a;
}mem_chain_v;

typedef struct {
    bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;



//==================FPGA_INPUT==================//
typedef struct {
    int read_id;
    int16_t chain_id;
    int16_t seed_id;
}id_struct;

typedef struct {
    int q_len;
    int r_len;
    uint8_t *qs;
    uint8_t *rs;
    int score_init;
}sw_cal_struct;

typedef struct {
    int16_t task_id;
    id_struct id;
    sw_cal_struct sw_cal[2];
    int length;//结构体占用的字节数
} fpga_data_input;

//===================E N D =======================//


//======================FPGA OUTPUT========================//

typedef struct {
    int16_t max;
    int8_t max_i;
    int16_t max_j;

    int16_t g_score;
    int8_t  g_score_i;
    int16_t g_score_j;
}fpga_score;

typedef struct{
    int16_t task_id;
    int read_id;
    int16_t chain_id;
    int16_t seed_id;
    int8_t flag[2];//用于表示是否做了拓展，如果有拓展则为1,如果没有则为0,flag[0]表示左边，flag[1]表示右边
    fpga_score result[2];// 0 is the left result ,1 is the right
}fpga_sw_output;


//=========================代码上下部分传输的数据========//

struct OneReadInputType{
    int l_seq;
    mem_chain_v* v;
    int64_t ** ramx_data;
    fpga_sw_output** output_buf;
    char *query;
};


struct TaskReadInputType{
    int read_num;
    int16_t task_id;
    int num_seed;
    int count;
    int64_t n_processed;
    bseq1_t *seqs;
    OneReadInputType** one_read;
};

typedef moodycamel::BlockingConcurrentQueue<fpga_data_input*> BlockingQueueStep1;
typedef bb::ByteBuffer ByteBuffer;
void mem_align1_core_pre(const mem_opt_t *opt, const bwaidx_t *idx, int l_seq, char* seq, TaskReadInputType* task_input, int temp, int *num_seed, void* buf);
void mem_align1_core_suf(const mem_opt_t *opt, const mem_chain_v *chn, int l_query, int64_t **rmax_data, fpga_sw_output** fpga_output, mem_alnreg_v *av);
void DataTobuffer(fpga_data_input* Data, ByteBuffer &buffer);
void writeDataToFPGA(Byte* buffer);
void TasktoFpga(const mem_opt_t *opt, const bwaidx_t *idx, char *seq, mem_chain_v *chn, int read_id, int task_id, int l_seq, int64_t **rmax );
#endif //BWA_FPGA_FPGA_SW_H
