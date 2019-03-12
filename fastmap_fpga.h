//
// Created by jie on 3/9/19.
//

#ifndef BWA_FPGA_FASTMAP_FPGA_H
#define BWA_FPGA_FASTMAP_FPGA_H

#include <iostream>
#include <functional>

#define BLOCKQUEUE_STEP0_SIZE 5
#define BLOCKQUEUE_STEP2_SIZE 100
#define BLOCKQUEUE_STEP3_SIZE 5

#define PACKAGE_SIZE 1*128*1024 //sizeof(char)
#define SEND_DELAY 1000000//超时时间1s

#define FPGA_OUT_PACKAGE_SIZE 1024 //the output of FPGA is 1k

struct ThreadPoolParameter{
    unsigned int Step0ThreadNum;
    unsigned int Step1ThreadNum;
    unsigned int SendToFpgaNum;
    unsigned int GetFromFpgaNum;
    unsigned int Step2ThreadNum;
    unsigned int Step3ThreadNum;
    std::function<void(void*)> Step0Func; // get 10000read
    std::function<void(void*)> Step1Func; // process 1 read, read to FPGA
    std::function<void(void*)> SendToFpgaFunc; // FPGA return to deal
    std::function<void(void*)> GetFromFpgaFunc; // FPGA return to deal
    std::function<void(void*)> Step2Func; // FPGA return to deal
    std::function<void(void*)> Step3Func;

    ThreadPoolParameter(int f0, int f1,  int s, int g, int f2, int f3, std::function<void(void*)> f0F,std::function<void(void*)> f1F,std::function<void(void*)> sF,std::function<void(void*)> gF,std::function<void(void*)> f2F, std::function<void(void*)> f3F)
    {
        Step0ThreadNum = f0;
        Step1ThreadNum = f1;
        SendToFpgaNum = s;
        GetFromFpgaNum = g;
        Step2ThreadNum = f2;
        Step3ThreadNum = f3;
        Step0Func = f0F;
        Step1Func = f1F;
        SendToFpgaFunc = sF;
        GetFromFpgaFunc = gF;
        Step2Func = f2F;
        Step3Func = f3F;
    };

    ThreadPoolParameter(){
        Step0ThreadNum = 0;
        Step1ThreadNum = 0;
        SendToFpgaNum = 0;
        GetFromFpgaNum = 0;
        Step2ThreadNum = 0;
        Step3ThreadNum = 0;
        Step0Func = 0;
        Step1Func = 0;
        SendToFpgaFunc = 0;
        GetFromFpgaFunc = 0;
        Step2Func = 0;
        Step3Func = 0;
    }
};

void bwa_mem_fpga(void);

#endif 
