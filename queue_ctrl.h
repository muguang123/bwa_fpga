//#ifndef SRC_DUPLEX_QUEUE_DUPLEX_QUEUE_H_
//#define SRC_DUPLEX_QUEUE_DUPLEX_QUEUE_H_
#ifndef QUEUE_CTRL_H_
#define QUEUE_CTRL_H_

#define DS_PSIZE_EXP	17	// package size 128KiB
#define DS_PCAP_EXP		14	// package capacity 16K (16K - 1 actual available), 2GB total

#define US_PSIZE_EXP	10	// package size 4KiB
#define US_PCAP_EXP		18	// package capacity 64K (64K - 1 actual available)

#define DS_PSIZE		(1L << DS_PSIZE_EXP)
#define DS_PCAP			((1L << DS_PCAP_EXP) - 1)
#define DS_ADDR_LIMIT	(DS_BASE_ADDR + (1L << (DS_PSIZE_EXP + DS_PCAP_EXP)))

#define US_PSIZE		(1L << US_PSIZE_EXP)
#define US_PCAP			((1L << US_PCAP_EXP) - 1)
#define US_ADDR_LIMIT	(US_BASE_ADDR + (1L << (US_PSIZE_EXP + US_PCAP_EXP)))
#define INTS_NUM	(1024 * 1024 / 4)

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <iostream>
#include <unistd.h>

//using namespace std;

#define DS_BASE_ADDR	0x00000000
#define US_BASE_ADDR	0x90000000
//#define US_BASE_ADDR	0x00000000

#define MAP_SIZE (32*1024UL)
#define MAP_MASK (MAP_SIZE - 1)

#define DS_TAIL_REG		0
#define DS_HEAD_REG     1
#define DS_CNT_REG    	2
#define RSV0_REG        3
#define US_TAIL_REG	    4
#define US_HEAD_REG		5
#define US_CNT_REG		6

#define SOP0_CNT        8
#define SOP1_CNT        9
#define SOP2_CNT        10
#define SOP3_CNT        11
#define RESULT_CNT      12

static void *regBase;
static int regFd;
static int dmaH2C, dmaC2H;

static uint dsTail = DS_BASE_ADDR;
static uint usHead = US_BASE_ADDR;

static void writeReg(uint regNum, uint data)
{
    *(regNum + (uint *)regBase) = data;
}
static uint readReg(uint regNum)
{
    return *(regNum + (uint *)regBase);
}
static void *printDebugInfo(void* args){

    while(1){
        FILE *file;
        file = fopen("gen_test.log","a+");
        if (!file) perror("gen_test.log"), exit(1);
        usleep(1000000);
        fprintf(file, "---------------------------------\n");
        fprintf(file, "down_tail=%08x\n",readReg(DS_TAIL_REG));
        fprintf(file, "down_head=%08x\n",readReg(DS_HEAD_REG));
        fprintf(file, "down_cnt=%d\n",readReg(DS_CNT_REG));
        fprintf(file, "up_tail=%08x\n",readReg(US_TAIL_REG));
        fprintf(file, "up_head=%08x\n",readReg(US_HEAD_REG));
        fprintf(file, "up_cnt=%d\n",readReg(US_CNT_REG));
        fprintf(file, "SOP0_CNT=%d\n",readReg(SOP0_CNT));
        fprintf(file, "SOP1_CNT=%d\n",readReg(SOP1_CNT));
        fprintf(file, "SOP2_CNT=%d\n",readReg(SOP2_CNT));
        fprintf(file, "SOP3_CNT=%d\n",readReg(SOP3_CNT));
        fprintf(file, "RESULT_CNT=%d\n",readReg(RESULT_CNT));
        fprintf(file, "down_tail=%08x\n",readReg(DS_TAIL_REG));
        fclose(file);
    }

}



static int GetDsPkgCnt()
{
    return readReg(DS_CNT_REG);
}
static int IsDsWriteable(int pkgCnt)
{
    int cnt = GetDsPkgCnt();
    return cnt + pkgCnt <= DS_PCAP;
}
static void updateDsTail(int pkgCnt)
{
    for(; pkgCnt > 0; pkgCnt--)
    {
        dsTail += DS_PSIZE;
        if(dsTail >= DS_ADDR_LIMIT)
            dsTail = DS_BASE_ADDR;
    }
    writeReg(DS_TAIL_REG, dsTail);
}
static void WriteDsPkg(void *buf, uint pkgCnt)
{
    off_t off = lseek(dmaH2C, dsTail, SEEK_SET);
    int rc = write(dmaH2C, buf, DS_PSIZE * pkgCnt);
    if (rc < 0) {
        std::cerr<<"err when write DDR at address at "<<dsTail<< std::endl;
    }
/*    else{
        printf(" writePkg at DDR address : %p",dsTail);
    }*/
    updateDsTail(pkgCnt);
}

static int GetUsPkgCnt()
{
    return readReg(US_CNT_REG);
}
static int IsUsReadable(int pkgCnt)
{
    int cnt = GetUsPkgCnt();
    return cnt >= pkgCnt;
}
static void updateUsHead(int pkgCnt)
{
    for(; pkgCnt > 0; pkgCnt--)
    {
        usHead += US_PSIZE;
        if(usHead >= US_ADDR_LIMIT)
            usHead = US_BASE_ADDR;
    }
    writeReg(US_HEAD_REG, usHead);
}
static void ReadUsPkg(void *buf, uint pkgCnt)
{
    off_t off = lseek(dmaC2H, usHead, SEEK_SET);
    int rc = read(dmaC2H, buf, US_PSIZE * pkgCnt);
    if(rc < 0){
        std::cerr<<"err when read DDR at address at "<<usHead<< std::endl;
    }
/*    else{
        printf("readPkg at DDR address : %d",usHead); 
    }*/
    updateUsHead(pkgCnt);
}
static int OpenQueueCtrl()
{
    regFd = open("/dev/xdma0_user", O_RDWR | O_SYNC);
    regBase = mmap(0, MAP_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, regFd, 0);

    dmaH2C = open("/dev/xdma0_h2c_0", O_RDWR);
    dmaC2H = open("/dev/xdma0_c2h_0", O_RDWR | O_NONBLOCK);

    if(regFd >= 0 && dmaH2C >= 0 && dmaC2H >= 0){
        writeReg(RSV0_REG, 1);
        writeReg(DS_TAIL_REG, dsTail);
        writeReg(US_HEAD_REG, usHead);
        while(GetDsPkgCnt()||GetUsPkgCnt()) ;
        writeReg(RSV0_REG, 0);
        return 1;
    }
    else
        return 0;
}
static void CloseQueueCtrl()
{
    if(regFd >= 0)
        close(regFd);
    if(dmaH2C >= 0)
        close(dmaH2C);
    if(dmaC2H >= 0)
        close(dmaC2H);
}

#endif
