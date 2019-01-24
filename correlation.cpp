#include <complex>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include "correlation.hpp"

#define w_c 5 // common frequency


std::complex<double> gaussian(double db){
  double No = pow(10.0,-db/10.0);
  std::complex<double> I(0.0,1.0);
  double random1 = ((double)rand())/RAND_MAX;
  double random2 = ((double)rand())/RAND_MAX;
  double z1=sqrt(-2.0*log(random1)) * sin(2.0*M_PI*random2);
  double z2=sqrt(-2.0*log(random2)) * cos(2.0*M_PI*random1);
  return sqrt(No/2)*(z1 + z2*I);// + I*sqrt(oneside)*z2;
}

double mod1_co(double x){
  return x - (double)((int)x);
}

std::complex<double> autoCorrelation(std::complex<double> *code,int k,int length,int data=1){
  std::complex<double> sum = std::complex<double>(0.0);
  double ddata = (double)data;
  for(int i=0;i<length-k;i++){
    sum = sum + code[(i+k)%length]*conj(code[i%length]);
  }

  for(int i=length-k;i<length;i++){
    sum = sum + ddata*code[(i+k)%length]*conj(code[i%length]);
  }
	
  sum = sum/(double)length;
	
  return sum;
}

std::complex<double> crossCorrelation(std::complex<double> *code1,std::complex<double> *code2,int k,int length,int data=1){
  std::complex<double> sum = std::complex<double>(0.0);
  double ddata = (double)data;
  for(int i=0;i<length-k;i++){
    sum = sum + code1[(i+k)%length]*conj(code2[i%length]);
  }

  for(int i=length-k;i<length;i++){
    sum = sum + ddata*code1[(i+k)%length]*conj(code2[i%length]);
  }

  sum = sum/(double)length;
	
  return sum;
}

void match(std::complex<double> *code,std::complex<double> *channel,int *bit,int codelength,int bitlength,int interval,int late,std::complex<double> phase){
  //Assume the Power P = 1
  
  int channnelLocater;
  std::complex<double> sum;
  for(int i=0;i<bitlength;i++){
    sum = std::complex<double>(0.0,0.0);
    for(int j=0;j<codelength;j++){
      //それぞれについて総和をとる
      channnelLocater = late + interval*(i*codelength + j);
      for(int k=0;k<interval;k++){
	sum += channel[channnelLocater+k]*conj(code[j]*phase)/(double)codelength;
      }
    }
    //実部が正なら1 負なら−1
    if(sum.real() >= 0){
      bit[i] = 1;
    }else{
      bit[i] = -1;
    }
  }
}

void match_symbol(std::complex<double> *code,std::complex<double> *channel,std::complex<double> *symbol,int codelength,int symbollength,int interval,int late,std::complex<double> phase){
  //Assume the Power P = 1 and T = 1
  // width = 1/(codelength * interval)
  
  int channnelLocater;
  double width = 1.0/(double)(codelength * interval);
  std::complex<double> I(0.0,1.0);
  std::complex<double> sum;
  for(int i=0;i<symbollength;i++){
    //printf("bit num %d\n",i);
    sum = std::complex<double>(0.0,0.0);
    for(int j=0;j<codelength;j++){
      //printf("code num %d\n",j);
      //それぞれについて総和をとる
      channnelLocater = late + interval*(i*codelength + j);
      for(int k=0;k<interval;k++){
	//printf("interval %d\n",k);
	sum += width * channel[channnelLocater+k]*conj(code[j]*phase);
	//noise term
	//sum += gaussian(db)*exp(I*((double)(w_c*(interval*j+k)*width)))*conj(code[j]*phase)*width;
      }
    }
    symbol[i] = sum;
  }
}

void match_noise(std::complex<double> *code,std::complex<double> *channel,int *bit,int codelength,int bitlength,int interval,int late,std::complex<double> phase,double db){
  //Assume the Power P = 1 and T = 1
  // width = 1/(codelength * interval)
  
  int channnelLocater;
  double width = 1.0/(double)(codelength * interval);
  double coef = 1/sqrt(2);
  std::complex<double> I(0.0,1.0);
  std::complex<double> sum;
  for(int i=0;i<bitlength;i++){
    //printf("bit num %d\n",i);
    sum = std::complex<double>(0.0,0.0);
    for(int j=0;j<codelength;j++){
      //printf("code num %d\n",j);
      //それぞれについて総和をとる
      channnelLocater = late + interval*(i*codelength + j);
      for(int k=0;k<interval;k++){
	//printf("interval %d\n",k);
	sum += width * channel[channnelLocater+k]*conj(code[j]*phase);
	//noise term
	//sum += gaussian(db)*exp(I*((double)(w_c*(interval*j+k)*width)))*conj(code[j]*phase)*width;
      }
    }
    sum = sum + gaussian(db);
    //実部が正なら1 負なら−1
    if(sum.real() >= 0){
      bit[i] = 1;
    }else{
      bit[i] = -1;
    }
  }
}


void symbolMapping(int *symbol,std::complex<double> *transmit,int mode, int symbol_l,int transmit_l){
  int size = symbol_l;
  int map_size = transmit_l;
  std::complex<double> I(0.0,1.0);
  double power;
  int iq = (int)mode/2.0;
  switch(mode){
  case 2:
    power = 1.0/std::sqrt(2);
    break;
  case 4:
    power = 1.0/std::sqrt(10);
    break;
  default:
    return;
  }
  //#pragma omp parallel for 
  for(int n=0;n<map_size;n++){
    transmit[n] = 0.0 + 0.0*I;
  }
  //#pragma omp parallel for 
  for(int n=0;n<map_size;n++){
    transmit[n] = symbol_location(mode,symbol,mode*n) + symbol_location(mode,symbol,mode*n+iq)*I;
    transmit[n] = transmit[n]*power;
  }
  
}

double symbol_location(int mode,int *binary,int start){
  int iq = mode/2;
  int val=0;
  int tmp = 1;
  for(int i=0;i<iq;i++){
    val = val + tmp*binary[start+i];
    tmp = tmp*2;
  }
  if(mode == 2){
    // 4QAM, QPSK
    switch(val){
    case 0:
      return 1.0;
    case 1:
      return -1.0;
    default:
      return 0.0;
    }
  }else if(mode == 4){
    // 16 QAM
    switch(val){
    case 0:
      return 3.0;
    case 1:
      return -3.0;
    case 2:
      return 1.0;
    case 3:
      return -1.0;
    default:
      return 0.0;
    }
  }
  return 0.0;
}

void symbolFromDistance(int d,int *binary,int start,int mode){
  int iq = mode/2;
  if(mode == 2){
    switch(d){
    case 1:
      binary[start]=0;
      return;
    case -1:
      binary[start]=1;
      return;
    default:
      binary[start]=-1;
      std::cout << "Error at Symbol From Distance" << std::endl;
      return;
    }
  }else if(mode == 4){
    switch(d){
    case 3:
      binary[start]=0;
      binary[start+1]=0;
      return;
    case 1:
      binary[start]=0;
      binary[start+1]=1;
      return;
    case -1:
      binary[start]=1;
      binary[start+1]=1;
      return;
    case -3:
      binary[start]=1;
      binary[start+1]=0;
      return;
    default:
      binary[start]=-1;
      binary[start+1]=-1;
      return;
    }    
  }
}

void symbolRecovering(int *symbol,std::complex<double> *receive,int mode,int receive_l){
  double power;
  std::complex<double> I(0.0,1.0);
  int iq = (int)((double)mode/2.0);
  double shift = std::pow(2.0,iq)-1.0;
  int real_v, imag_v,tmp_r,tmp_i;
  switch(mode){
  case 2:
    power = std::sqrt(2);
    break;
  case 4:
    power = std::sqrt(10);
    break;
  default:
    std::cout << "Error at Symbol Recovering" << std::endl;
    return;
  }
  //std::cout << iq << std::endl;
  //#pragma omp parallel for private(real_v,imag_v,tmp_r,tmp_i)
  for(int n=0;n<receive_l;n++){
    tmp_r = (power*receive[n].real() + shift + 1.0)/2.0;
    tmp_i = (power*receive[n].imag() + shift + 1.0)/2.0;
    //std::cout << receive[n] << std::endl;
    real_v = (int)tmp_r;
    imag_v = (int)tmp_i;
    if(real_v < 0)
      real_v = 0;
    if(real_v >= (int)std::pow(2.0,iq)-1)
      real_v = (int)std::pow(2.0,iq)-1;
    if(imag_v < 0)
      imag_v = 0;
    if(imag_v >= (int)std::pow(2.0,iq)-1)
      imag_v = (int)std::pow(2.0,iq)-1;
    real_v = 2*(real_v - shift/2.0);
    imag_v = 2*(imag_v - shift/2.0);
    //if(std::abs(real_v)!=1 || std::abs(imag_v)!=1)
    //std::cout << real_v << " " << imag_v << std::endl;
    symbolFromDistance(real_v,symbol,mode*n,mode);
    symbolFromDistance(imag_v,symbol,mode*n+iq,mode);
    
  } 
}


void addGaussianNoise(std::complex<double>  *coefficient,double db, int mode,int coef_l){
  // TWO SIDE
  double power = std::pow(10,-db/10.0)/(2.0*log2((double)mode));
  double rand1,rand2,rand3,rand4,rand_x,rand_y;
  std::complex<double> I(0.0,1.0);
  double ave=0.0;
  //#pragma omp parallel for private (rand1,rand2,rand3,rand4,rand_x,rand_y)
  for(int x=0;x<coef_l;x++){
    // POWER N0/2
    rand1 = ((double)rand())/RAND_MAX;
    rand2 = ((double)rand())/RAND_MAX;
    rand3 = ((double)rand())/RAND_MAX;
    rand4 = ((double)rand())/RAND_MAX;
    //std::cout << rand1 << std::endl;
    rand_x = std::sqrt(power)*std::sqrt(-2.0*std::log(rand1))*std::cos(2.0*M_PI*rand2);
    rand_y = std::sqrt(power)*std::sqrt(-2.0*std::log(rand3))*std::sin(2.0*M_PI*rand4);
    //std::cout << std::log(rand1) << std::endl;
    coefficient[x] += rand_x + rand_y*I;
  }
  //std::cout << "random: power;" << ave/(double)cof_l << std::endl;
}




