#include <iostream>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <math.h>

#include "code.hpp"

//Gold////////
void init(int,int *);
void setPoly(int,int *,int);
int output(int,int *,int *);
double bitConvert(int );
void shift(int,int *,int *);
//////////////


void weyl(std::complex<double> *code,int p,double m,int length){
  double rootP = sqrt(p);
  rootP = rootP - (int)rootP;
  std::complex<double> I = std::complex<double>(0.0,1.0);
  for(int n=1; n<=length; n++) {
    //2*pi*n*m*sqrt()
    code[n-1] = exp(2.0*M_PI*I*(double)n*m*rootP);
  }
}

void weyl(std::complex<double> *code,double q,int length){
  q = q - (int)q;
  std::complex<double> I = std::complex<double>(0.0,1.0);
  for(int n=1; n<=length; n++) {
    //2*pi*n*m*sqrt()
    code[n-1] = exp(2.0*M_PI*I*(double)n*q);
  }
}

void kohda(std::complex<double> *code,int k, int number, double init,int length){
  // Quantization
  double kd = (double)k;
  double bit;
  if(std::abs(init) > 1){
    init = init - int(init);
  }
  double tmp = init;
  for(int n=0;n<length;n++){
    tmp = cos(mod1(std::abs(kd*acos(tmp))));
    bit = pow(2,number-1)*std::abs(tmp);
    bit = bit - (int)bit;
    if(bit < 0.5){
      code[n] = -1;
    }else{
      code[n] = 1;
    }
  }
}

void eigenvectorSequencePositive(std::complex<double> *code,int k,int l,int N){
  std::complex<double> I = std::complex<double>(0.0,1.0);
  for(int n=0;n<N;n++){
    code[(n*l)%N] = exp(2.0 * M_PI * I * mod1(((double)(n*k))/N));
  }
}

void eigenvectorSequenceNegative(std::complex<double> *code,int k,int l,int N){
  std::complex<double> I = std::complex<double>(0.0,1.0);
  for(int n=0;n<N;n++){
    code[(n*l)%N] = exp(2.0 * M_PI * I * mod1(((double)(n*(2.0*k+1)))/(2.0*N)));
  }
}

void eigenvectorSequenceNegative2(std::complex<double> *code,int k,int N){
  std::complex<double> I = std::complex<double>(0.0,1.0);
  double head = 1.0;
  double iter = -1.0;
  for(int n=0;n<N;n++){
    code[n] = head * exp(2.0 * M_PI * I * mod1(((double)(n*(2.0*k+1.0)))/(2.0*N)));
    head = head*iter;
  }
}

void doubleWeyl(std::complex<double> *code,double q,int half_length){
  q = q - (int)q;
  std::complex<double> I = std::complex<double>(0.0,1.0);
  for(int n=1; n<=half_length; n++) {
    //2*pi*n*m*sqrt()
    code[n-1] = exp(2.0*M_PI*I*(double)n*q);
  }

  for(int n=1; n<=half_length; n++) {
    //2*pi*n*m*sqrt()
    code[half_length+n-1] = exp(2.0*M_PI*I*(double)n*q);
  }
}

void half_weyl_chaos(std::complex<double> *code,int l,int maxLength,double q,int length){
  gold(code,l,length,maxLength);
  std::complex<double> I = std::complex<double>(0.0,1.0);
  //半分だけ使う
  for(int i=0;i<length/2;i++){
    code[2*i] = code[i];
  }
  for(int i=0;i<length/2;i++){
    code[2*i+1] = exp(2.0*M_PI*I*((double)i)*q);
  }
}



double mod1(double x) {
  return x - (int)x;
}



void gold_weyl(std::complex<double> *code,double q,int l,int length,int maxLength=1023){
  gold(code,l,length,maxLength);
  std::complex<double> I = std::complex<double>(0.0,1.0);
  for(int n=0;n<length;n++){
    code[n] = code[n]*exp(2.0*M_PI*I*((double)(n+1))*q);
  }
  return;
}

void weylPrimitiveRoot(std::complex<double> *code,double q,int root,int length){
  std::complex<double> I = std::complex<double>(0.0,1.0);
  int prim = 1;
  for(int n=1; n<=length; n++) {
    //2*pi*n*m*sqrt()
    code[n-1] = exp(2.0*M_PI*I*((double)prim)*q);
    prim = (root*prim)%(length+1);
  }
}

void addAngle(std::complex<double> *code,double theta, int length){
  std::complex<double> I = std::complex<double>(0.0,1.0);
  std::complex<double> rotation = exp(2.0*M_PI*I*theta);
  for(int i=0;i<length;i++){
    code[i] = code[i] * rotation;
  }
}



void primitiveRoot(std::complex<double> *code,int q,int k,int p){
  code[0] = 1.0;
  int prim = 1.0;
  std::complex<double> I = std::complex<double>(0.0,1.0);
  for(int n=1;n<p;n++){
    code[n] = exp(2.0 * M_PI * I * (double)(k* prim) / (double)p);
    prim = (prim * q)%p;
  }
}

void primitiveRoot_weyl(std::complex<double> *code,int q,int k,double init,int p){
  primitiveRoot(code,q,k,p);
  std::complex<double> I = std::complex<double>(0.0,1.0);
  for(int n=0;n<p;n++){
    code[n] = code[n] * exp(2.0 * M_PI * I * ((double)n) * init);
    std::cout << code[n] << std::endl;
  }
}

void LSF(std::complex<double> *code,double r,int length){
  std::complex<double> M(0.0,0.0);
  std::complex<double> tmp;
  M = code[0];
  for(int i=1;i<length;i++){
    tmp = code[i];
    code[i] = code[i] + r*M;
    M = r*M + tmp;
  }
  for(int i=0;i<length;i++){
    code[i] = code[i]/std::abs(code[i]);
  }
}

void realLSF(std::complex<double> *code,double r,int length){
  std::complex<double> M(0.0,0.0);
  std::complex<double> tmp;
  M = code[0];
  for(int i=1;i<length;i++){
    tmp = code[i];
    code[i] = code[i] + r*M;
    M = r*M + tmp;
  }
  std::complex<double> num = 0;
  for(int i=0;i<length;i++){
    code[i] = code[i].real();
    num = num + code[i]*code[i];
  }
  double norm = sqrt(std::abs(num.real())/length);
  for(int i=0;i<length;i++){
    code[i] = code[i]/norm;
  }
}

void FZC(std::complex<double> *code,int M,double p,double q,double r,int length){
  std::complex<double> I = std::complex<double>(0.0,1.0);
  for(int n=1;n<=length;n++){
    code[n-1] = pow(-1,n*M)*exp(I*M_PI*((pow(M,p)*pow(n,q)+pow(n,r))/(double)length));
  }
}

void FZC(std::complex<double> *code,int M,int length){
  std::complex<double> I = std::complex<double>(0.0,1.0);
  for(int n=1;n<=length;n++){
    if(M%2 == 0){
      code[n-1] = exp(-1.0*I*(M_PI*(double)(M*n*n)/length));
    }else{
      code[n-1] = exp(-1.0*I*(M_PI*(double)(M*n*(n+1))/length));
    }
  }
}

void spSequence(std::complex<double> *code,int k,int length){
  int L = length/2 - 1;
  std::complex<double> I(0.0,1.0);
  for(int l=0;l<length;l++){
    code[l] = pow(-1,l)*exp(2.0*I*M_PI*(l*k/(double)L));
  }
}


double vdc(int n, double base = 2.0){
  double vdc = 0, denom = 1;
  while (n){
    vdc += fmod(n, base) / (denom *= base);
    n /= base; // note: conversion from 'double' to 'int'
  }
  return vdc;
}

void vdcSequence(std::complex<double> *code,double base,int length){
  std::complex<double> I(0.0,1.0);
  for(int i=0;i<length;i++){
    code[i] = exp(2.0 * M_PI * I * vdc(i,base));
  }

}

void addWeyl(std::complex<double> *code,double init,int length){
  std::complex<double> I(0.0,1.0);
  for(int n=0;n<length;n++){
    code[n] = code[n]*exp(2.0*M_PI*I*((double)(n+1))*init);
  }
}

void boolMap(std::complex<double> *code,double a,double theta0,int length){
  std::complex<double> I = std::complex<double>(0.0,1.0);
  double theta = theta0 - (int)theta0;
  for(int i=0; i<length; i++) {
    theta = atan(tan(2.0*theta)/(2.0*a));
    code[i] = exp(2.0*I*theta);
  }
}

void chaos(std::complex<double> *code,int m,double init,int length){
  std::complex<double> I = std::complex<double>(0.0,1.0);
  trigonometric tri(2.0*M_PI*init);
  for(int i=0; i<length; i++) {
    code[i] = tri.cosValue + tri.sinValue*I;
    tri.mAngle(m);
  }
}

void cheby_opt(std::complex<double> *code,double r,int m,double init,int length){
  double init_tmp = init;
  double tmp;
  double norm = 0.0;
  for(int n=0;n<length;n++){
    trigonometric tri(2.0*M_PI*init_tmp);
    tmp = 0;
    for(int j=0;j<length;j++){
      tri.mAngle(m);
      tmp = tmp + std::pow(r,j)*tri.cosValue;
    }
    code[n] = tmp;
    init_tmp = mod1(init_tmp*(double)m);
  }
  for(int n=0;n<length;n++){
    norm = norm + (code[n]*std::conj(code[n])).real();
  }
  for(int n=0;n<length;n++){
    code[n] = code[n]*std::sqrt(length/norm);
  }
}

void chaos_weyl(std::complex<double> *code,int m,double init,double q,int length){
  //First, create the chaos sequence.
  std::complex<double> I(0.0,1.0);
  init = mod1(init);
  trigonometric tri(2.0*M_PI*init);
  for(int n=0;n<length;n++){
    code[n] = tri.cosValue + tri.sinValue*I;
    tri.mAngle(m);
  }
  //Second, apply weyl sequence to "code"
  for(int n=0;n<length;n++){
    code[n] = code[n]*exp(2.0*M_PI*I*(double)(n+1)*q);
  }
}

void mSequence(int bitlength,int *code,int *poly,int *box,int mode=0){
  init(bitlength,box);
  setPoly(bitlength,poly,mode);
  for(int i=0; i<=pow(2,bitlength)-1; i++) {
    code[i] = bitConvert(box[0]);
    shift(bitlength,box,poly);
  }


  return;

}

void gold_IQ(std::complex<double> *code,int l,int length,int maxLength=1023){
  int bitlength = (int)(log10(maxLength+1)/log10(2));
  int poly[bitlength];
  int box[bitlength];
  int mseq1[maxLength];
  int mseq2[maxLength];
  int code1[maxLength];
  int code2[maxLength];
  std::complex<double> I(0.0,1.0);

  mSequence(bitlength,mseq1,poly,box,0);
  mSequence(bitlength,mseq2,poly,box,1);


  for(int i=0; i<length; i++) {
    code1[i] = mseq1[i]*mseq2[(i+11)%length];

  }

  for(int i=0; i<length; i++) {
    code2[i] = mseq1[i]*mseq2[(i+5)%length];
  }

  for(int i=0; i<length; i++) {
    code[i] = (double)(code1[i]*code2[(i+7+l)%length])/sqrt(2) + ((double)code1[i]*code2[(i+l+19)%length])*I/sqrt(2);
  }


  return;
}

void gold(std::complex<double> *code,int l,int length,int maxLength=1023){
    int bitlength = (int)(log10(maxLength+1)/log10(2));
  int poly[bitlength];
  int box[bitlength];
  int mseq1[maxLength];
  int mseq2[maxLength];
  std::complex<double> I(0.0,1.0);

  mSequence(bitlength,mseq1,poly,box,0);
  mSequence(bitlength,mseq2,poly,box,1);


  for(int i=0; i<length; i++) {
    code[i] = mseq1[i]*mseq2[(i+l)%length];
  }
}

void mSeq(std::complex<double> *code,int length,int mode){
  int bitlength = (int)(log10(length+1)/log10(2));
  int poly[bitlength];
  int box[bitlength];
  int seq[length];

  mSequence(bitlength,seq,poly,box,mode);
  for(int i=0; i<length; i++) {
    code[i] = (double)seq[i];
  }

}

trigonometric::trigonometric(double theta){
  this->sinValue = sin(theta);
  this->cosValue = cos(theta);
}

trigonometric::~trigonometric(){
}

double trigonometric::sinM(int m,double sinT,double cosT){
  if(m == 1) {
    return sinT;
  }else{
    return sinM(m-1,sinT,cosT)*cosM(1,sinT,cosT) + cosM(m-1,sinT,cosT)*sinM(1,sinT,cosT);
  }
}

double trigonometric::cosM(int m,double sinT,double cosT){
  if(m == 1) {
    return cosT;
  }else{
    return cosM(m-1,sinT,cosT)*cosM(1,sinT,cosT) - sinM(m-1,sinT,cosT)*sinM(1,sinT,cosT);
  }
}

void trigonometric::mAngle(int m){
  double cosMAngle,sinMAngle,r;
  sinMAngle = sinM(m,sinValue,cosValue);
  cosMAngle = cosM(m,sinValue,cosValue);
  r = sqrt(sinMAngle*sinMAngle + cosMAngle*cosMAngle);
  sinValue = sinMAngle/r;
  cosValue = cosMAngle/r;
}

////////////////////////////////
//Gold
/////////////////////////////

void init(int bitlength,int *box){
  for(int i=0; i<bitlength; i++) {
    box[i] = 1;
  }
}

void weyl01(std::complex<double> *code,double init,int length){
  double tmp = init;
  std::complex<double> I(0.0,1.0);
  for(int i=0;i<length;i++){
    code[i] = tmp;
    tmp = mod1(tmp+init);
    code[i] += tmp*I;
    tmp = mod1(tmp+init);
  }
}

void transform1(std::complex<double> *code,int length){
  double r = std::sqrt(length);
  double realP[2*length];
  std::complex<double> I(0.0,1.0);
  double tmp = r;
  for(int n=0;n<length;n++){
    realP[2*n] = std::real(code[n]);
    realP[2*n+1] = std::imag(code[n]);
  }
  for(int n=0;n<length-1;n++){
    code[n] = tmp * cos(2.0*M_PI*realP[2*n]) + tmp * sin(2.0*M_PI*realP[2*n]) * cos(2.0*M_PI*realP[2*n+1])*I;
    tmp = tmp * sin(2.0*M_PI*realP[2*n]) * sin(2.0*M_PI*realP[2*n+1]);

  }
   code[length-1] = tmp * cos(2.0*M_PI*realP[2*length-2]) + tmp * sin(2.0*M_PI*realP[2*length-2]) * I;
}


void setPoly(int bitlength,int *poly,int mode=0){
  switch(bitlength) {
  case 4:
    poly[0]=1;
    poly[1]=1;
    poly[2]=0;
    poly[3]=0;
    break;
  case 5:
    if(mode == 0) {
      poly[0]=1;
      poly[1]=0;
      poly[2]=1;
      poly[3]=0;
      poly[4]=0;
      break;
    }else{
      poly[0]=1;
      poly[1]=1;
      poly[2]=1;
      poly[3]=1;
      poly[4]=0;
      break;
    }
  case 6:
    if(mode==0) {
      poly[0]=1;
      poly[1]=1;
      poly[2]=0;
      poly[3]=0;
      poly[4]=0;
      poly[5]=0;
      break;
    }else{
      poly[0]=1;
      poly[1]=1;
      poly[2]=1;
      poly[3]=0;
      poly[4]=0;
      poly[5]=1;
      break;
    }
  case 7:
    if(mode==0) {
      poly[0]=1;
      poly[1]=1;
      poly[2]=0;
      poly[3]=0;
      poly[4]=0;
      poly[5]=0;
      poly[6]=0;
      break;
    }else{
      poly[0]=1;
      poly[1]=0;
      poly[2]=0;
      poly[3]=1;
      poly[4]=0;
      poly[5]=0;
      poly[6]=0;
      break;
    }
  case 8:
    if(mode==0) {
      poly[0]=1;
      poly[1]=1;
      poly[2]=1;
      poly[3]=0;
      poly[4]=0;
      poly[5]=0;
      poly[6]=0;
      poly[7]=1;
      break;
    }else{
      poly[0]=1;
      poly[1]=1;
      poly[2]=1;
      poly[3]=1;
      poly[4]=0;
      poly[5]=0;
      poly[6]=1;
      poly[7]=1;
      break;
    }
  case 9:
    if(mode==0) {
      poly[0]=1;
      poly[1]=0;
      poly[2]=0;
      poly[3]=0;
      poly[4]=0;
      poly[5]=1;
      poly[6]=0;
      poly[7]=0;
      poly[8]=0;
      break;
    }else{
      poly[0]=1;
      poly[1]=0;
      poly[2]=0;
      poly[3]=1;
      poly[4]=0;
      poly[5]=1;
      poly[6]=1;
      poly[7]=0;
      poly[8]=0;
      break;
    }
  case 10:
    if(mode==0) {
      poly[0]=1;
      poly[1]=0;
      poly[2]=0;
      poly[3]=0;
      poly[4]=0;
      poly[5]=0;
      poly[6]=0;
      poly[7]=1;
      poly[8]=0;
      poly[9]=0;
      break;
    }else{
      poly[0]=1;
      poly[1]=0;
      poly[2]=1;
      poly[3]=0;
      poly[4]=0;
      poly[5]=0;
      poly[6]=0;
      poly[7]=1;
      poly[8]=1;
      poly[9]=0;
      break;
    }
  default:
    exit(1);
    break;
  }
}

int output(int bitlength,int *box,int *poly){
  int sum = 0;
  for(int i=0; i<bitlength; i++) {
    sum = sum + box[i]*poly[bitlength-1-i];
  }
  return sum%2;
}

double bitConvert(int a){
  if(a==1)
    return 1.0;

  return -1.0;
}

void shift(int bitlength,int *box,int *poly){
  //まず出力を計算
  int out;
  for(int i=0; i<bitlength; i++) {
    out = output(bitlength,box,poly);
  }
  //シフトする
  for(int i=bitlength-1; i>0; i--) {
    box[i] = box[i-1];
  }
  box[0] = out;
}
