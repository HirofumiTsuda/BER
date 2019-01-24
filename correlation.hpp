#include <complex>

std::complex<double> autoCorrelation(std::complex<double> *,int,int,int);

std::complex<double>  crossCorrelation(std::complex<double> *,std::complex<double> *,int,int,int);

void match(std::complex<double> *,std::complex<double> *,int *,int ,int ,int ,int,std::complex<double>);

void match_noise(std::complex<double> *,std::complex<double> *,int *,int ,int ,int ,int,std::complex<double>,double);

void match_symbol(std::complex<double> *,std::complex<double> *,std::complex<double> *,int ,int ,int ,int ,std::complex<double> );

double correlation(std::complex<double> *,int *,double ,int, std::complex<double> *,int *,double ,int ,int );

double symbol_location(int,int *,int);

void symbolFromDistance(int ,int *,int ,int );

void symbolMapping(int *,std::complex<double> *,int, int,int);

void symbolRecovering(int *,std::complex<double> *,int,int);

void addGaussianNoise(std::complex<double>  *,double, int ,int);
