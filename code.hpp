#include <iostream>
#include <complex>

void weyl(std::complex<double> *,int ,double,int);

void weyl(std::complex<double> *,double ,int);

void doubleWeyl(std::complex<double> *,double ,int);

void weylPrimitiveRoot(std::complex<double> *,int ,int ,int);

void primitiveRoot(std::complex<double> *,int ,int ,int);

void eigenvectorSequencePositive(std::complex<double> *,int,int,int);

void eigenvectorSequenceNegative(std::complex<double> *,int,int,int);

void eigenvectorSequenceNegative2(std::complex<double> *,int,int);

void boolMap(std::complex<double> *,double ,double,int);

void chaos(std::complex<double> *,int,double,int);

void cheby_opt(std::complex<double> *,double ,int,double,int);

void kohda(std::complex<double> *,int, int, double, int);

void gold_IQ(std::complex<double> *,int,int,int);

void gold(std::complex<double> *,int,int,int);

void gold_weyl(std::complex<double> *,double,int,int,int);

void FZC(std::complex<double> *,int,double,double,double,int);

void FZC(std::complex<double> *,int,int);

void addWeyl(std::complex<double> *,double,int);

void chaos_weyl(std::complex<double> *,int,double,double,int);

void half_weyl_chaos(std::complex<double> *,int ,int,double,int);

void primitiveRoot_weyl(std::complex<double> *,int ,int ,double, int);

void addAngle(std::complex<double> *,double, int);

void spSequence(std::complex<double> *,int ,int );

double vdc(int , double);

void vdcSequence(std::complex<double> *,double ,int);

void LSF(std::complex<double> *,double ,int);

void realLSF(std::complex<double> *,double ,int);

void mSequence(int ,int *,int *,int *,int);
void mSeq(std::complex<double> *,int,int);

void weyl01(std::complex<double> *,double,int);

void transform1(std::complex<double> *,int length);

double mod1(double x);

class trigonometric
{
private:
	double sinM(int,double,double);
	double cosM(int,double,double);
public:
	double sinValue,cosValue;
	 trigonometric(double);
	 ~trigonometric();
	 void mAngle(int);

};
