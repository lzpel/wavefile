#include <iostream>
#include "wavefile.h"
#include "fftsg.h"
int main() {
    double*p,*t;
    long l;
    waveread("../test.wav",&l,&p);
    FFT fft(8192);
    fft.alloc(t);
    fft.copy(t,p);
	printf("%f\t%f\t%f\t%f\t%f\t\n",t[0],t[1],t[2],t[3],t[4]);
    fft.rdft(t);
	printf("%f\t%f\t%f\t%f\t%f\t\n",t[0],t[1],t[2],t[3],t[4]);
	fft.irdft(t);
	printf("%f\t%f\t%f\t%f\t%f\t\n",t[0],t[1],t[2],t[3],t[4]);
    fft.free(t);
    fft.firlpf(NULL,NULL,300,1200);
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
void    hogehoge(int fftsize,int length,int step,double*p){
	double*t;
	FFT fft(fftsize);
	fft.alloc(t);
	for(int i=0;i<=length/step;i++){
		fft.copy(t,p+FFT::min(step*i,length-fftsize));
		fft.rdft(t);

		fft.irdft(t);
	}
	fft.copy(t,p);

}