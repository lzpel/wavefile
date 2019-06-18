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
    std::cout << "Hello, World!" << std::endl;
    return 0;
}