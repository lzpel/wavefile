#include <iostream>
#include "wavefile.h"
#include "fftsg.h"

int main() {
	double *p, *p2, *t;
	long l;
	waveread("../test.wav", &l, &p);
	FFT fft(1024);
	fft.alloc(t);
	fft.fir(t,2,4,fft.size);
	fft.print("out0.csv",t,1);
	fft.rdft(t);
	fft.spectrum_cos(t);
	fft.print("out1.csv",t,1);
	fft.free(t);
	return 0;
}