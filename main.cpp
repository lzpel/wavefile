#include <iostream>
#include "wavefile.h"
#include "fftsg.h"

int main() {
	double *p, *t;
	long l;
	waveread("../test.wav", &l, &p);
	FFT fft(1024);
	fft.print("out0.csv",p,1,l);
	fft.alloc(t);
	fft.fir(t,44100/200,44100/400);
	fft.conv(p,p,l,t);
	fft.free(t);
	fft.print("out1.csv",p,1,l);
	fft.zerocrosslenarray(44100/200,p,p,l);
	fft.print("out2.csv",p,1,l);
	return 0;
}