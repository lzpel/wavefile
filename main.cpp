#include <iostream>
#include "wavefile.h"
#include "fftsg.h"

void test1(){
	double *p, *t;
	signed l, r;
	waveload("../test.wav", &l, &r, &p);
	FFT fft(1024);
	fft.print("out10.csv", p, 1, l);
	fft.alloc(t);
	fft.fir(t, r / 400, r / 200);
	fft.conv(p, p, l, t, true);
	wavesave("out10.wav", l, r, p);
	for(int i=0;i<(l/fft.size);i++){
		double* pi=p+i*fft.size;
		int t0=fft.min(r/400, fft.zerocrosslen(r/400,r/200,pi));
		fft.interpolate(t,fft.size,pi,t0*2);
		//fft.window(t);窓を使うと非周期成分が見えない
		fft.rdft(t);
		fft.spectrum_log(t);
		fft.copy(pi,t);
	}
	fft.print("out11.csv", p, 1, l);
	fft.free(t);
}
int main() {
	test1();
	return 0;
}