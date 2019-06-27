#include <iostream>
#include "wavefile.h"
#include "fftsg.h"

void test1(int f0){
	static const int n=3;
	double *p, *p2;
	signed l, r;
	waveload("../test.wav", &l, &r, &p);
	FFT fft(1024);
	fft.print("out10.csv", p, l);
	fft.alloc(p2,l);
	fft.fir(fft.buf0, r / 400, r / 200);
	fft.conv(p2, p, l, fft.buf0, true);
	wavesave("out10.wav", l, r, p2);
	for(int i=0;i<l;i++){
		if(i%fft.size==0) {
			signed t0;
			if (f0) {
				t0 = r / f0;
			} else {
				t0 = fft.zerocrosslen(r / 400, r / 200, p2 + i);
				t0 = MAX(r / 400, t0);
			}
			fft.interpolate(fft.buf0, fft.size, p + i, t0 * n);
			//fft.window(fft.buf0);窓を使うと非周期成分が見えない
			fft.rdft(fft.buf0);
			fft.spectrum_log(fft.buf0);
			fft.spectrum_log_periodseparate(fft.buf1, fft.buf2, fft.buf0, n);
			fft.print("out11.csv",fft.buf0);
			fft.print("out12.csv",fft.buf1);
			fft.print("out13.csv",fft.buf2);
		}
	}
	fft.free(p2);
}
void test2(){
	FFT fft(1024);
	double *t,power=0;
	fft.alloc(t);
	fft.white(t);
	printf("%f\n",fft.power(t,0));
	fft.rdft(t);
	printf("%f\n",fft.power(t,1));
	fft.spectrum_amp(t);
	printf("%f\n",fft.power(t,1));
	fft.free(t);
}
void test3(){

}
int main() {
	test1(0);
	test2();
	return 0;
}