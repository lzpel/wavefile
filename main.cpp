#include <iostream>
#include "wavefile.h"
#include "fftsg.h"

void test1(int f0){
	static const int n=3;
	double *p, *p2, *t;
	signed l, r;
	waveload("../test.wav", &l, &r, &p);
	FFT fft(1024);
	fft.print("out10.csv", p, 1, l);
	fft.alloc(p2,l);
	fft.alloc(t);
	fft.fir(t, r / 400, r / 200);
	fft.conv(p2, p, l, t, true);
	wavesave("out10.wav", l, r, p2);
	for(int i=0;i<(l/fft.size);i++){
		signed pos=i*fft.size,t0;
		if(f0){
			t0=r/f0;
		}else{
			t0=fft.zerocrosslen(r/400,r/200,p2+pos);
		}
		fft.interpolate(t,fft.size,p+pos,t0*n);
		//fft.window(t);窓を使うと非周期成分が見えない
		fft.rdft(t);
		fft.spectrum_log(t);
		fft.spectrum_aperiodicity(t,n);
		fft.copy(p+pos,t);
	}
	fft.print("out11.csv", p, 1, l);
	fft.free(t);
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
	fft.spectrum_cos(t);
	printf("%f\n",fft.power(t,1));
	fft.free(t);
}
int main() {
	//test1(306);
	test2();
	return 0;
}