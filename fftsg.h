//
// Created by misumi3104 on 2019/06/19.
//

#ifndef TEST_FFTSG_H
#define TEST_FFTSG_H

#include <math.h>

void cdft(int, int, double *, int *, double *);

void rdft(int, int, double *, int *, double *);

void ddct(int, int, double *, int *, double *);

void ddst(int, int, double *, int *, double *);

void dfct(int, double *, double *, int *, double *);

void dfst(int, double *, double *, int *, double *);

class FFT {
public:
	double *buf, *w;
	signed size;
	signed *ip;

	FFT(int s) {
		size = s;
		ip = new int[2 + (int)sqrt(s / 2)];
		buf = new double[s];
		w = new double[s / 2];
		ip[0] = 0;
	}

	~FFT() {
		delete[] ip;
		delete[] buf;
		delete[] w;
	}

	void rdft(double *a) {
		::rdft(size, 1, a, ip, w);
	}

	void irdft(double *a) {
		::rdft(size, -1, a, ip, w);
		for (int j = 0; j < size; j++) a[j] *= 2.0 / size;
	}

	void alloc(double*&p) {
		p=new double[size];
	}

	void free(double *p) {
		delete[] p;
	}

	void copy(double *into, const double *from) {
		for (int i = 0; i < size; ++i)into[i] = from[i];
	}
};

#endif //TEST_WAVE_H