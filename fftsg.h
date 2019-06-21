//
// Created by misumi3104 on 2019/06/19.
//

#ifndef TEST_FFTSG_H
#define TEST_FFTSG_H

#include <math.h>
#include <cstring>

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
		ip = new int[2 + (int) sqrt(s / 2)];
		buf = new double[s];
		w = new double[s / 2];
		ip[0] = 0;
	}

	~FFT() {
		delete[] ip;
		delete[] buf;
		delete[] w;
	}

	static int min(const int a, const int b) {
		return a > b ? a : b;
	}

	void rdft(double *a) {
		::rdft(size, 1, a, ip, w);
	}

	void irdft(double *a) {
		::rdft(size, -1, a, ip, w);
		for (int j = 0; j < size; j++) a[j] *= 2.0 / size;
	}

	void alloc(double *&p) {
		p = new double[size];
	}

	void free(double *p) {
		delete[] p;
	}

	void copy(double *into, const double *from) {
		copy(into, from, size);
	}

	static void copy(double *into, const double *from, int len) {
		for (int i = 0; i < len; ++i)into[i] = from[i];
	}

	void print(const char *fn, const double *p, bool row) {
		print(fn, p, row, size);
	}

	static void print(const char *fn, const double *p, bool row, int len) {
		const char *fmt = row ? "%f\n" : "%f ";
		if (fn) {
			FILE *f = fopen(fn, "w");
			for (int i = 0; i < len; i++)fprintf(f, fmt, p[i]);
			fclose(f);
		} else {
			for (int i = 0; i < len; i++)printf(fmt, p[i]);
		}
	}

	void spectrum_cos(double *w) {
		w[0] = w[1] = 0;
		for (int i = 1; i < size / 2; ++i) {
			w[i * 2] = sqrt(w[i * 2] * w[i * 2] + w[i * 2 + 1] * w[i * 2 + 1]);
			w[i * 2 + 1] = 0;
		}
	}

	void conv(double *d, double *s, const int sn, const double *f, bool fh) {
		conv(d, s, sn, f, fh, size);
	}

	static void conv(double *d, double *s, const int sn, const double *f, bool fh, const int fn) {
		for (int i = 0; i < sn - fn; ++i) {
			double sum = 0;
			for (int j = 0; j < fn; ++j)sum += s[i + j] * f[j];
			d[i] = sum;
		}
		for (int i = sn - fn; i < sn; ++i)d[i] = 0;
		if (fh) {
			for (int i = sn - 1; i >= fn / 2; --i)d[i] = d[i - fn / 2];
			for (int i = fn / 2 - 1; i >= 0; --i)d[i] = 0;
		}
	}

	static inline double sinc(int i, int len) {
		//i%(len/2)==0の時return0;
		return (i == 0) ? (2.0 / len) : sin(i * M_PI * 2.0 / len) / (i * M_PI);
	}

	void fir(double *p, int lenmin, int lenmax) {
		fir(p, lenmin, lenmax, size);
	}

	static void fir(double *p, const int lenmin, const int lenmax, const int order) {
		//本当はi<=orderにして対称性が欲しいがフィルタ次数を奇数にしたくない。
		if (lenmin) {
			for (int i = 0; i < order; i++) {
				p[i] = sinc(i - order / 2, lenmin);
			}
		} else {
			for (int i = 0; i < order; i++) {
				p[i] = (i == order / 2) ? 1 : 0;
			}
		}
		if (lenmax) {
			for (int i = 0; i < order; i++) {
				p[i] -= sinc(i - order / 2, lenmax);
			}
		}
	}

	static signed zerocrosslen(int lenmax, double *p) {
		int count = 0, pos[3]={0};
		for (int i = 0; i < lenmax; ++i)if ((p[i] * p[i + 1] < 0) && (count < 3))pos[count++]=i;
		return pos[2]?pos[2]-pos[0]:0;
	}

	static void zerocrosslenarray(int lenmax, double *pd, double *ps, int pn) {
		for (int i = 0; i < pn - lenmax; i++) {
			pd[i] = zerocrosslen(lenmax, ps + i);
		}
	}
};

#endif //TEST_WAVE_H
