//
// Created by misumi3104 on 2019/06/19.
//

#ifndef TEST_FFTSG_H
#define TEST_FFTSG_H

#include <math.h>
#include <float.h>
#include <cstring>

#ifndef MIN
#define MIN(a, b) (a)>(b)?(b):(a)
#endif
#ifndef MAX
#define MAX(a, b) (a)>(b)?(a):(b)
#endif

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

	void rdft(double *a) {
		::rdft(size, 1, a, ip, w);
	}

	void irdft(double *a) {
		::rdft(size, -1, a, ip, w);
		//sum_0^s f_t^2 =(sum_0^s/2 F_a^2)/(s/2)
		for (int j = 0; j < size; j++) a[j] *= 2.0 / size;
	}

	void alloc(double *&p) {
		alloc(p, size);
	}

	static void alloc(double *&p, int s) {
		zero(p = new double[s],s);
	}

	static void free(double *p) {
		delete[] p;
	}

	double power(double *d, bool spectrum) {
		return power(d, spectrum, size);
	}

	static double power(double *d, bool spectrum, int len) {
		double s=0,a=0;
		if(spectrum){
			for (int i = 1; i < len; ++i)s+=d[i]*d[i]/len*2;
		}else{
			for (int i = 0; i < len; ++i)a+=d[i]/len;
			for (int i = 0; i < len; ++i)s+=(d[i]-a)*(d[i]-a);
		}
		return s;
	}

	void copy(double *into, const double *from) {
		copy(into, from, size);
	}

	static void copy(double *into, const double *from, int len) {
		for (int i = 0; i < len; ++i)into[i] = from[i];
	}

	static void interpolate(double *into, int intolen, const double *from, int fromlen) {
		if (intolen > fromlen) {
			for (int i = 0; i < intolen; ++i) {
				double pos = 1.0 * fromlen * i / intolen;
				into[i] = from[int(pos)] * ((int) pos + 1 - pos) + from[int(pos) + 1] * (pos - (int) pos);
			}
		} else {
			for (int i = 0; i < intolen; ++i) {
				into[i] = from[fromlen * i / intolen];
			}
		}
	}

	void window(double *p) {
		for (int i = 0; i < size; ++i)p[i] *= 1 - cos(2 * M_PI * i / size);
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
		for (int i = 1; i < size / 2; ++i) w[i] = w[i * 2] * w[i * 2] + w[i * 2 + 1] * w[i * 2 + 1];
		for (int i = 0; i < size / 2; ++i) {
			w[i + size / 2] = 0;
			w[i] = w[i] ? sqrt(w[i]) : 0;
		}
	}

	void spectrum_log(double *w) {
		w[0] = w[1] = 0;
		for (int i = 1; i < size / 2; ++i) w[i] = w[i * 2] * w[i * 2] + w[i * 2 + 1] * w[i * 2 + 1];
		for (int i = 0; i < size / 2; ++i) {
			w[i + size / 2] = 0;
			w[i] = w[i] ? log(w[i]) / 2 : 0;
		}
	}

	void spectrum_aperiodicity(double *w, int n) {
		w[0] = w[size / 2] = 0;
		for (int i = 1; i < size / 2 / n; i++) {
			const int npos = i * n - (n - 1) / 2;
			double min = DBL_MAX, max = DBL_MIN;
			for (int j = npos; j < npos + n; j++) {
				min = MIN(min, w[j]);
				max = MAX(max, w[j]);
			}
			w[i] = max;
			w[i + size / 2] = min;
		}
		for (int i = size / 2 / n; i < size / 2; i++) {
			w[i] = w[i + size / 2] = 0;
		}
	}

	void zero(double *t) {
		zero(t, size);
	}

	static void zero(double *t, int s) {
		for (int i = 0; i < s; i++)t[i] = 0;
	}

	void white(double *d) {
		white(d, size);
	}

	static void white(double *d, int n) {
		for (int i = 0; i < n / 2; i++) {
			double r1 = (double) rand() / RAND_MAX, r2 = (double) rand() / RAND_MAX;
			d[i * 2 + 0] = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
			d[i * 2 + 1] = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2);
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

	static void fir(double *p, const int lenmax, const int lenmin, const int order) {
		//本当はi<=orderにして対称性が欲しいがフィルタ次数を奇数にしたくない。
		if (lenmax) {
			for (int i = 0; i < order; ++i)p[i] = sinc(i - order / 2, lenmax);
		} else {
			for (int i = 0; i < order; ++i)p[i] = (i == order / 2) ? 1 : 0;
		}
		if (lenmin) {
			for (int i = 0; i < order; i++)p[i] -= sinc(i - order / 2, lenmin);
		}
	}

	static signed zerocrosslen(int lenmin, int lenmax, double *p) {
		for (int i = 0; i < lenmax / 2; ++i) {
			if (p[i] * p[i + 1] < 0) {
				for (int j = i + lenmin; j < i + lenmax; j++) {
					if ((p[j] * p[j + 1] < 0) && (p[i] * p[j] > 0))return j - i;
				}
			}
		}
		return 0;
	}

	static void zerocrosslenarray(int lenmin, int lenmax, double *pd, double *ps, int pn) {
		for (int i = 0; i < pn - lenmax * 1.5; i++)pd[i] = zerocrosslen(lenmin, lenmax, ps + i);
	}
};

#endif //TEST_WAVE_H
