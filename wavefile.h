//
// Created by misumi3104 on 2019/06/19.
//

#ifndef TEST_WAVEFILE_H
#define TEST_WAVEFILE_H
#define WAVEFILE_LOAD(FP, T, CH, L, D) for(int i=0;i<(L);i++){T m,s=0;if(CH==1){fread(&s,sizeof(T),1,FP);}else{for(int j=0;j<CH;j++){fread(&m,sizeof(T),1,FP);s+=m;}s/=CH;}(D)[i]=s;}
#define WAVEFILE_SAVE(FP, T, CH, L, D) for(int i=0;i<(L);i++){T s=(T)D[i];for(int j=0;j<CH;j++)fwrite(&s,sizeof(T),1,FP);}

#include<stdio.h>

bool waveriff(const char *riff, char *p) {
	for (int i = 0; riff[i]; i++) {
		if (p[i] != riff[i])return false;
	}
	return true;
}

signed waveload(const char *fn, signed *audio_length, double **audio_buffer) {
	FILE *f;
	char name[4];
	long temp;
	if ((f = fopen(fn, "rb")) == NULL) {
		return 1;
	}
	fread(name, 4, 1, f);
	if (!waveriff("RIFF", name)) {
		fclose(f);
		return 2;
	}
	fread(&temp, 4, 1, f);
	fread(name, 4, 1, f);
	if (!waveriff("WAVE", name)) {
		fclose(f);
		return 3;
	}
	short format_num, format_chnum, format_blockalign, format_bitssample;
	long format_samplerate, format_byterate;
	while (true) {
		char chunk_name[4];
		long chunk_size;
		fread(chunk_name, 4, 1, f);
		fread(&chunk_size, 4, 1, f);
		if (feof(f)) {
			break;
		} else if (waveriff("fmt ", chunk_name)) {
			fread(&format_num, 2, 1, f);
			fread(&format_chnum, 2, 1, f);
			fread(&format_samplerate, 4, 1, f);
			fread(&format_byterate, 4, 1, f);
			fread(&format_blockalign, 2, 1, f);
			fread(&format_bitssample, 2, 1, f);
		} else if (waveriff("data", chunk_name)) {
			*audio_buffer = new double[*audio_length = chunk_size / format_blockalign];
			if (format_bitssample == 8)WAVEFILE_LOAD(f, char, format_chnum, *audio_length, *audio_buffer);
			if (format_bitssample == 16)WAVEFILE_LOAD(f, short, format_chnum, *audio_length, *audio_buffer);
			if (format_bitssample == 32)WAVEFILE_LOAD(f, long, format_chnum, *audio_length, *audio_buffer);
		} else {
			fseek(f, chunk_size, SEEK_CUR);
		}
	}
	fclose(f);
	return 0;
}

signed wavesave(const char *fn, signed audio_length, double *audio_buffer) {
	FILE *f;
	char name[4];
	long temp;
	fpos_t fposlist[3];
	if ((f = fopen(fn, "wb")) == NULL) {
		return 1;
	}
	fwrite("RIFF", 4, 1, f);
	fgetpos(f, &fposlist[0]);
	fwrite(&temp, 4, 1, f);
	fwrite("WAVE", 4, 1, f);
	short format_num = 1, format_chnum = 1, format_blockalign, format_bitssample=16;
	long format_samplerate = 44100, format_byterate;
	temp = 16;
	format_blockalign=format_chnum*format_bitssample/8;
	format_byterate=format_samplerate*format_blockalign;
	fwrite("fmt ", 4, 1, f);
	fwrite(&temp, 4, 1, f);
	fwrite(&format_num, 2, 1, f);
	fwrite(&format_chnum, 2, 1, f);
	fwrite(&format_samplerate, 4, 1, f);
	fwrite(&format_byterate, 4, 1, f);
	fwrite(&format_blockalign, 2, 1, f);
	fwrite(&format_bitssample, 2, 1, f);
	fwrite("data", 4, 1, f);
	fgetpos(f, &fposlist[1]);
	fwrite(&temp, 4, 1, f);
	if (format_bitssample == 8)WAVEFILE_SAVE(f, char, format_chnum, audio_length, audio_buffer);
	if (format_bitssample == 16)WAVEFILE_SAVE(f, short, format_chnum, audio_length, audio_buffer);
	if (format_bitssample == 32)WAVEFILE_SAVE(f, long, format_chnum, audio_length, audio_buffer);
	fgetpos(f, &fposlist[2]);
	fsetpos(f, &fposlist[1]);
	temp = fposlist[2] - fposlist[1];
	fwrite(&temp, 4, 1, f);
	fsetpos(f, &fposlist[0]);
	temp = fposlist[2] - fposlist[0];
	fwrite(&temp, 4, 1, f);
	fclose(f);
	return 0;
}

#undef  WAVEFILE_LOAD
#undef  WAVEFILE_SAVE
#endif //TEST_WAVEFILE_H
