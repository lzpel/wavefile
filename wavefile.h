//
// Created by misumi3104 on 2019/06/19.
//

#ifndef TEST_WAVEFILE_H
#define TEST_WAVEFILE_H
#define WAVEFILE_LOAD(FP, T, CH, L, D) for(int i=0;i<(L);i++){T m,s=0;if(CH==1){fread(&s,sizeof(T),1,FP);}else{for(int j=0;j<CH;j++){fread(&m,sizeof(T),1,FP);s+=m;}s/=CH;}(D)[i]=s;}

#include<stdio.h>

bool waveriff(const char *riff, char *p) {
	for (int i = 0; riff[i]; i++) {
		if (p[i] != riff[i])return false;
	}
	return true;
}

signed waveread(const char *fn,long*audio_length,double**audio_buffer) {
	FILE *f;
	char name[4];
	long temp;
	if ((f=fopen(fn, "rb"))==NULL) {
		return 1;
	}
	fread(name, 4, 1, f);
	if (!waveriff("RIFF", name)) {
		fclose(f);
		return 2;
	}
	fread(&temp,4,1,f);
	fread(name, 4, 1, f);
	if (!waveriff("WAVE", name)) {
		fclose(f);
		return 3;
	}
	short   format_num,format_chnum,format_blockalign,format_bitssample;
	long    format_samplerate,format_byterate;
	while(true){
		char    chunk_name[4];
		long    chunk_size;
		fread(chunk_name,4,1,f);
		fread(&chunk_size,4,1,f);
		if(feof(f)) {
			break;
		}else if(waveriff("fmt ",chunk_name)){
			fread(&format_num,2,1,f);
			fread(&format_chnum,2,1,f);
			fread(&format_samplerate,4,1,f);
			fread(&format_byterate,4,1,f);
			fread(&format_blockalign,2,1,f);
			fread(&format_bitssample,2,1,f);
		}else if(waveriff("data",chunk_name)){
			*audio_buffer=new double[*audio_length=chunk_size/format_blockalign];
			switch(format_blockalign/format_chnum){
				case 1:
					WAVEFILE_LOAD(f,char,format_chnum,*audio_length,*audio_buffer);
					break;
				case 2:
					WAVEFILE_LOAD(f,short,format_chnum,*audio_length,*audio_buffer);
					break;
				case 4:
					WAVEFILE_LOAD(f,long,format_chnum,*audio_length,*audio_buffer);
					break;
			}
		}else{
			fseek(f,chunk_size,SEEK_CUR);
		}
	}
}

#undef  WAVEFILE_LOAD
#endif //TEST_WAVEFILE_H
