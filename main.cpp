#include <iostream>
#include "wavefile.h"
int main() {
    double*p;
    long l;
    waveread("../test.wav",&l,&p);
    for(int i=0;i<l;i++){
        printf("%f\n",p[i]);
    }
    std::cout << "Hello, World!" << std::endl;
    return 0;
}