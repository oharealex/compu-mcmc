#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

void average(float*, float*);

int main(){
  float mean;
  float var;
  average(&mean,&var);
  cout<<"mean "<<" "<<mean<<endl;
  cout<<"variance "<<" "<<var<<endl;

  return 0;
}

void average(float *mean, float *var) {
  char str[16];
  int i;
  float x;
  i = 0;
  *var=0;
  *mean=0;
  while(fgets(str, sizeof str, stdin) != NULL) {  // reads newline too
    x=atof(str);
    *mean+=x;
    *var=*var+x*x;
    i++;
  }
  *var=(*var-pow(*mean,2)/((float)i))/((float)(i-1));
  *mean=*mean/((float)i);
  return ;
}

