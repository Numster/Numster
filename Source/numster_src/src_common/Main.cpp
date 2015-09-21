#include <iostream>
#include <fstream>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>

using namespace std;

double LOCAL_NOISE_FACTOR=0.5;
double NUC_NOISE_FACTOR=0.5;
double GLOBAL_NOISE_FACTOR=0.25;
double OVERLAP_RATE=0.3;
double MERGE_DISTANCE=0.5;
double EUCLIDEAN_DISTANCE=0.1;
int MN_INIT_CLUSTER_SIZE=35;
int HO_INIT_CLUSTER_SIZE=3;
int SEARCH_LOWER_BOUND=140;
int SEARCH_UPPER_BOUND=180;
int DEFAULT_NUC_LOWER_BOUND=80;
int DEFAULT_NUC_UPPER_BOUND=200;

namespace mn {
  extern void process(const char*,const char*,int);
}

namespace ho {
  extern void process(const char*,const char*,int);
}

void message() {
  cout<<"usage"<<endl;
  cout<<"  numster arguments"<<endl;
  cout<<"mandatory arguments"<<endl;
  cout<<"  -i: path of input file"<<endl;
  cout<<"  -m: 1 for mnase data. 2 for hydroxyl radical data"<<endl;
  cout<<"  -o: path of output file"<<endl;
  cout<<"optional arguments"<<endl;
  cout<<"  -a: distance threshold for cluster aggregation. default 0.5"<<endl;
  cout<<"  -b: interquartile range value for noisy cluster removal. default 0.5"<<endl;
  cout<<"  -d: euclidean distance threshold for conflicting templates. default 0.1"<<endl;
  cout<<"  -e: distance to mode tag in initial cluster. default 35 for -m 1. default 3 for -m 2"<<endl;
  cout<<"  -n: interquartile range value for noisy nucleosome removal. default 0.5"<<endl;
  cout<<"  -t: thread size. default 4"<<endl;
  cout<<"  -l: left boundary of nucleosomal dna range. default 80 for -m 1. default 140 for -m 2"<<endl;
  cout<<"  -r: right boundary of nucleosomal dna range. default 200 for -m 1. default 180 for -m 2"<<endl;
  cout<<"  -v: maximum overlapping rate of nucleosomes. default 0.3"<<endl;
  cout<<"usage examples"<<endl;
  cout<<"  numster -m 1 -i test/mn_tag.txt -o test/mn_nuc.txt"<<endl;
  cout<<"  numster -m 2 -t 1 -i test/ho_tag.txt -o test/ho_nuc.txt"<<endl;
}

int main(int argc,const char* argv[]) {
  const char* option=NULL;
  const char* infile_name=NULL;
  const char* outfile_name=NULL;
  int td_size=4;
  int method=0;
  int cluster_size=0;
  int lower_bound=-1;
  int upper_bound=-1;

  for(int i=1;i<argc;++i) {
    if(argv[i][0]=='-') {
      option=argv[i++]+1;
      if(strcmp(option,"m")==0)
	method=atoi(argv[i]);
      else if(strcmp(option,"b")==0)
	LOCAL_NOISE_FACTOR=atof(argv[i]);
      else if(strcmp(option,"n")==0)
	NUC_NOISE_FACTOR=atof(argv[i]);
      else if(strcmp(option,"t")==0)
	td_size=atoi(argv[i]);
      else if(strcmp(option,"e")==0)
	cluster_size=atoi(argv[i]);
      else if(strcmp(option,"l")==0)
	lower_bound=atoi(argv[i]);
      else if(strcmp(option,"r")==0)
	upper_bound=atoi(argv[i]);
      else if(strcmp(option,"a")==0)
	MERGE_DISTANCE=atof(argv[i]);
      else if(strcmp(option,"d")==0)
	EUCLIDEAN_DISTANCE=atof(argv[i]);
      else if(strcmp(option,"v")==0)
	OVERLAP_RATE=atof(argv[i]);
      else if(strcmp(option,"i")==0)
	infile_name=argv[i];
      else if(strcmp(option,"o")==0)
	outfile_name=argv[i];
      else {
	cerr<<"option error"<<endl;
	message();
	exit(1);
      }
    }
    else {
      cerr<<"option error"<<endl;
      message();
      exit(1);
    }
  }

  if(infile_name==NULL || strlen(infile_name)==0) {
    cerr<<"no input file name"<<endl;
    message();
    exit(1);
  }
  if(outfile_name==NULL || strlen(outfile_name)==0) {
    cerr<<"no output file name"<<endl;
    message();
    exit(1);
  }

  if(method==1) {
    if(cluster_size>0)
      MN_INIT_CLUSTER_SIZE=cluster_size;
    if(lower_bound>-1)
      DEFAULT_NUC_LOWER_BOUND=lower_bound;
    if(upper_bound>-1)
      DEFAULT_NUC_UPPER_BOUND=upper_bound;

    mn::process(infile_name,outfile_name,td_size);
  }
  else if(method==2) {
    if(cluster_size>0)
      HO_INIT_CLUSTER_SIZE=cluster_size;
    if(lower_bound>-1)
      SEARCH_LOWER_BOUND=lower_bound;
    if(upper_bound>-1)
      SEARCH_UPPER_BOUND=upper_bound;

    ho::process(infile_name,outfile_name,td_size);
  }
  else {
    cerr<<"no proper method"<<endl;
    message();
    exit(1);
  }

  return 0;
}
