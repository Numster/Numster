#include "HoNucPos.h"
#include "HoTagClust.h"

#include <iostream>
#include <fstream>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>

extern double NUC_NOISE_FACTOR;

using namespace std;

namespace ho {
bool cmp_tag_chromosome(Tag tag1,Tag tag2) {
  return tag1._chromosome<tag2._chromosome;
}

vector<Tag>* read_tag(const char* f_name) {
  ifstream f_input;
  f_input.open(f_name);
  if(!f_input.is_open()) {
    cerr<<"file open error: "<<f_name<<endl;
    exit(1);
  }

  unsigned short chromosome;
  unsigned int position;
  char strand;
  unsigned short read;
  char buf[256];
  vector<Tag>* tag_list=new vector<Tag>();
  int count=0;

  f_input.getline(buf,256); //skip the head line
  while(f_input.getline(buf,256)) {
    if(sscanf(buf,"%hu\t%u\t%c\t%hu",&chromosome,&position,&strand,&read)!=4)
      continue;
    Tag tag={chromosome,position,strand,read,0.0};
    tag_list->push_back(tag);
    ++count;
  }
  f_input.close();
  //cout<<"("<<count<<"==";

  return tag_list;
}

void write_nucleosome(vector<Nuc>* nuc_list,const char* f_name) {
  ofstream f_output;
  f_output.open(f_name);
  if(!f_output.is_open()) {
    cerr<<"file open error: "<<f_name<<endl;
    exit(1);
  }

  //f_output<<"Chr\tCenter\tSize"<<endl;
  f_output<<"Chr\tStart\tEnd\tOccupancy"<<endl;
  for(vector<Nuc>::iterator nuc_iter=nuc_list->begin();nuc_iter!=nuc_list->end();++nuc_iter) {
    f_output<<nuc_iter->_chromosome<<"\t";
    f_output<<nuc_iter->_start<<"\t";
    f_output<<nuc_iter->_end<<"\t";
    f_output<<(nuc_iter->_f_occupancy+nuc_iter->_r_occupancy)<<endl;
  }
  f_output.close();
}

void* position_nucleosome(void* data) {
  vector<Tag>* tag_list=(vector<Tag>*) data;
  vector<Nuc>* nuc_list=new vector<Nuc>();
  sort(tag_list->begin(),tag_list->end(),cmp_tag_chromosome);
  int chr_start=tag_list->at(0)._chromosome;
  int chr_end=tag_list->back()._chromosome;

  vector<Tag>::iterator tag_iter=tag_list->begin();
  for(int i=chr_start;i<=chr_end;++i) {
    cout<<"chromosome "<<i<<" started"<<endl;
    vector<Tag> temp_tag_list;

    for(;tag_iter!=tag_list->end();++tag_iter) {
      if(tag_iter->_chromosome==i)
	temp_tag_list.push_back(*tag_iter);
      else
	break;
    }
    if(temp_tag_list.size()==0) {
      cerr<<"!!!! no tags for chromosome "<<i<<" !!!!"<<endl;
      continue;
    }

    TagClust tag_cluster(&temp_tag_list,temp_tag_list.size());
    tag_cluster.determine_center(true); //forward strand cluster initialization
    tag_cluster.determine_center(false); //reverse strand cluster initialization
    //cout<<"Processing 1"<<endl;
    tag_cluster.create_local_cluster(true); //forward strand clustering
    tag_cluster.create_local_cluster(false); //reverse stran clustering
    //cout<<"Processing 2"<<endl;
    tag_cluster.create_global_cluster(); //minimizing f_r cluster difference
    //cout<<"Processing 3"<<endl;
    NucPos nuc_position(i,tag_cluster.get_f_cluster_list(),tag_cluster.get_r_cluster_list(),temp_tag_list.size());
    nuc_position.position_nucleosome();

    nuc_list->insert(nuc_list->end(),nuc_position.get_nucleosome()->begin(),nuc_position.get_nucleosome()->end());

    cout<<"chromosome "<<i<<" ended with "<<nuc_position.get_nuc_count()<<" nucleosomes"<<endl;
  }

  pthread_exit((void*)nuc_list);
}

int cal_occ_quantile(vector<int>* occ,double q) {
  sort(occ->begin(),occ->end());

  int total=0;
  for(int i=0;i<occ->size();++i)
    total+=(*occ)[i];
  double q_value=total*q;
  int value=0;

  total=0;
  for(int i=0;i<occ->size();++i) {
    total+=(*occ)[i];
    if(total>=q_value) {
      value=(*occ)[i];
      break;
    }
  }
  return value;
}

void clean_noisy_nucleosome(vector<Nuc>* nuc_list) {
  vector<int> occ;
  for(vector<Nuc>::iterator nuc_iter=nuc_list->begin();nuc_iter!=nuc_list->end();++nuc_iter)
    occ.push_back(nuc_iter->_f_occupancy+nuc_iter->_r_occupancy);

  double q1=cal_occ_quantile(&occ,0.25);
  double q3=cal_occ_quantile(&occ,0.75);

  double max_cutoff=q3-q1;
  double cutoff=q1-NUC_NOISE_FACTOR*(q3-q1);

  if(cutoff>max_cutoff)
    cutoff=max_cutoff;

  int removed_nuc_count=0;
  int original_nuc_count=nuc_list->size();
  for(vector<Nuc>::iterator nuc_iter=nuc_list->begin();nuc_iter!=nuc_list->end();) {
    if((nuc_iter->_f_occupancy+nuc_iter->_r_occupancy)<=cutoff) {
      nuc_iter=nuc_list->erase(nuc_iter);
      ++removed_nuc_count;
    }
    else
      ++nuc_iter;
  }
}

void process(const char* infile_name,const char* outfile_name,int td_size) {
  time_t t0,t1;
  t0=time(NULL);

  vector<Tag>* tag_list=read_tag(infile_name);
  //cout<<tag_list->size()<<")"<<" locations"<<endl<<endl;
  sort(tag_list->begin(),tag_list->end(),cmp_tag_chromosome);
  int chr_size=tag_list->back()._chromosome;

  if(td_size>chr_size)
    td_size=chr_size;

  int *chr_count=new int[td_size];
  for(int i=0;i<td_size;++i)
    chr_count[i]=0;

  int chr=1;
  bool flag=false;
  while(true) {
    for(int i=0;i<td_size;++i) {
      if(chr>chr_size) {
	flag=true;
	break;
      }
      chr_count[i]+=1;
      ++chr;
    }
    if(flag==true)
      break;
  }

  pthread_t *tds=new pthread_t[td_size];
  int **chrs=new int*[td_size];

  chr=1;
  for(int i=0;i<td_size;++i) {
    chrs[i]=new int[2];
    chrs[i][0]=chr;
    for(int j=0;j<chr_count[i];++j)
      ++chr;
    chrs[i][1]=chr-1;
  }

  vector<Tag> **td_tag_list=new vector<Tag>*[td_size];
  for(int i=0;i<td_size;++i)
    td_tag_list[i]=new vector<Tag>;

  vector<Tag>::iterator tag_iter;
  for(tag_iter=tag_list->begin();tag_iter!=tag_list->end();++tag_iter) {
    int temp_chr=tag_iter->_chromosome;
    for(int i=0;i<td_size;++i) {
      if(temp_chr>=chrs[i][0] && temp_chr<=chrs[i][1]) {
	td_tag_list[i]->push_back(*tag_iter);
	break;
      }
    }
  }

  delete [] chr_count;
  chr_count=0;
  for(int i=0;i<td_size;++i) {
    delete [] chrs[i];
    chrs[i]=0;
  }
  delete [] chrs;
  chrs=0;
  delete tag_list;
  tag_list=0;

  for(int i=0;i<td_size;++i)
    pthread_create(&tds[i],NULL,position_nucleosome,(void*)td_tag_list[i]);

  vector<Nuc>* nuc_list=new vector<Nuc>();
  vector<Nuc>* temp_nuc_list=new vector<Nuc>();
  for(int i=0;i<td_size;++i) {
    pthread_join(tds[i],(void**)&temp_nuc_list);
    nuc_list->insert(nuc_list->end(),temp_nuc_list->begin(),temp_nuc_list->end());
  }

  for(int i=0;i<td_size;++i) {
    delete td_tag_list[i];
    td_tag_list[i]=0;
  }
  delete [] td_tag_list;
  td_tag_list=0;

  delete [] tds;
  tds=0;

  clean_noisy_nucleosome(nuc_list);

  cout<<endl<<"** total nucleosome count: "<<nuc_list->size()<<endl;

  write_nucleosome(nuc_list,outfile_name);

  //t1=time(NULL);
  //cout<<"** time: "<<t1-t0<<" seconds"<<endl<<endl;

  delete temp_nuc_list;
  temp_nuc_list=0;
  delete nuc_list;
  nuc_list=0;
}
}
