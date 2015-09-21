#include "HoNucPos.h"
#include <math.h>

extern double OVERLAP_RATE;

namespace ho {
extern bool cmp_cluster_center(Cluster,Cluster);
extern bool cmp_cluster_strength(Cluster,Cluster);
extern bool cmp_cluster_id(Cluster,Cluster);
extern bool cmp_tag_start(Tag tag1,Tag tag2);

bool cmp_nuc_start(Nuc nuc1,Nuc nuc2) {
  return nuc1._start<nuc2._start;
}

bool cmp_nuc_occupancy(Nuc nuc1,Nuc nuc2) {
  return nuc1._occupancy<nuc2._occupancy;
}

double cal_occ(vector<Tag>* tag_list,unsigned int location) {
  sort(tag_list->begin(),tag_list->end(),cmp_tag_start);
  double total=0.0;
  int count=0;
  for(int i=0;i<tag_list->size();++i) {
    if((*tag_list)[i]._start==location) {
      total+=(*tag_list)[i]._read;
      ++count;
      if(i>0) {total+=(*tag_list)[i-1]._read;++count;}
      if(i<tag_list->size()-1) {total+=(*tag_list)[i+1]._read;++count;}
      break;
    }
  }

  if(count==0) {
    cerr<<"occupancy error"<<endl;
    exit(1);
  }

  return total/count;
}

NucPos::NucPos(unsigned short chromosome,vector<Cluster>* f_cluster_list,vector<Cluster>* r_cluster_list,double size) {
  sort(f_cluster_list->begin(),f_cluster_list->end(),cmp_cluster_center);
  sort(r_cluster_list->begin(),r_cluster_list->end(),cmp_cluster_center);
  vector<Cluster> f_tmp_list(*f_cluster_list);
  vector<Cluster> r_tmp_list(*r_cluster_list);
  for(int i=0;i<f_tmp_list.size();++i)
    sort(f_tmp_list[i]._point.begin(),f_tmp_list[i]._point.end(),cmp_tag_start);
  for(int i=0;i<r_tmp_list.size();++i)
    sort(r_tmp_list[i]._point.begin(),r_tmp_list[i]._point.end(),cmp_tag_start);

  sort(f_cluster_list->begin(),f_cluster_list->end(),cmp_cluster_strength);
  sort(r_cluster_list->begin(),r_cluster_list->end(),cmp_cluster_strength);

  double f_occ=0.0;
  double r_occ=0.0;

  _nuc_list=new vector<Nuc>();

  for(vector<Cluster>::reverse_iterator f_cluster_iter=f_cluster_list->rbegin();f_cluster_iter!=f_cluster_list->rend();++f_cluster_iter) {
    vector<unsigned int> value;
    for(int i=0;i<f_cluster_iter->_point.size();++i)
      for(int j=0;j<f_cluster_iter->_point[i]._read;++j)
	value.push_back(f_cluster_iter->_point[i]._start);
    double f_std=cal_center_std(&value,f_cluster_iter->_center);
    value.clear();
    f_occ=cal_occ(&(f_cluster_iter->_point),f_cluster_iter->_peak); //mode

    if(f_cluster_iter->_id==1) {
      if(f_cluster_iter->_peak<NUC_HALF_SIZE-1)
	continue;
      Nuc nuc={chromosome,f_cluster_iter->_peak+1-NUC_HALF_SIZE,f_cluster_iter->_peak+1+NUC_HALF_SIZE,f_cluster_iter->_strength,0,f_cluster_iter->_strength,f_occ,false,false,f_std,f_occ,0};
      _nuc_list->push_back(nuc);
    }
    else if(f_cluster_iter->_id==2) {
      if(f_cluster_iter->_peak<NUC_HALF_SIZE+6)
	continue;
      Nuc nuc={chromosome,f_cluster_iter->_peak-6-NUC_HALF_SIZE,f_cluster_iter->_peak-6+NUC_HALF_SIZE,f_cluster_iter->_strength,0,f_cluster_iter->_strength,f_occ,false,false,f_std,f_occ,0};
      _nuc_list->push_back(nuc);
    }
    else if(f_cluster_iter->_id==3) {
      if(f_cluster_iter->_peak<NUC_HALF_SIZE+6)
	continue;

      double f_occ1=0;
      double f_occ2=0;
      for(int i=0;i<f_tmp_list.size();++i) {
	if(f_cluster_iter->_peak-10>f_tmp_list[i]._end)
	  continue;
	if(f_cluster_iter->_peak+10<f_tmp_list[i]._start)
	  break;
	for(int j=0;j<f_tmp_list[i]._point.size();++j) {
	  if(f_tmp_list[i]._point[j]._start>=f_cluster_iter->_peak-10 && f_tmp_list[i]._point[j]._start<=f_cluster_iter->_peak)
	    f_occ1+=f_tmp_list[i]._point[j]._read;
	  if(f_tmp_list[i]._point[j]._start>=f_cluster_iter->_peak && f_tmp_list[i]._point[j]._start<=f_cluster_iter->_peak+10)
	    f_occ2+=f_tmp_list[i]._point[j]._read;
	}
      }
      f_occ1/=10;
      f_occ2/=10;

      Nuc nuc1={chromosome,f_cluster_iter->_peak-6-NUC_HALF_SIZE,f_cluster_iter->_peak-6+NUC_HALF_SIZE,f_cluster_iter->_strength,0,f_cluster_iter->_strength,f_occ1,false,false,f_std,f_occ1,0};
      Nuc nuc2={chromosome,f_cluster_iter->_peak+1-NUC_HALF_SIZE,f_cluster_iter->_peak+1+NUC_HALF_SIZE,f_cluster_iter->_strength,0,f_cluster_iter->_strength,f_occ2,false,false,f_std,f_occ2,0};
      _nuc_list->push_back(nuc1);
      _nuc_list->push_back(nuc2);
    }
  }

  for(vector<Cluster>::reverse_iterator r_cluster_iter=r_cluster_list->rbegin();r_cluster_iter!=r_cluster_list->rend();++r_cluster_iter) {
    vector<unsigned int> value;
    for(int i=0;i<r_cluster_iter->_point.size();++i)
      for(int j=0;j<r_cluster_iter->_point[i]._read;++j)
	value.push_back(r_cluster_iter->_point[i]._start);
    double r_std=cal_center_std(&value,r_cluster_iter->_center);
    value.clear();
    r_occ=cal_occ(&(r_cluster_iter->_point),r_cluster_iter->_peak); //mode

    if(r_cluster_iter->_id==1) {
      if(r_cluster_iter->_peak<NUC_HALF_SIZE+1)
	continue;
      Nuc nuc={chromosome,r_cluster_iter->_peak-1-NUC_HALF_SIZE,r_cluster_iter->_peak-1+NUC_HALF_SIZE,0,r_cluster_iter->_strength,r_cluster_iter->_strength,r_occ,false,false,r_std,0,r_occ};
      _nuc_list->push_back(nuc);
    }
    else if(r_cluster_iter->_id==2) {
      if(r_cluster_iter->_peak<NUC_HALF_SIZE-6)
	continue;
      Nuc nuc={chromosome,r_cluster_iter->_peak+6-NUC_HALF_SIZE,r_cluster_iter->_peak+6+NUC_HALF_SIZE,0,r_cluster_iter->_strength,r_cluster_iter->_strength,r_occ,false,false,r_std,0,r_occ};
      _nuc_list->push_back(nuc);
    }
    else if(r_cluster_iter->_id==3) {
      if(r_cluster_iter->_peak<NUC_HALF_SIZE+1)
	continue;

      double r_occ1=0;
      double r_occ2=0;
      for(int i=0;i<r_tmp_list.size();++i) {
	if(r_cluster_iter->_peak-10>r_tmp_list[i]._end)
	  continue;
	if(r_cluster_iter->_peak+10<r_tmp_list[i]._start)
	  break;
	for(int j=0;j<r_tmp_list[i]._point.size();++j) {
	  if(r_tmp_list[i]._point[j]._start>=r_cluster_iter->_peak-10 && r_tmp_list[i]._point[j]._start<=r_cluster_iter->_peak)
	    r_occ1+=r_tmp_list[i]._point[j]._read;
	  if(r_tmp_list[i]._point[j]._start>=r_cluster_iter->_peak && r_tmp_list[i]._point[j]._start<=r_cluster_iter->_peak+10)
	    r_occ2+=r_tmp_list[i]._point[j]._read;
	}
      }
      r_occ1/=10;
      r_occ2/=10;

      Nuc nuc1={chromosome,r_cluster_iter->_peak-1-NUC_HALF_SIZE,r_cluster_iter->_peak-1+NUC_HALF_SIZE,0,r_cluster_iter->_strength,r_cluster_iter->_strength,r_occ1,false,false,r_std,0,r_occ1};
      Nuc nuc2={chromosome,r_cluster_iter->_peak+6-NUC_HALF_SIZE,r_cluster_iter->_peak+6+NUC_HALF_SIZE,0,r_cluster_iter->_strength,r_cluster_iter->_strength,r_occ2,false,false,r_std,0,r_occ2};
      _nuc_list->push_back(nuc1);
      _nuc_list->push_back(nuc2);
    }
  }
}

NucPos::~NucPos() {
  delete _nuc_list;
  _nuc_list=0;
}

double NucPos::cal_center_std(vector<unsigned int>* value,unsigned int center) {
  double dev=0.0;
  for(int i=0;i<value->size();++i)
    dev+=((*value)[i]-center)*((*value)[i]-center);
  return sqrt(dev/value->size());
}

void NucPos::position_nucleosome() {
  sort(_nuc_list->begin(),_nuc_list->end(),cmp_nuc_start);

  vector<Nuc>* nuc_list=new vector<Nuc>();

  unsigned int upper_limit=0;
  vector<Nuc> temp_nuc_list;
  vector<Nuc> new_nuc_list;
  vector<Nuc>::iterator nuc_iter=_nuc_list->begin();
  temp_nuc_list.push_back(*nuc_iter);
  upper_limit=nuc_iter->_end;

  for(++nuc_iter;nuc_iter!=_nuc_list->end();++nuc_iter) {
    if(nuc_iter==(_nuc_list->end()-1)) {
      if(nuc_iter->_start<upper_limit)
	temp_nuc_list.push_back(*nuc_iter);

      define_maximal_nucleosome(&temp_nuc_list,&new_nuc_list);

      nuc_list->insert(nuc_list->end(),new_nuc_list.begin(),new_nuc_list.end());
      temp_nuc_list.clear();
      new_nuc_list.clear();

      if(nuc_iter->_start>=upper_limit) {
	temp_nuc_list.push_back(*nuc_iter);

	define_maximal_nucleosome(&temp_nuc_list,&new_nuc_list);

	nuc_list->insert(nuc_list->end(),new_nuc_list.begin(),new_nuc_list.end());
	temp_nuc_list.clear();
	new_nuc_list.clear();
      }
    }
    else if(nuc_iter->_start<upper_limit) {
      temp_nuc_list.push_back(*nuc_iter);
      upper_limit=nuc_iter->_end;
    }
    else {
      define_maximal_nucleosome(&temp_nuc_list,&new_nuc_list);

      nuc_list->insert(nuc_list->end(),new_nuc_list.begin(),new_nuc_list.end());
      temp_nuc_list.clear();
      new_nuc_list.clear();
      temp_nuc_list.push_back(*nuc_iter);
      upper_limit=nuc_iter->_end;
    }
  }

  double original_nuc_size=_nuc_list->size();
  delete _nuc_list;
  _nuc_list=nuc_list;
  _fuzzy_rate=(original_nuc_size-_nuc_list->size())/original_nuc_size;
}

void NucPos::define_maximal_nucleosome(vector<Nuc>* temp_nuc_list,vector<Nuc>* new_nuc_list) { //get the max nuc
  sort(temp_nuc_list->begin(),temp_nuc_list->end(),cmp_nuc_occupancy);

  unsigned short nuc_size=0;
  unsigned int overlap_size=0;
  int overlap_index=-1;

  for(vector<Nuc>::reverse_iterator temp_nuc_iter=temp_nuc_list->rbegin();temp_nuc_iter!=temp_nuc_list->rend();++temp_nuc_iter) {
    overlap_index=-1;
    for(int i=0;i<new_nuc_list->size();++i) {
      if((temp_nuc_iter->_end-temp_nuc_iter->_start)<((*new_nuc_list)[i]._end-(*new_nuc_list)[i]._start))
	nuc_size=temp_nuc_iter->_end-temp_nuc_iter->_start;
      else
	nuc_size=(*new_nuc_list)[i]._end-(*new_nuc_list)[i]._start;

      if((*new_nuc_list)[i]._start>=temp_nuc_iter->_start && (*new_nuc_list)[i]._start<temp_nuc_iter->_end) { //left-from overlap
	overlap_size=temp_nuc_iter->_end-(*new_nuc_list)[i]._start;
	if(double(overlap_size)/double(nuc_size)>OVERLAP_RATE) {
	  overlap_index=i;
	  break;
	}
      }
      if(temp_nuc_iter->_start>=(*new_nuc_list)[i]._start && temp_nuc_iter->_start<(*new_nuc_list)[i]._end) { //right-from overlap
	overlap_size=(*new_nuc_list)[i]._end-temp_nuc_iter->_start;
	if(double(overlap_size)/double(nuc_size)>OVERLAP_RATE) {
	  overlap_index=i;
	  break;
	}
      }
    }

    if(overlap_index==-1) {
      new_nuc_list->push_back(*temp_nuc_iter);
    }
    else {
      //(*new_nuc_list)[overlap_index]._occupancy+=temp_nuc_iter->_occupancy;
      (*new_nuc_list)[overlap_index]._overlap=true;
    }
  }
}

vector<Nuc>* NucPos::get_nucleosome() {
  sort(_nuc_list->begin(),_nuc_list->end(),cmp_nuc_start);
  return _nuc_list;
}

int NucPos::get_nuc_count() {
  return _nuc_list->size();
}

double NucPos::get_fuzzy_rate() {
  return _fuzzy_rate;
}
}
