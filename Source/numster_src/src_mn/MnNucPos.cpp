#include "MnNucPos.h"
#include <math.h>

extern double OVERLAP_RATE;

namespace mn {
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
  sort(f_cluster_list->begin(),f_cluster_list->end(),cmp_cluster_strength);
  sort(r_cluster_list->begin(),r_cluster_list->end(),cmp_cluster_strength);

  int id=0;
  unsigned int lower_limit=0;
  unsigned int upper_limit=0;
  double f_occ=0.0;
  double r_occ=0.0;

  for(vector<Cluster>::iterator f_cluster_iter=f_cluster_list->begin();f_cluster_iter!=f_cluster_list->end();++f_cluster_iter)
    f_cluster_iter->_id=-1;
  for(vector<Cluster>::iterator r_cluster_iter=r_cluster_list->begin();r_cluster_iter!=r_cluster_list->end();++r_cluster_iter)
    r_cluster_iter->_id=-1;

  _nuc_list=new vector<Nuc>();

  for(vector<Cluster>::reverse_iterator f_cluster_iter=f_cluster_list->rbegin();f_cluster_iter!=f_cluster_list->rend();++f_cluster_iter) {
    lower_limit=f_cluster_iter->_center+NUC_MIN_WIDTH;
    upper_limit=f_cluster_iter->_center+NUC_MAX_WIDTH;

    for(vector<Cluster>::reverse_iterator r_cluster_iter=r_cluster_list->rbegin();r_cluster_iter!=r_cluster_list->rend();++r_cluster_iter) {
      if(r_cluster_iter->_center>=lower_limit && r_cluster_iter->_center<=upper_limit) {
	if(r_cluster_iter->_id!=-1)
	  continue;
	f_cluster_iter->_id=id;
	r_cluster_iter->_id=id;

	vector<unsigned int> value;
	for(int i=0;i<f_cluster_iter->_point.size();++i)
	  for(int j=0;j<f_cluster_iter->_point[i]._read;++j)
	    value.push_back(f_cluster_iter->_point[i]._start);
	double f_std=cal_center_std(&value,f_cluster_iter->_center);
	value.clear();
	for(int i=0;i<r_cluster_iter->_point.size();++i)
	  for(int j=0;j<r_cluster_iter->_point[i]._read;++j)
	    value.push_back(r_cluster_iter->_point[i]._start);
	double r_std=cal_center_std(&value,r_cluster_iter->_center);
	value.clear();

	//int height=f_cluster_iter->_max_read;
	//if(r_cluster_iter->_max_read>f_cluster_iter->_max_read)
	//height=r_cluster_iter->_max_read;
	//occ=sqrt(height+pow(f_cluster_iter->_strength+r_cluster_iter->_strength,2.0)/size); //original

	//f_occ=cal_occ(&(f_cluster_iter->_point),f_cluster_iter->_center); //median
	//r_occ=cal_occ(&(r_cluster_iter->_point),r_cluster_iter->_center); //median
	f_occ=cal_occ(&(f_cluster_iter->_point),f_cluster_iter->_peak); //mode
	r_occ=cal_occ(&(r_cluster_iter->_point),r_cluster_iter->_peak); //mode

	Nuc nuc={chromosome,f_cluster_iter->_center,r_cluster_iter->_center,f_cluster_iter->_strength,r_cluster_iter->_strength,
		 (f_cluster_iter->_strength+r_cluster_iter->_strength),(f_occ+r_occ)/2.0,false,true,(f_std+r_std)/2.0,f_occ,r_occ};
	_nuc_list->push_back(nuc);
	++id;
	break;
      }
    }
  }

  if(UNPAIRED_NUC==true) {
    for(vector<Cluster>::iterator f_cluster_iter=f_cluster_list->begin();f_cluster_iter!=f_cluster_list->end();++f_cluster_iter) {
      if(f_cluster_iter->_id!=-1)
	continue;

      unsigned int nuc_width=DEFAULT_NUC_WIDTH;
      if(DEFAULT_NUC_WIDTH<NUC_MIN_WIDTH)
	nuc_width=NUC_MIN_WIDTH;
      if(DEFAULT_NUC_WIDTH>NUC_MAX_WIDTH)
	nuc_width=NUC_MAX_WIDTH;

      vector<unsigned int> value;
      for(int i=0;i<f_cluster_iter->_point.size();++i)
	for(int j=0;j<f_cluster_iter->_point[i]._read;++j)
	  value.push_back(f_cluster_iter->_point[i]._start);
      double f_std=cal_center_std(&value,f_cluster_iter->_center);
      value.clear();

      //occ=sqrt(f_cluster_iter->_max_read+pow(f_cluster_iter->_strength,2.0)/size);
      //f_occ=cal_occ(&(f_cluster_iter->_point),f_cluster_iter->_center); //median
      f_occ=cal_occ(&(f_cluster_iter->_point),f_cluster_iter->_peak); //mode

      Nuc nuc={chromosome,f_cluster_iter->_center,f_cluster_iter->_center+nuc_width,f_cluster_iter->_strength,0,f_cluster_iter->_strength,f_occ/2.0,false,false,f_std,f_occ,0};
      _nuc_list->push_back(nuc);
    }

    for(vector<Cluster>::iterator r_cluster_iter=r_cluster_list->begin();r_cluster_iter!=r_cluster_list->end();++r_cluster_iter) {
      if(r_cluster_iter->_id!=-1)
	continue;

      unsigned int nuc_width=DEFAULT_NUC_WIDTH;
      if(DEFAULT_NUC_WIDTH<NUC_MIN_WIDTH)
	nuc_width=NUC_MIN_WIDTH;
      if(DEFAULT_NUC_WIDTH>NUC_MAX_WIDTH)
	nuc_width=NUC_MAX_WIDTH;

      vector<unsigned int> value;
      for(int i=0;i<r_cluster_iter->_point.size();++i)
	for(int j=0;j<r_cluster_iter->_point[i]._read;++j)
	  value.push_back(r_cluster_iter->_point[i]._start);
      double r_std=cal_center_std(&value,r_cluster_iter->_center);
      value.clear();

      //occ=sqrt(r_cluster_iter->_max_read+pow(r_cluster_iter->_strength,2.0)/size);
      //r_occ=cal_occ(&(r_cluster_iter->_point),r_cluster_iter->_center); //median
      r_occ=cal_occ(&(r_cluster_iter->_point),r_cluster_iter->_peak); //mode

      Nuc nuc={chromosome,r_cluster_iter->_center-nuc_width,r_cluster_iter->_center,0,r_cluster_iter->_strength,r_cluster_iter->_strength,r_occ/2.0,false,false,r_std,0,r_occ};
      _nuc_list->push_back(nuc);
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

      if(COMBINED_NUC==true)
	define_weighted_nucleosome(&temp_nuc_list,&new_nuc_list);
      else
	define_maximal_nucleosome(&temp_nuc_list,&new_nuc_list);

      nuc_list->insert(nuc_list->end(),new_nuc_list.begin(),new_nuc_list.end());
      temp_nuc_list.clear();
      new_nuc_list.clear();

      if(nuc_iter->_start>=upper_limit) {
	temp_nuc_list.push_back(*nuc_iter);

	if(COMBINED_NUC==true)
	  define_weighted_nucleosome(&temp_nuc_list,&new_nuc_list);
	else
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
      if(COMBINED_NUC==true)
	define_weighted_nucleosome(&temp_nuc_list,&new_nuc_list);
      else
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

void NucPos::define_weighted_nucleosome(vector<Nuc>* temp_nuc_list,vector<Nuc>* new_nuc_list) { //get the weighted center
  sort(temp_nuc_list->begin(),temp_nuc_list->end(),cmp_nuc_occupancy);

  unsigned short nuc_size=0;
  unsigned int overlap_size=0;
  unsigned int max_overlap_size=0;
  int overlap_index=-1;
  int overlap_type=-1;
  unsigned int center=0;
  unsigned int new_nuc_center=0;
  unsigned int temp_nuc_center=0;

  for(vector<Nuc>::reverse_iterator temp_nuc_iter=temp_nuc_list->rbegin();temp_nuc_iter!=temp_nuc_list->rend();++temp_nuc_iter) {
    max_overlap_size=0;
    overlap_index=-1;
    overlap_type=-1;
    for(int i=0;i<new_nuc_list->size();++i) {
      if((temp_nuc_iter->_end-temp_nuc_iter->_start)<((*new_nuc_list)[i]._end-(*new_nuc_list)[i]._start))
	nuc_size=temp_nuc_iter->_end-temp_nuc_iter->_start;
      else
	nuc_size=(*new_nuc_list)[i]._end-(*new_nuc_list)[i]._start;
 
      if((*new_nuc_list)[i]._start>=temp_nuc_iter->_start && (*new_nuc_list)[i]._start<temp_nuc_iter->_end) { //left-from overlap
	overlap_size=temp_nuc_iter->_end-(*new_nuc_list)[i]._start;
	if(overlap_size>max_overlap_size && (double(overlap_size)/double(nuc_size)>OVERLAP_RATE)) {
	  overlap_index=i;
	  max_overlap_size=overlap_size;
	  new_nuc_center=(*new_nuc_list)[i]._start+(unsigned int)(((*new_nuc_list)[i]._end-(*new_nuc_list)[i]._start)/2.0);
	  temp_nuc_center=temp_nuc_iter->_start+(unsigned int)((temp_nuc_iter->_end-temp_nuc_iter->_start)/2.0);

	  if(new_nuc_center>temp_nuc_center)
	    overlap_type=1;
	  else
	    overlap_type=2;
	}
      }
      if(temp_nuc_iter->_start>=(*new_nuc_list)[i]._start && temp_nuc_iter->_start<(*new_nuc_list)[i]._end) { //right-from overlap
	overlap_size=(*new_nuc_list)[i]._end-temp_nuc_iter->_start;
	if(overlap_size>max_overlap_size && ((double(overlap_size)/double(nuc_size))>OVERLAP_RATE)) {
	  overlap_index=i;
	  max_overlap_size=overlap_size;
	  new_nuc_center=(*new_nuc_list)[i]._start+(unsigned int)(((*new_nuc_list)[i]._end-(*new_nuc_list)[i]._start)/2.0);
	  temp_nuc_center=temp_nuc_iter->_start+(unsigned int)((temp_nuc_iter->_end-temp_nuc_iter->_start)/2.0);

	  if(new_nuc_center<temp_nuc_center)
	    overlap_type=3;
	  else
	    overlap_type=4;
	}
      }
    }

    if(overlap_index==-1) {
      new_nuc_list->push_back(*temp_nuc_iter);
    }
    else {
      if((*new_nuc_list)[overlap_index]._paired==true) {
	if((overlap_type==1 || overlap_type==2) && temp_nuc_iter->_f_strength==0)
	  continue;
	if((overlap_type==3 || overlap_type==4) && temp_nuc_iter->_r_strength==0)
	  continue;
      }
      if((*new_nuc_list)[overlap_index]._paired==false && temp_nuc_iter->_paired==false) {
	if((*new_nuc_list)[overlap_index]._overlap==false) {
	  if((*new_nuc_list)[overlap_index]._f_strength==0 && temp_nuc_iter->_f_strength==0)
	    continue;
	  if((*new_nuc_list)[overlap_index]._r_strength==0 && temp_nuc_iter->_r_strength==0)
	    continue;
	}
	else {
	  if((overlap_type==1 || overlap_type==2) && temp_nuc_iter->_f_strength==0)
	    continue;
	  if((overlap_type==3 || overlap_type==4) && temp_nuc_iter->_r_strength==0)
	    continue;
	}
      }

      new_nuc_center=(*new_nuc_list)[overlap_index]._start+(unsigned int)(((*new_nuc_list)[overlap_index]._end-(*new_nuc_list)[overlap_index]._start)/2.0);
      temp_nuc_center=temp_nuc_iter->_start+(unsigned int)((temp_nuc_iter->_end-temp_nuc_iter->_start)/2.0);

      if(overlap_type==1) {
	center=temp_nuc_center+(unsigned int)((new_nuc_center-temp_nuc_center)*
		 ((*new_nuc_list)[overlap_index]._strength/double((*new_nuc_list)[overlap_index]._strength+temp_nuc_iter->_strength)));
      }
      else if(overlap_type==2) {
	center=new_nuc_center+(unsigned int)((temp_nuc_center-new_nuc_center)*
	  (temp_nuc_iter->_strength/double((*new_nuc_list)[overlap_index]._strength+temp_nuc_iter->_strength)));
      }
      else if(overlap_type==3) {
	center=new_nuc_center+(unsigned int)((temp_nuc_center-new_nuc_center)*
		 (temp_nuc_iter->_strength/double((*new_nuc_list)[overlap_index]._strength+temp_nuc_iter->_strength)));
      }
      else if(overlap_type==4) {
	center=temp_nuc_center+(unsigned int)((new_nuc_center-temp_nuc_center)*
	  ((*new_nuc_list)[overlap_index]._strength/double((*new_nuc_list)[overlap_index]._strength+temp_nuc_iter->_strength)));
      }
      else {
	cerr<<"nuc overlap error"<<endl;
	exit(1);
      }

      (*new_nuc_list)[overlap_index]._start=center-(unsigned int)(DEFAULT_NUC_WIDTH/2.0);
      (*new_nuc_list)[overlap_index]._end=center+(unsigned int)(DEFAULT_NUC_WIDTH/2.0);
      (*new_nuc_list)[overlap_index]._overlap=true;

      if(DEFAULT_NUC_WIDTH%2==1) {
	if((*new_nuc_list)[overlap_index]._occupancy>temp_nuc_iter->_occupancy) {
	  if(overlap_type==1 || overlap_type==4)
	    (*new_nuc_list)[overlap_index]._end=center+(unsigned int)(DEFAULT_NUC_WIDTH/2.0)+1;
	  else
	    (*new_nuc_list)[overlap_index]._start=center-(unsigned int)(DEFAULT_NUC_WIDTH/2.0)-1;
	}
	else {
	  if(overlap_type==2 || overlap_type==3)
	    (*new_nuc_list)[overlap_index]._end=center+(unsigned int)(DEFAULT_NUC_WIDTH/2.0)+1;
	  else
	    (*new_nuc_list)[overlap_index]._start=center-(unsigned int)(DEFAULT_NUC_WIDTH/2.0)-1;
	}
      }
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
