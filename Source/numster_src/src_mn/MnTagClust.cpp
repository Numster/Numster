#include "MnTagClust.h"

#include <math.h>
#include <stdlib.h>

extern double LOCAL_NOISE_FACTOR;
extern double GLOBAL_NOISE_FACTOR;
extern double OVERLAP_RATE;
extern int MN_INIT_CLUSTER_SIZE;
extern double MERGE_DISTANCE;
extern int DEFAULT_NUC_LOWER_BOUND;
extern int DEFAULT_NUC_UPPER_BOUND;

namespace mn {
bool cmp_cluster_center(Cluster cluster1,Cluster cluster2) {
  return cluster1._center<cluster2._center;
}

bool cmp_cluster_strength(Cluster cluster1,Cluster cluster2) {
  if(cluster1._max_read==cluster2._max_read) {
    if(cluster1._strength==cluster2._strength) {
      if(cluster1._point.size()==cluster2._point.size())
	return abs(int(cluster1._center-cluster1._peak))>abs(int(cluster2._center-cluster2._peak));
      else
	return cluster1._point.size()>cluster2._point.size();
    }
    else
      return cluster1._strength<cluster2._strength;
  }
  else
    return cluster1._max_read<cluster2._max_read;
}

bool cmp_tag_start(Tag tag1,Tag tag2) {
  return tag1._start<tag2._start;
}

bool cmp_tag_strength(Tag tag1,Tag tag2) {
  return tag1._strength<tag2._strength;
}

bool cmp_tag_read(Tag tag1,Tag tag2) {
  return tag1._read<tag2._read;
}

TagClust::TagClust(vector<Tag>* tag_list,double size) {
  _f_tag_list=new vector<Tag>();
  _r_tag_list=new vector<Tag>();

  sort(tag_list->begin(),tag_list->end(),cmp_tag_start);

  for(vector<Tag>::iterator tag_iter=tag_list->begin();tag_iter!=tag_list->end();++tag_iter) {
    if(tag_iter->_strand=='F')
      _f_tag_list->push_back(*tag_iter);
    else
      _r_tag_list->push_back(*tag_iter);
  }

  cal_tag_score(size);
}

TagClust::~TagClust() {
  delete _f_tag_list;
  delete _r_tag_list;
  delete _f_cluster_list;
  delete _r_cluster_list;
  _f_tag_list=0;
  _r_tag_list=0;
  _f_cluster_list=0;
  _r_cluster_list=0;
}

void TagClust::cal_tag_score(double size) {
  unsigned int lower_limit=0;
  unsigned int upper_limit=0;
  unsigned int total_read=0;
  vector<Tag>::iterator tmp_tag_iter;
  vector<unsigned int> tmp_tag_start;
  double std=0.0;

  for(int i=0;i<_f_tag_list->size();++i) {
    lower_limit=(*_f_tag_list)[i]._start-MN_INIT_CLUSTER_SIZE;
    upper_limit=(*_f_tag_list)[i]._start+MN_INIT_CLUSTER_SIZE;
    if(lower_limit<0)
      lower_limit=0;
    if(upper_limit>=_f_tag_list->size())
      upper_limit=_f_tag_list->size()-1;

    total_read=0;
    tmp_tag_start.clear();
    for(int j=(i-1);j>=0;--j) {
      if((*_f_tag_list)[j]._start<lower_limit)
	break;
      else {
	total_read+=(*_f_tag_list)[j]._read;
	tmp_tag_start.push_back((*_f_tag_list)[j]._start);
      }
    }
    for(int j=(i+1);j<_f_tag_list->size();++j) {
      if((*_f_tag_list)[j]._start>upper_limit)
	break;
      else {
	total_read+=(*_f_tag_list)[j]._read;
	tmp_tag_start.push_back((*_f_tag_list)[j]._start);
      }
    }
    std=cal_center_std(&tmp_tag_start,(*_f_tag_list)[i]._start);
    (*_f_tag_list)[i]._strength=(*_f_tag_list)[i]._read+(total_read+(1.0-std/(2.0*MN_INIT_CLUSTER_SIZE)))/sqrt(size);
  }

  for(int i=0;i<_r_tag_list->size();++i) {
    lower_limit=(*_r_tag_list)[i]._start-MN_INIT_CLUSTER_SIZE;
    upper_limit=(*_r_tag_list)[i]._start+MN_INIT_CLUSTER_SIZE;
    if(lower_limit<0)
      lower_limit=0;
    if(upper_limit>=_r_tag_list->size())
      upper_limit=_r_tag_list->size()-1;

    total_read=0;
    tmp_tag_start.clear();
    for(int j=(i-1);j>=0;--j) {
      if((*_r_tag_list)[j]._start<lower_limit)
	break;
      else {
	total_read+=(*_r_tag_list)[j]._read;
	tmp_tag_start.push_back((*_r_tag_list)[j]._start);
      }
    }
    for(int j=(i+1);j<_r_tag_list->size();++j) {
      if((*_r_tag_list)[j]._start>upper_limit)
	break;
      else {
	total_read+=(*_r_tag_list)[j]._read;
	tmp_tag_start.push_back((*_r_tag_list)[j]._start);
      }
    }
    std=cal_center_std(&tmp_tag_start,(*_r_tag_list)[i]._start);
    (*_r_tag_list)[i]._strength=(*_r_tag_list)[i]._read+(total_read+(1.0-std/(2.0*MN_INIT_CLUSTER_SIZE)))/sqrt(size);
  }
}

double TagClust::cal_center_std(vector<unsigned int>* value,unsigned int center) {
  double dev=0.0;
  for(int i=0;i<value->size();++i)
    dev+=((*value)[i]-center)*((*value)[i]-center);
  if(value->size()==0)
    return 0.0;
  else
    return sqrt(dev/value->size());
}

void TagClust::determine_center(bool fwd) {
  vector<Tag>* tag_list=0;
  if(fwd==true)
    tag_list=_f_tag_list;
  else
    tag_list=_r_tag_list;

  unsigned int tag_start=0;
  unsigned int center_size=0;
  unsigned int lower_limit=0;
  unsigned int upper_limit=0;
  int cluster_half_width=MN_INIT_CLUSTER_SIZE;

  vector<int> center;
  vector<Tag> tmp_tag_list(*tag_list);

  unsigned int start=0;
  unsigned int end=0;
  bool flag=false;

  sort(tag_list->begin(),tag_list->end(),cmp_tag_strength);
  sort(tmp_tag_list.begin(),tmp_tag_list.end(),cmp_tag_start);

  vector<Cluster> *cluster_list=new vector<Cluster>();

  for(vector<Tag>::reverse_iterator iter=tag_list->rbegin();iter!=tag_list->rend();++iter) {
    tag_start=iter->_start;
    center_size=center.size();

    if(center_size<tag_start) {
      center.resize(tag_start+cluster_half_width);
      center_size=center.size();
    }
    else {
      if(center[tag_start-1]==-1)
	continue;
    }
    center[tag_start-1]=1;

    lower_limit=tag_start-cluster_half_width;
    if(tag_start<(unsigned int)(1+cluster_half_width))
      lower_limit=1;
    start=lower_limit;
    for(unsigned int i=tag_start-1;i>(lower_limit-1);--i) {
      if(center[i-1]==-1) {
	start=i+1;
	break;
      }
      center[i-1]=-1;
    }

    upper_limit=tag_start+cluster_half_width;
    end=upper_limit;
    for(unsigned int i=tag_start;i<upper_limit;++i) {
      if(center[i]==-1) {
	end=i;
	break;
      }
      center[i]=-1;
    }

    Cluster cluster;
    cluster._id=-1;
    cluster._peak=tag_start;
    cluster._start=start;
    cluster._end=end;
    for(vector<Tag>::iterator tmp_iter=tmp_tag_list.begin();tmp_iter!=tmp_tag_list.end();) {
      if(tmp_iter->_start>end)
	break;

      if(tmp_iter->_start>=start) {
	cluster._point.push_back(*tmp_iter);
	tmp_iter=tmp_tag_list.erase(tmp_iter);
      }
      else
	++tmp_iter;
    }
    if(flag==false && cluster._point.size()==0) {
      flag=true;
      cerr<<"!! warning: there exist multiple tags with the same start position"<<endl;
    }
    if(cluster._point.size()>0)
      cluster_list->push_back(cluster);
  }

  cal_statistic(cluster_list);

  if(fwd==true)
    _f_cluster_list=cluster_list;
  else
    _r_cluster_list=cluster_list;
}

void TagClust::create_local_cluster(bool fwd) {
  vector<Cluster>* cluster_list=0;
  if(fwd==true)
    cluster_list=_f_cluster_list;
  else
    cluster_list=_r_cluster_list;

  double left_distance=0;
  double right_distance=0;

  sort(cluster_list->begin(),cluster_list->end(),cmp_cluster_center);
  vector<Cluster> tmp_cluster_list(*cluster_list);

  /*
  for(int i=0;i<cluster_list->size();++i) {
    cout<<i<<","<<(*cluster_list)[i]._start<<","<<(*cluster_list)[i]._end<<","<<(*cluster_list)[i]._max_read<<","<<(*cluster_list)[i]._left_read<<","<<(*cluster_list)[i]._right_read<<endl;
  }
  */

  for(int i=0,j=0;i<tmp_cluster_list.size();++i,++j) {
    left_distance=tmp_cluster_list[i]._peak-tmp_cluster_list[i]._start;
    right_distance=tmp_cluster_list[i]._end-tmp_cluster_list[i]._peak;

    //left-merge from right-side
    if(j>0 && right_distance!=0 && (((*cluster_list)[j-1]._right_read>(*cluster_list)[j]._max_read && left_distance/right_distance<MERGE_DISTANCE) || 
				    ((*cluster_list)[j-1]._right_read==(*cluster_list)[j]._max_read && left_distance/right_distance<(MERGE_DISTANCE/2.0)))) { //merge
      (*cluster_list)[j-1]._point.insert((*cluster_list)[j-1]._point.end(),(*cluster_list)[j]._point.begin(),(*cluster_list)[j]._point.end());
      (*cluster_list)[j-1]._right_read=(*cluster_list)[j]._right_read;
      (*cluster_list)[j-1]._end=(*cluster_list)[j]._end;
      cluster_list->erase(cluster_list->begin()+j);
      --j;
    }//right-merge from left-side
    else if(i<tmp_cluster_list.size()-1 && left_distance!=0 && (((*cluster_list)[j+1]._left_read>(*cluster_list)[j]._max_read && right_distance/left_distance<MERGE_DISTANCE) 
|| ((*cluster_list)[j+1]._left_read==(*cluster_list)[j]._max_read && right_distance/left_distance<(MERGE_DISTANCE/2.0)))) { //merge
      (*cluster_list)[j+1]._point.insert((*cluster_list)[j+1]._point.begin(),(*cluster_list)[j]._point.begin(),(*cluster_list)[j]._point.end());
      (*cluster_list)[j+1]._left_read=(*cluster_list)[j]._left_read;
      (*cluster_list)[j+1]._start=(*cluster_list)[j]._start;
      cluster_list->erase(cluster_list->begin()+j);
      --j;
    }
  }

  /*
  for(int i=0;i<cluster_list->size();++i) {
    cout<<i<<","<<(*cluster_list)[i]._start<<","<<(*cluster_list)[i]._end<<","<<(*cluster_list)[i]._max_read<<","<<(*cluster_list)[i]._left_read<<","<<(*cluster_list)[i]._right_read<<endl;
  }
  */

  cal_statistic(cluster_list);
  clean_noisy_cluster(cluster_list);
}

void TagClust::cal_corresponding_cluster(vector<Cluster>* f_cluster_list,vector<Cluster>* r_cluster_list,unsigned int upper_limit) {
  bool flag=false;
  unsigned int lower_limit=(f_cluster_list->front()._center+NUC_MIN_WIDTH);
  if(upper_limit==0)
    upper_limit=_r_cluster_list->back()._center;

  for(vector<Cluster>::iterator f_cluster_iter=f_cluster_list->begin();f_cluster_iter!=f_cluster_list->end();++f_cluster_iter)
    f_cluster_iter->_id=-1;

  for(vector<Cluster>::iterator r_cluster_iter=_r_cluster_list->begin();r_cluster_iter!=_r_cluster_list->end();++r_cluster_iter) {
    if(r_cluster_iter->_center>=lower_limit)
      flag=true;
    if(flag==true && r_cluster_iter->_center<=upper_limit) {
      r_cluster_iter->_id=-1;
      r_cluster_list->push_back(*r_cluster_iter);
    }
    if(flag==true && r_cluster_iter->_center>upper_limit)
      break;
  }

  sort(f_cluster_list->begin(),f_cluster_list->end(),cmp_cluster_strength);
  sort(r_cluster_list->begin(),r_cluster_list->end(),cmp_cluster_strength);

  int id=0;
  for(int i=(f_cluster_list->size()-1);i>=0;--i) {
    lower_limit=((*f_cluster_list)[i]._center+NUC_MIN_WIDTH);
    upper_limit=((*f_cluster_list)[i]._center+NUC_MAX_WIDTH);

    for(int j=(r_cluster_list->size()-1);j>=0;--j) {
      if((*r_cluster_list)[j]._center>=lower_limit && (*r_cluster_list)[j]._center<=upper_limit) {
	if((*r_cluster_list)[j]._id!=-1)
	  continue;
	(*f_cluster_list)[i]._id=id;
	(*r_cluster_list)[j]._id=id;
	++id;
	break;
      }
    }
  }

  maximize_cluster_slimiarity(id,f_cluster_list,r_cluster_list);
  maximize_cluster_slimiarity(id,r_cluster_list,f_cluster_list);
}

void TagClust::maximize_cluster_slimiarity(int size,vector<Cluster>* target_cluster_list,vector<Cluster>* corresponding_cluster_list) {
  sort(target_cluster_list->begin(),target_cluster_list->end(),cmp_cluster_center);
  sort(corresponding_cluster_list->begin(),corresponding_cluster_list->end(),cmp_cluster_center);

  vector<Tag> tmp_tag_list;
  unsigned int corresponding_center=0;
  unsigned short corresponding_strength=0;
  int left_distance=0;
  int right_distance=0;
  bool left_flag=false;
  bool right_flag=false;
  int original_diff=0;
  int tmp_diff=0;
  unsigned int tmp_center=0;

  for(int i=0;i<size;++i) {
    corresponding_strength=0;
    for(int j=0;j<corresponding_cluster_list->size();++j) {
      if((*corresponding_cluster_list)[j]._id==i) {
	corresponding_center=(*corresponding_cluster_list)[j]._center;
	corresponding_strength=(*corresponding_cluster_list)[j]._strength;
	break;
      }
    }

    if(corresponding_strength==0) {
      cerr<<"maximization error."<<endl;
      exit(1);
    }

    for(int j=0;j<target_cluster_list->size();++j) {
      if((*target_cluster_list)[j]._id==i) {
	while(true) {
	  left_flag=false;
	  right_flag=false;
	  left_distance=0;
	  right_distance=0;

	  if(j>0 && (*target_cluster_list)[j-1]._id==-1) {
	    left_distance=(*target_cluster_list)[j]._point.front()._start-(*target_cluster_list)[j-1]._point.back()._start;
	    left_flag=true;
	  }
	  if((j+1)<target_cluster_list->size() && (*target_cluster_list)[j+1]._id==-1) {
	    right_distance=(*target_cluster_list)[j+1]._point.front()._start-(*target_cluster_list)[j]._point.back()._start;
	    right_flag=true;
	  }

	  if(left_flag==false && right_flag==false)
	    break;

	  if(left_flag==true && right_flag==true && left_distance<=right_distance) {
	    right_flag=false;
	  }
	  else if(left_flag==true && right_flag==true && left_distance>right_distance) {
	    left_flag=false;
	  }

	  if(left_flag==true && right_flag==true) {
	    cerr<<"maximization error."<<endl;
	    exit(1);
	  }

	  if(left_flag==true) {
	    original_diff=abs((int)((*target_cluster_list)[j]._strength-corresponding_strength));
	    tmp_diff=abs((int)((*target_cluster_list)[j-1]._strength+(*target_cluster_list)[j]._strength-corresponding_strength));
	    if(tmp_diff<original_diff) {
	      tmp_tag_list.clear();
	      tmp_tag_list.insert(tmp_tag_list.end(),(*target_cluster_list)[j-1]._point.begin(),(*target_cluster_list)[j-1]._point.end());
	      tmp_tag_list.insert(tmp_tag_list.end(),(*target_cluster_list)[j]._point.begin(),(*target_cluster_list)[j]._point.end());
	      tmp_center=cal_center(&tmp_tag_list,(*target_cluster_list)[j-1]._strength+(*target_cluster_list)[j]._strength);

	      if(abs((int)(tmp_center-corresponding_center))<=NUC_MAX_WIDTH && abs((int)(tmp_center-corresponding_center))>=NUC_MIN_WIDTH) {
		(*target_cluster_list)[j]._point.insert((*target_cluster_list)[j]._point.begin(),
							(*target_cluster_list)[j-1]._point.begin(),(*target_cluster_list)[j-1]._point.end());
		(*target_cluster_list)[j]._center=tmp_center;
		(*target_cluster_list)[j]._strength=(*target_cluster_list)[j-1]._strength+(*target_cluster_list)[j]._strength;
		(*target_cluster_list)[j]._start=(*target_cluster_list)[j-1]._start;
		(*target_cluster_list)[j]._left_read=(*target_cluster_list)[j-1]._max_read;
		if((*target_cluster_list)[j-1]._max_read>(*target_cluster_list)[j]._max_read) {
		  (*target_cluster_list)[j]._max_read=(*target_cluster_list)[j-1]._max_read;
		  (*target_cluster_list)[j]._peak=(*target_cluster_list)[j-1]._peak;
		}
		target_cluster_list->erase(target_cluster_list->begin()+(j-1));
		--j;
	      }
	      else //outof range
		break;
	    }
	    else //no better
	      break;
	  }
	  else {
	    original_diff=abs((int)((*target_cluster_list)[j]._strength-corresponding_strength));
	    tmp_diff=abs((int)((*target_cluster_list)[j+1]._strength+(*target_cluster_list)[j]._strength-corresponding_strength));
	    if(tmp_diff<original_diff) {
	      tmp_tag_list.clear();
	      tmp_tag_list.insert(tmp_tag_list.end(),(*target_cluster_list)[j]._point.begin(),(*target_cluster_list)[j]._point.end());
	      tmp_tag_list.insert(tmp_tag_list.end(),(*target_cluster_list)[j+1]._point.begin(),(*target_cluster_list)[j+1]._point.end());
	      tmp_center=cal_center(&tmp_tag_list,(*target_cluster_list)[j]._strength+(*target_cluster_list)[j+1]._strength);

	      if(abs((int)(tmp_center-corresponding_center))<=NUC_MAX_WIDTH && abs((int)(tmp_center-corresponding_center))>=NUC_MIN_WIDTH) {
		(*target_cluster_list)[j]._point.insert((*target_cluster_list)[j]._point.end(),
							(*target_cluster_list)[j+1]._point.begin(),(*target_cluster_list)[j+1]._point.end());
		(*target_cluster_list)[j]._center=tmp_center;
		(*target_cluster_list)[j]._strength=(*target_cluster_list)[j]._strength+(*target_cluster_list)[j+1]._strength;
		(*target_cluster_list)[j]._end=(*target_cluster_list)[j+1]._start;
		(*target_cluster_list)[j]._right_read=(*target_cluster_list)[j+1]._max_read;
		if((*target_cluster_list)[j+1]._max_read>(*target_cluster_list)[j]._max_read) {
		  (*target_cluster_list)[j]._max_read=(*target_cluster_list)[j+1]._max_read;
		  (*target_cluster_list)[j]._peak=(*target_cluster_list)[j+1]._peak;
		}
		target_cluster_list->erase(target_cluster_list->begin()+(j+1));
	      }
	      else //outof range
		break;
	    }
	    else //no better
	      break;
	  }
	}
      }
    }
  }
}

void TagClust::create_global_cluster() {
  vector<Cluster> *f_cluster_list=new vector<Cluster>();
  vector<Cluster> *r_cluster_list=new vector<Cluster>();

  sort(_f_cluster_list->begin(),_f_cluster_list->end(),cmp_cluster_center);
  sort(_r_cluster_list->begin(),_r_cluster_list->end(),cmp_cluster_center);

  vector<Cluster> tmp_f_cluster_list;
  vector<Cluster> tmp_r_cluster_list;
  vector<Cluster>::iterator f_cluster_iter=_f_cluster_list->begin();
  tmp_f_cluster_list.push_back(*f_cluster_iter);

  int count=0;

  for(++f_cluster_iter;f_cluster_iter!=_f_cluster_list->end();++f_cluster_iter) {
    if(f_cluster_iter==(_f_cluster_list->end()-1)) {
      unsigned int upper_limit=0;
      if((f_cluster_iter->_center-(f_cluster_iter-1)->_center)<=(NUC_MAX_WIDTH-NUC_MIN_WIDTH)) {
	tmp_f_cluster_list.push_back(*f_cluster_iter);
	upper_limit=0;
      }
      else
	upper_limit=f_cluster_iter->_center+NUC_MIN_WIDTH-1;

      cal_corresponding_cluster(&tmp_f_cluster_list,&tmp_r_cluster_list,upper_limit);
      f_cluster_list->insert(f_cluster_list->end(),tmp_f_cluster_list.begin(),tmp_f_cluster_list.end());
      r_cluster_list->insert(r_cluster_list->end(),tmp_r_cluster_list.begin(),tmp_r_cluster_list.end());
      tmp_f_cluster_list.clear();
      tmp_r_cluster_list.clear();

      if((f_cluster_iter->_center-(f_cluster_iter-1)->_center)>(NUC_MAX_WIDTH-NUC_MIN_WIDTH)) {
	tmp_f_cluster_list.push_back(*f_cluster_iter);
	cal_corresponding_cluster(&tmp_f_cluster_list,&tmp_r_cluster_list,0);
	f_cluster_list->insert(f_cluster_list->end(),tmp_f_cluster_list.begin(),tmp_f_cluster_list.end());
	r_cluster_list->insert(r_cluster_list->end(),tmp_r_cluster_list.begin(),tmp_r_cluster_list.end());
      }
    }
    else if((f_cluster_iter->_center-(f_cluster_iter-1)->_center)<=(NUC_MAX_WIDTH-NUC_MIN_WIDTH))
      tmp_f_cluster_list.push_back(*f_cluster_iter);
    else {
      cal_corresponding_cluster(&tmp_f_cluster_list,&tmp_r_cluster_list,f_cluster_iter->_center+NUC_MIN_WIDTH-1);

      f_cluster_list->insert(f_cluster_list->end(),tmp_f_cluster_list.begin(),tmp_f_cluster_list.end());
      r_cluster_list->insert(r_cluster_list->end(),tmp_r_cluster_list.begin(),tmp_r_cluster_list.end());
      tmp_f_cluster_list.clear();
      tmp_r_cluster_list.clear();
      tmp_f_cluster_list.push_back(*f_cluster_iter);
    }
  }

  delete _f_cluster_list;
  delete _r_cluster_list;
  _f_cluster_list=f_cluster_list;
  _r_cluster_list=r_cluster_list;

  sort(_f_tag_list->begin(),_f_tag_list->end(),cmp_tag_start);
  sort(_r_tag_list->begin(),_r_tag_list->end(),cmp_tag_start);
  clean_no_matching_cluster(_f_cluster_list,true);
  clean_no_matching_cluster(_r_cluster_list,false);
}

unsigned int TagClust::cal_center(vector<Tag>* tag_list,unsigned short strength) { //median
  sort(tag_list->begin(),tag_list->end(),cmp_tag_start);

  double median_read=strength/2.0;
  unsigned short read=0;
  unsigned int total_read=0;
  unsigned int center=0;
  for(vector<Tag>::iterator tag_iter=tag_list->begin();tag_iter!=tag_list->end();++tag_iter) {
    read=tag_iter->_read;
    total_read+=read;
    if(total_read>=median_read) {
      center=tag_iter->_start;
      break;
    }
  }
  return center;
}

void TagClust::cal_statistic(vector<Cluster>* cluster_list) {
  unsigned short read=0;
  unsigned short max_read=0;
  unsigned int total_read=0;

  vector<Tag>* point_list;
  for(vector<Cluster>::iterator cluster_iter=cluster_list->begin();cluster_iter!=cluster_list->end();++cluster_iter) {
    read=0;
    max_read=0;
    total_read=0;
    point_list=&(cluster_iter->_point);

    sort(point_list->begin(),point_list->end(),cmp_tag_start);

    for(vector<Tag>::iterator tag_iter=point_list->begin();tag_iter!=point_list->end();++tag_iter) {
      read=tag_iter->_read;
      if(read>max_read)
	max_read=read;
      total_read+=read;
    }

    cluster_iter->_max_read=max_read;
    cluster_iter->_left_read=max_read;
    cluster_iter->_right_read=max_read;
    cluster_iter->_strength=total_read;
    cluster_iter->_center=cal_center(point_list,total_read);
  }
}

int TagClust::cal_quantile(vector<int>* strength_list,double q) {
  sort(strength_list->begin(),strength_list->end());

  int total=0;
  for(int i=0;i<strength_list->size();++i)
    total+=(*strength_list)[i];
  double q_value=total*q;
  int value=0;

  total=0;
  for(int i=0;i<strength_list->size();++i) {
    total+=(*strength_list)[i];
    if(total>=q_value) {
      value=(*strength_list)[i];
      break;
    }
  }
  return value;
}

void TagClust::clean_no_matching_cluster(vector<Cluster>* cluster_list,bool fwd) {
  vector<int> max_read;
  for(vector<Cluster>::iterator cluster_iter=cluster_list->begin();cluster_iter!=cluster_list->end();++cluster_iter)
      max_read.push_back(cluster_iter->_max_read);

  double q1=cal_quantile(&max_read,0.25);
  double q3=cal_quantile(&max_read,0.75);
  double cutoff=q1-GLOBAL_NOISE_FACTOR*(q3-q1);
  if(q1==q3)
    cutoff=0.0;

  int removed_cluster_count=0;
  int original_cluster_count= cluster_list->size();

  if(UNPAIRED_NUC==true) {
    for(vector<Cluster>::iterator cluster_iter=cluster_list->begin();cluster_iter!=cluster_list->end();) {
      if(cluster_iter->_id==-1) {
	vector<Tag>* tag_list=0;
	unsigned int lower_limit=0;
	unsigned int upper_limit=0;
	if(fwd==true) {
	  tag_list=_r_tag_list;
	  lower_limit=cluster_iter->_center+(unsigned int)(DEFAULT_NUC_WIDTH-(DEFAULT_NUC_WIDTH-DEFAULT_NUC_LOWER_BOUND));
	  upper_limit=cluster_iter->_center+(unsigned int)(DEFAULT_NUC_WIDTH+(DEFAULT_NUC_UPPER_BOUND-DEFAULT_NUC_WIDTH));
	}
	else {
	  tag_list=_f_tag_list;
	  lower_limit=cluster_iter->_center-(unsigned int)(DEFAULT_NUC_WIDTH+(DEFAULT_NUC_UPPER_BOUND-DEFAULT_NUC_WIDTH));
	  upper_limit=cluster_iter->_center-(unsigned int)(DEFAULT_NUC_WIDTH-(DEFAULT_NUC_WIDTH-DEFAULT_NUC_LOWER_BOUND));
	}
	bool flag=false;
	for(vector<Tag>::iterator tag_iter=tag_list->begin();tag_iter!=tag_list->end();++tag_iter) {
	  if(tag_iter->_start>=lower_limit && tag_iter->_start<=upper_limit) {
	    flag=true;
	    break;
	  }
	  if(tag_iter->_start>upper_limit)
	    break;
	}

	if(flag==false || cluster_iter->_max_read<=cutoff) {
	  cluster_iter=cluster_list->erase(cluster_iter);
	  ++removed_cluster_count;
	}
	else
	  ++cluster_iter;
      }
      else
	++cluster_iter;
    }
  }
  else {
    for(vector<Cluster>::iterator cluster_iter=cluster_list->begin();cluster_iter!=cluster_list->end();) {
      if(cluster_iter->_id==-1) {
	cluster_iter=cluster_list->erase(cluster_iter);
	++removed_cluster_count;
      }
      else
	++cluster_iter;
    }
  }
  //cout<<"   - no matching cluster count: "<<removed_cluster_count<<"("<<removed_cluster_count/double(original_cluster_count)*100<<"%)"<<endl;
}

void TagClust::clean_noisy_cluster(vector<Cluster>* cluster_list) {
  vector<int> strength;
  for(vector<Cluster>::iterator cluster_iter=cluster_list->begin();cluster_iter!=cluster_list->end();++cluster_iter)
    strength.push_back(cluster_iter->_strength);

  double q1=cal_quantile(&strength,0.25);
  double q3=cal_quantile(&strength,0.75);

  double max_cutoff=q3-q1;
  double cutoff=0;
  if(OVERLAP_RATE<=0.3)
    cutoff=q1-LOCAL_NOISE_FACTOR*(q3-q1);
  else if(OVERLAP_RATE<=0.35)
    cutoff=q1-(LOCAL_NOISE_FACTOR-0.05)*(q3-q1);
  else if(OVERLAP_RATE<=0.4)
    cutoff=q1-(LOCAL_NOISE_FACTOR-0.1)*(q3-q1);
  else if(OVERLAP_RATE<=0.45)
    cutoff=q1-(LOCAL_NOISE_FACTOR-0.2)*(q3-q1);
  else
    cutoff=q1-(LOCAL_NOISE_FACTOR-0.25)*(q3-q1);

  if(cutoff>max_cutoff)
    cutoff=max_cutoff;

  int removed_cluster_count=0;
  int original_cluster_count= cluster_list->size();
  for(vector<Cluster>::iterator cluster_iter=cluster_list->begin();cluster_iter!=cluster_list->end();) {
    if(cluster_iter->_strength<=cutoff) {
      cluster_iter=cluster_list->erase(cluster_iter);
      ++removed_cluster_count;
    }
    else
      ++cluster_iter;
  }
  //cout<<"   - noisy cluster count: "<<removed_cluster_count<<"("<<removed_cluster_count/double(original_cluster_count)*100<<"%)"<<endl;
}

vector<Cluster>* TagClust::get_f_cluster_list() {
  return _f_cluster_list;
}

vector<Cluster>* TagClust::get_r_cluster_list() {
  return _r_cluster_list;
}
}
