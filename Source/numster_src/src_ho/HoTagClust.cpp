#include "HoTagClust.h"

#include <math.h>
#include <stdlib.h>

extern double LOCAL_NOISE_FACTOR;
extern double GLOBAL_NOISE_FACTOR;
extern double OVERLAP_RATE;
extern int HO_INIT_CLUSTER_SIZE;
extern double MERGE_DISTANCE;
extern double EUCLIDEAN_DISTANCE;
extern int SEARCH_LOWER_BOUND;
extern int SEARCH_UPPER_BOUND;

namespace ho {
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
    lower_limit=(*_f_tag_list)[i]._start-HO_INIT_CLUSTER_SIZE;
    upper_limit=(*_f_tag_list)[i]._start+HO_INIT_CLUSTER_SIZE;
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
    (*_f_tag_list)[i]._strength=(*_f_tag_list)[i]._read+(total_read+(1.0-std/(2.0*HO_INIT_CLUSTER_SIZE)))/sqrt(size);
  }

  for(int i=0;i<_r_tag_list->size();++i) {
    lower_limit=(*_r_tag_list)[i]._start-HO_INIT_CLUSTER_SIZE;
    upper_limit=(*_r_tag_list)[i]._start+HO_INIT_CLUSTER_SIZE;
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
    (*_r_tag_list)[i]._strength=(*_r_tag_list)[i]._read+(total_read+(1.0-std/(2.0*HO_INIT_CLUSTER_SIZE)))/sqrt(size);
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
  int cluster_half_width=HO_INIT_CLUSTER_SIZE;

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
    else if(i<tmp_cluster_list.size()-1 && left_distance!=0 && (((*cluster_list)[j+1]._left_read>(*cluster_list)[j]._max_read && right_distance/left_distance<MERGE_DISTANCE) || ((*cluster_list)[j+1]._left_read==(*cluster_list)[j]._max_read && right_distance/left_distance<(MERGE_DISTANCE/2.0)))) { //merge
      (*cluster_list)[j+1]._point.insert((*cluster_list)[j+1]._point.begin(),(*cluster_list)[j]._point.begin(),(*cluster_list)[j]._point.end());
      (*cluster_list)[j+1]._left_read=(*cluster_list)[j]._left_read;
      (*cluster_list)[j+1]._start=(*cluster_list)[j]._start;
      cluster_list->erase(cluster_list->begin()+j);
      --j;
    }
  }

  cal_statistic(cluster_list);
  clean_noisy_cluster(cluster_list);
}

int TagClust::classify_id(vector<Cluster> &f_cluster_list,vector<Cluster> &r_cluster_list,unsigned int f_peak,unsigned int r_peak,unsigned short f_peak_read,unsigned short r_peak_read) {
  double type1[]={0.3858,0.1739,0.2554,0.1847};
  double type2[]={0.1847,0.2554,0.1739,0.3858};
  double type3[]={0.3858,0.1739,0.1847,0.2554};
  double type4[]={0.2554,0.1847,0.1739,0.3858};
  double type5[]={0.1739,0.3858,0.2554,0.1847};
  double type6[]={0.1847,0.2554,0.3858,0.1739};

  double f_type1[4]={0};
  double r_type1[4]={0};
  double f_read1[4]={0};
  double r_read1[4]={0};

  double f_type2[4]={0};
  double r_type2[4]={0};
  double f_read2[4]={0};
  double r_read2[4]={0};

  f_read1[0]=f_peak_read;
  f_read2[1]=f_peak_read;
  r_read1[2]=r_peak_read;
  r_read2[3]=r_peak_read;

  for(int i=0;i<f_cluster_list.size();++i) {
    if(f_peak>f_cluster_list[i]._end+7)
      continue;
    if(r_peak+12<f_cluster_list[i]._start)
      break;

    if(f_peak>=f_cluster_list[i]._start+7 && f_peak<=f_cluster_list[i]._end+7) {
      for(vector<Tag>::iterator tag_iter=f_cluster_list[i]._point.begin();tag_iter!=f_cluster_list[i]._point.end();++tag_iter) {
	if(f_peak==tag_iter->_start+7) {
	  f_read2[0]=tag_iter->_read;
	  break;
	}
      }
    }

    if(f_peak+7>=f_cluster_list[i]._start && f_peak+7<=f_cluster_list[i]._end) {
      for(vector<Tag>::iterator tag_iter=f_cluster_list[i]._point.begin();tag_iter!=f_cluster_list[i]._point.end();++tag_iter) {
	if(f_peak+7==tag_iter->_start) {
	  f_read1[1]=tag_iter->_read;
	  break;
	}
      }
    }

    if(r_peak>=f_cluster_list[i]._start+2 && r_peak<=f_cluster_list[i]._end+2) {
      for(vector<Tag>::iterator tag_iter=f_cluster_list[i]._point.begin();tag_iter!=f_cluster_list[i]._point.end();++tag_iter) {
	if(r_peak==tag_iter->_start+2) {
	  r_read2[0]=tag_iter->_read;
	  break;
	}
      }
    }

    if(r_peak+5>=f_cluster_list[i]._start && r_peak+5<=f_cluster_list[i]._end) {
      for(vector<Tag>::iterator tag_iter=f_cluster_list[i]._point.begin();tag_iter!=f_cluster_list[i]._point.end();++tag_iter) {
	if(r_peak+5==tag_iter->_start) {
	  r_read1[0]=tag_iter->_read;
	  r_read2[1]=tag_iter->_read;
	  break;
	}
      }
    }

    if(r_peak+12>=f_cluster_list[i]._start && r_peak+12<=f_cluster_list[i]._end) {
      for(vector<Tag>::iterator tag_iter=f_cluster_list[i]._point.begin();tag_iter!=f_cluster_list[i]._point.end();++tag_iter) {
	if(r_peak+12==tag_iter->_start) {
	  r_read1[1]=tag_iter->_read;
	  break;
	}
      }
    }
  }

  for(int i=0;i<r_cluster_list.size();++i) {
    if(f_peak>r_cluster_list[i]._end+12)
      continue;
    if(r_peak+7<r_cluster_list[i]._start)
      break;

    if(f_peak>=r_cluster_list[i]._start+12 && f_peak<=r_cluster_list[i]._end+12) {
      for(vector<Tag>::iterator tag_iter=r_cluster_list[i]._point.begin();tag_iter!=r_cluster_list[i]._point.end();++tag_iter) {
	if(f_peak==tag_iter->_start+12) {
	  f_read2[2]=tag_iter->_read;
	  break;
	}
      }
    }

    if(f_peak>=r_cluster_list[i]._start+5 && f_peak<=r_cluster_list[i]._end+5) {
      for(vector<Tag>::iterator tag_iter=r_cluster_list[i]._point.begin();tag_iter!=r_cluster_list[i]._point.end();++tag_iter) {
	if(f_peak==tag_iter->_start+5) {
	  f_read1[2]=tag_iter->_read;
	  f_read2[3]=tag_iter->_read;
	  break;
	}
      }
    }

    if(f_peak+2>=r_cluster_list[i]._start && f_peak+2<=r_cluster_list[i]._end) {
      for(vector<Tag>::iterator tag_iter=r_cluster_list[i]._point.begin();tag_iter!=r_cluster_list[i]._point.end();++tag_iter) {
	if(f_peak+2==tag_iter->_start) {
	  f_read1[3]=tag_iter->_read;
	  break;
	}
      }
    }

    if(r_peak>=r_cluster_list[i]._start+7 && r_peak<=r_cluster_list[i]._end+7) {
      for(vector<Tag>::iterator tag_iter=r_cluster_list[i]._point.begin();tag_iter!=r_cluster_list[i]._point.end();++tag_iter) {
	if(r_peak==tag_iter->_start+7) {
	  r_read2[2]=tag_iter->_read;
	  break;
	}
      }
    }

    if(r_peak+7>=r_cluster_list[i]._start && r_peak+7<=r_cluster_list[i]._end) {
      for(vector<Tag>::iterator tag_iter=r_cluster_list[i]._point.begin();tag_iter!=r_cluster_list[i]._point.end();++tag_iter) {
	if(r_peak+7==tag_iter->_start) {
	  r_read1[3]=tag_iter->_read;
	  break;
	}
      }
    }
  }

  double f_total1=f_read1[0]+f_read1[1]+f_read1[2]+f_read1[3];
  double f_total2=f_read2[0]+f_read2[1]+f_read2[2]+f_read2[3];
  double r_total1=r_read1[0]+r_read1[1]+r_read1[2]+r_read1[3];
  double r_total2=r_read2[0]+r_read2[1]+r_read2[2]+r_read2[3];
  for(int i=0;i<4;++i) {
    f_type1[i]=f_read1[i]/f_total1;
    f_type2[i]=f_read2[i]/f_total2;
    r_type1[i]=r_read1[i]/r_total1;
    r_type2[i]=r_read2[i]/r_total2;
  }

  double f_dist1=sqrt((f_type1[0]-type1[0])*(f_type1[0]-type1[0])+(f_type1[1]-type1[1])*(f_type1[1]-type1[1])+
		    (f_type1[2]-type1[2])*(f_type1[2]-type1[2])+(f_type1[3]-type1[3])*(f_type1[3]-type1[3]));
  double f_dist2=sqrt((f_type2[0]-type2[0])*(f_type2[0]-type2[0])+(f_type2[1]-type2[1])*(f_type2[1]-type2[1])+
		    (f_type2[2]-type2[2])*(f_type2[2]-type2[2])+(f_type2[3]-type2[3])*(f_type2[3]-type2[3]));
  double f_dist3=sqrt((f_type1[0]-type3[0])*(f_type1[0]-type3[0])+(f_type1[1]-type3[1])*(f_type1[1]-type3[1])+
		    (f_type1[2]-type3[2])*(f_type1[2]-type3[2])+(f_type1[3]-type3[3])*(f_type1[3]-type3[3]));
  double f_dist4=sqrt((f_type1[0]-type4[0])*(f_type1[0]-type4[0])+(f_type1[1]-type4[1])*(f_type1[1]-type4[1])+
		    (f_type1[2]-type4[2])*(f_type1[2]-type4[2])+(f_type1[3]-type4[3])*(f_type1[3]-type4[3]));
  double f_dist5=sqrt((f_type2[0]-type5[0])*(f_type2[0]-type5[0])+(f_type2[1]-type5[1])*(f_type2[1]-type5[1])+
		    (f_type2[2]-type5[2])*(f_type2[2]-type5[2])+(f_type2[3]-type5[3])*(f_type2[3]-type5[3]));
  double f_dist6=sqrt((f_type2[0]-type6[0])*(f_type2[0]-type6[0])+(f_type2[1]-type6[1])*(f_type2[1]-type6[1])+
		    (f_type2[2]-type6[2])*(f_type2[2]-type6[2])+(f_type2[3]-type6[3])*(f_type2[3]-type6[3]));

  double r_dist1=sqrt((r_type1[0]-type1[0])*(r_type1[0]-type1[0])+(r_type1[1]-type1[1])*(r_type1[1]-type1[1])+
		    (r_type1[2]-type1[2])*(r_type1[2]-type1[2])+(r_type1[3]-type1[3])*(r_type1[3]-type1[3]));
  double r_dist2=sqrt((r_type2[0]-type2[0])*(r_type2[0]-type2[0])+(r_type2[1]-type2[1])*(r_type2[1]-type2[1])+
		    (r_type2[2]-type2[2])*(r_type2[2]-type2[2])+(r_type2[3]-type2[3])*(r_type2[3]-type2[3]));
  double r_dist3=sqrt((r_type2[0]-type3[0])*(r_type2[0]-type3[0])+(r_type2[1]-type3[1])*(r_type2[1]-type3[1])+
		    (r_type2[2]-type3[2])*(r_type2[2]-type3[2])+(r_type2[3]-type3[3])*(r_type2[3]-type3[3]));
  double r_dist4=sqrt((r_type2[0]-type4[0])*(r_type2[0]-type4[0])+(r_type2[1]-type4[1])*(r_type2[1]-type4[1])+
		    (r_type2[2]-type4[2])*(r_type2[2]-type4[2])+(r_type2[3]-type4[3])*(r_type2[3]-type4[3]));
  double r_dist5=sqrt((r_type1[0]-type5[0])*(r_type1[0]-type5[0])+(r_type1[1]-type5[1])*(r_type1[1]-type5[1])+
		    (r_type1[2]-type5[2])*(r_type1[2]-type5[2])+(r_type1[3]-type5[3])*(r_type1[3]-type5[3]));
  double r_dist6=sqrt((r_type1[0]-type6[0])*(r_type1[0]-type6[0])+(r_type1[1]-type6[1])*(r_type1[1]-type6[1])+
		    (r_type1[2]-type6[2])*(r_type1[2]-type6[2])+(r_type1[3]-type6[3])*(r_type1[3]-type6[3]));

  double f_min=f_dist1;
  double f_min1=f_dist1;
  double f_min2=f_dist2;
  int f_primary_secondary=10;
  if(f_dist2<f_min) {
    f_min=f_dist2;
    f_primary_secondary=20;
  }
  if(f_dist3<=f_min) {
    f_min=f_dist3;
    f_primary_secondary=10;
  }
  if(f_dist4<f_min) {
    f_min=f_dist4;
    f_primary_secondary=10;
  }
  if(f_dist5<f_min) {
    f_min=f_dist5;
    f_primary_secondary=20;
  }
  if(f_dist6<f_min) {
    f_min=f_dist6;
    f_primary_secondary=20;
  }

  if(f_dist3<f_min1)
    f_min1=f_dist3;
  if(f_dist4<f_min1)
    f_min1=f_dist4;
  if(f_dist5<f_min2)
    f_min2=f_dist5;
  if(f_dist6<f_min2)
    f_min2=f_dist6;

  double r_min=r_dist1;
  double r_min1=r_dist2;
  double r_min2=r_dist1;
  int r_primary_secondary=2;
  if(r_dist2<=r_min) {
    r_min=r_dist2;
    r_primary_secondary=1;
  }
  if(r_dist3<r_min) {
    r_min=r_dist3;
    r_primary_secondary=1;
  }
  if(r_dist4<=r_min) {
    r_min=r_dist4;
    r_primary_secondary=1;
  }
  if(r_dist5<r_min) {
    r_min=r_dist5;
    r_primary_secondary=2;
  }
  if(r_dist6<r_min) {
    r_min=r_dist6;
    r_primary_secondary=2;
  }

  if(r_dist3<r_min1)
    r_min1=r_dist3;
  if(r_dist4<r_min1)
    r_min1=r_dist4;
  if(r_dist5<r_min2)
    r_min2=r_dist5;
  if(r_dist6<r_min2)
    r_min2=r_dist6;

  double f_diff=f_min1-f_min2;
  if(f_diff<0) f_diff*=-1;
  double r_diff=r_min1-r_min2;
  if(r_diff<0) r_diff*=-1;
  if(f_diff<EUCLIDEAN_DISTANCE)
    f_primary_secondary=30;
  if(r_diff<EUCLIDEAN_DISTANCE)
    r_primary_secondary=3;

  return (f_primary_secondary+r_primary_secondary);
}

void TagClust::create_global_cluster() {
  sort(_f_cluster_list->begin(),_f_cluster_list->end(),cmp_cluster_center);
  sort(_r_cluster_list->begin(),_r_cluster_list->end(),cmp_cluster_center);

  vector<Cluster> f_cluster_list(*_f_cluster_list);
  vector<Cluster> r_cluster_list(*_r_cluster_list);
  for(int i=0;i<f_cluster_list.size();++i)
    sort(f_cluster_list[i]._point.begin(),f_cluster_list[i]._point.end(),cmp_tag_start);
  for(int i=0;i<r_cluster_list.size();++i)
    sort(r_cluster_list[i]._point.begin(),r_cluster_list[i]._point.end(),cmp_tag_start);

  sort(_f_cluster_list->begin(),_f_cluster_list->end(),cmp_cluster_strength);

  for(vector<Cluster>::reverse_iterator f_cluster_iter=_f_cluster_list->rbegin();f_cluster_iter!=_f_cluster_list->rend();++f_cluster_iter) {
    unsigned int f_peak=f_cluster_iter->_peak;
    bool flag=false;
    Cluster *r_cluster=0;

    for(vector<Cluster>::iterator r_cluster_iter=_r_cluster_list->begin();r_cluster_iter!=_r_cluster_list->end();++r_cluster_iter) { //start for
      unsigned int r_peak=r_cluster_iter->_peak;
      if(r_peak>(f_peak+SEARCH_UPPER_BOUND))
	break;
      if(r_cluster_iter->_id!=-1)
	continue;
      if(r_peak>(f_peak+SEARCH_LOWER_BOUND)) { //start if
	if(flag==false) {
	  r_cluster=&(*r_cluster_iter);
	  flag=true;
	}
	else {
	  if(cmp_cluster_strength(*r_cluster,*r_cluster_iter)==true) {
	    r_cluster=&(*r_cluster_iter);
	  }
	}
      } //end if
    } //end for

    if(flag==true) {
      int id=classify_id(f_cluster_list,r_cluster_list,f_cluster_iter->_peak,r_cluster->_peak,f_cluster_iter->_max_read,r_cluster->_max_read);

      if(id==11) { //primary,primary
	f_cluster_iter->_id=1;
	r_cluster->_id=1;
      }
      else if(id==12) { //primary,secondary
	f_cluster_iter->_id=1;
	r_cluster->_id=2;
      }
      else if(id==21) { //secondary,primary
	f_cluster_iter->_id=2;
	r_cluster->_id=1;
      }
      else if(id==22) { //secondary,secondary
	f_cluster_iter->_id=2;
	r_cluster->_id=2;
      }
      else if(id==31) {
	f_cluster_iter->_id=3;
	r_cluster->_id=1;
      }
      else if(id==32) {
	f_cluster_iter->_id=3;
	r_cluster->_id=2;
      }
      else if(id==13) {
	f_cluster_iter->_id=1;
	r_cluster->_id=3;
      }
      else if(id==23) {
	f_cluster_iter->_id=2;
	r_cluster->_id=3;
      }
      else if(id==33) {
	f_cluster_iter->_id=3;
	r_cluster->_id=3;
      }
      else {
	cout<<"nuc id error"<<endl;
	exit(0);
      }
    }
  }

  sort(_f_tag_list->begin(),_f_tag_list->end(),cmp_tag_start);
  sort(_r_tag_list->begin(),_r_tag_list->end(),cmp_tag_start);
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

void TagClust::clean_no_matching_cluster(vector<Cluster>* cluster_list) {
  for(vector<Cluster>::iterator cluster_iter=cluster_list->begin();cluster_iter!=cluster_list->end();) {
    if(cluster_iter->_id==-1)
      cluster_iter=cluster_list->erase(cluster_iter);
    else
      ++cluster_iter;
  }
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
  int original_cluster_count=cluster_list->size();
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
