#ifndef MN_TAG_CLUST_H
#define MN_TAG_CLUST_H

#include "../src_common/Option.h"
#include "../src_common/Struct.h"

#include <algorithm>
#include <fstream>
#include <vector>

namespace mn {
class TagClust {
 private:
  vector<Tag>* _f_tag_list;
  vector<Tag>* _r_tag_list;
  vector<Cluster>* _f_cluster_list;
  vector<Cluster>* _r_cluster_list;

  void cal_tag_score(double);
  void cal_statistic(vector<Cluster>*);
  double cal_score(vector<Cluster>*,vector<Cluster>*);
  void clean_no_matching_cluster(vector<Cluster>*,bool);
  void clean_noisy_cluster(vector<Cluster>*);
  void print_cluster(bool,vector<Cluster>*);
  void cal_corresponding_cluster(vector<Cluster>*,vector<Cluster>*,unsigned int);
  void maximize_cluster_slimiarity(int,vector<Cluster>*,vector<Cluster>*);
  unsigned int cal_center(vector<Tag>*,unsigned short);
  double cal_center_std(vector<unsigned int>*,unsigned int);
  int cal_quantile(vector<int>*,double);

 public:
  TagClust(vector<Tag>*,double);
  ~TagClust();
  void determine_center(bool);
  void create_local_cluster(bool);
  void create_global_cluster();
  void write_cluster(bool,const char*);
  vector<Cluster>* get_f_cluster_list();
  vector<Cluster>* get_r_cluster_list();
};
}

#endif //MN_TAG_CLUST_H
