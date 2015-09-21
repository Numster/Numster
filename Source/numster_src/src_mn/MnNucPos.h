#ifndef MN_NUC_POS_H
#define MN_NUC_POS_H

#include "../src_common/Option.h"
#include "../src_common/Struct.h"

#include <algorithm>
#include <fstream>
#include <vector>

namespace mn {
class NucPos {
 private:
  vector<Nuc>* _nuc_list;
  double _fuzzy_rate;
  void define_maximal_nucleosome(vector<Nuc>*,vector<Nuc>*);
  void define_weighted_nucleosome(vector<Nuc>*,vector<Nuc>*);
  double cal_center_std(vector<unsigned int>*,unsigned int);

 public:
  NucPos(unsigned short,vector<Cluster>*,vector<Cluster>*,double);
  ~NucPos();
  void position_nucleosome();
  void write_nucleosome(const char*);
  vector<Nuc>* get_nucleosome();
  int get_nuc_count();
  double get_fuzzy_rate();
};
}

#endif //MN_NUC_POS_H
