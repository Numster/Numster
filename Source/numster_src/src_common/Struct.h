#ifndef STRUCT_H
#define STRUCT_H

#include <iostream>
#include <vector>

using namespace std;

struct Tag {
  unsigned short _chromosome;
  unsigned int _start;
  char _strand;
  unsigned short _read;
  double _strength;
};

struct Cluster {
  int _id;
  unsigned int _center;
  unsigned int _peak;
  unsigned int _start;
  unsigned int _end;
  unsigned short _max_read;
  unsigned short _left_read;
  unsigned short _right_read;
  unsigned int _strength;
  vector<Tag> _point;
};

struct Nuc {
  unsigned short _chromosome;
  unsigned int _start;
  unsigned int _end;
  unsigned int _f_strength;
  unsigned int _r_strength;
  unsigned int _strength;
  double _occupancy;
  bool _overlap;
  bool _paired;
  double _std;
  double _f_occupancy;
  double _r_occupancy;
};

#endif //STRUCT_H
