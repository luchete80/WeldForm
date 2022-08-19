#ifndef NASTRAN_VOL_READER_H_
#define NASTRAN_VOL_READER_H_

#include "NastranReader.h"
#include "matvec.h"

namespace SPH {

class NASElement{
  int nodecount;
  std::vector <int> node;
  Vec3_t 	centroid;	
  
};

class Hexa: public NASElement{
  Hexa(std::vector<int> nodes){}
};

class Prism: public NASElement{
  Prism(std::vector<int> nodes){}
};

class Pyramdid: public NASElement{
  Pyramdid(std::vector<int> nodes){}
};

using namespace std;

class NastranVolReader:
public NastranReader {
  int line_start_elem;
  std::vector <NASElement*> elem;
  public:
  NastranVolReader(char* fName){read(fName);}
	inline void read(const char *fName);
  };
};

#include "NastranVolReader.cpp"

#endif