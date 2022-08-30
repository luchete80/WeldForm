#ifndef NASTRAN_VOL_READER_H_
#define NASTRAN_VOL_READER_H_

#include "NastranReader.h"
#include "matvec.h"

namespace SPH {

class NastranVolReader;

class NASElement{
  public:
  int nodecount;
  double vol;
  NastranVolReader *mesh;
  std::vector <int> node;
  int id;
  Vec3_t 	centroid;	
  void allocate(std::vector<int>nid){
    node.resize(nodecount);
    for (int i=0;i<nodecount;i++)node[i]=nid[i];
  }
  
  void CalcCentroids();
  virtual void CalcVol(){};
  virtual ~NASElement(){}
};

class Hexa: public NASElement{
  public:
  Hexa(NastranVolReader *m, std::vector<int> nodes){
    mesh = m;
    nodecount=8;
    if (nodes.size()!=8) cout << "CHEXA Node Count ERROR"<<endl;
    allocate(nodes);
    CalcCentroids();
    CalcVol();
  }
  void CalcVol();
};

class Tetra: public NASElement{
  public:
  Tetra(NastranVolReader *m, std::vector<int> nodes){
     mesh = m;
    nodecount=4;
    if (nodes.size()!=4) cout << "CTETRA Node Count ERROR"<<endl;
    allocate(nodes);
    CalcCentroids();
    CalcVol();
  }
  void CalcVol();
  ~Tetra(){}
};


class Prism: public NASElement{
  public:
  Prism(NastranVolReader *m, std::vector<int> nodes){
     mesh = m;
    nodecount=6;
    if (nodes.size()!=6) cout << "CPRISM Node Count ERROR"<<endl;
    allocate(nodes);
    CalcCentroids();
    CalcVol();
  }
  void CalcVol(){
    vol = 0.;
  }
};

class Pyram: public NASElement{
  public:
  Pyram(NastranVolReader *m, std::vector<int> nodes){
     mesh = m;
    nodecount=5;
    if (nodes.size()!=5) cout << "CPYRAM Node Count ERROR"<<endl;
    allocate(nodes);
    CalcCentroids();
    CalcVol();
  }
  void CalcVol();
};

using namespace std;

class NastranVolReader:
public NastranReader {
  public: 
  int line_start_elem;
  std::vector<Vec3_t* > nod;
  std::vector <NASElement*> elem;
  public:
  NastranVolReader(char* fName){
    read(fName);
    double vol = 0.;
    for (int i=0;i<elem.size();i++)
      vol +=elem[i]->vol;
    cout << "Total vol: "<< vol << endl;
  }
	inline void read(const char *fName);
  ~NastranVolReader(){
    for (int i=0;i<nod.size();i++)
      delete nod[i];
    for (int i=0;i<elem.size();i++)
      delete elem[i];
    nod.clear();
    elem.clear();
  }
  };
};



#include "NastranVolReader.cpp"

#endif