#ifndef _CSV_READER_
#define _CSV_READER_

#include <vector>
#include <string>

struct Point {
public:
  Point(const double &x,const double &y,const double &z) {
    coords[0] = x;coords[1] = y;coords[2] = z;
  }
  double coords[3];
};

class CSVPointsReader{
public:
  CSVPointsReader (const char *);
  Point* getPoint(const int &i){return points[i];} 

protected:
  std::vector <Point*> points;
  int m_line_count;
  std::vector <std::string> m_line;  
};


#endif