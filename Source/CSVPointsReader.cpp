#include "CSVPointReader.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdio.h>

using namespace std;

CSVPointsReader::CSVPointsReader(const char *fname){
  string line;
  m_line_count = 0;
  int start, end;
  char dl = ' ';
  ifstream file(fname);
  if (file.is_open()) {
    while (getline(file, line)) {
      m_line.push_back(line);
      m_line_count++;
    }
  }   
}