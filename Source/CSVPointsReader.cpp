#include "CSVPointReader.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdio.h>

using namespace std;

CSVPointsReader::CSVPointsReader(const char *fname){
  string str;
  m_line_count = 0;
  int start, end;
  char dl = ' ';
  ifstream file(fname);
  if (file.is_open()) {
    while (getline(file, str)) {
      m_line.push_back(str);
      m_line_count++;

      string delimiter = ",";

      size_t found = str.find(delimiter);
      double coords[3];
      int i =0;
      //print all tokens except last one 
      while (found != string::npos)
      {
      //create token
          string token = str.substr(0, found);
      //create a substring from the found position to end of the string
          str = str.substr(found + delimiter.length(), str.length() - 1);
          cout << token << '\n';
          found = str.find(delimiter);
          
          coords[i] =  strtod(token.c_str(),NULL);
          cout << "double "<<coords[i]<<endl;
          i++;
      }
      m_points.push_back(new Point(coords[0],coords[1],coords[2]));

    }
  }   
}