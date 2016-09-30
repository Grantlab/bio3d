#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <utility>
#include "utils.h"
#include "convert.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::export('.read_crd')]]
List read_crd(std::string filename) {

  // out is a List object
  Rcpp::List out = Rcpp::List::create();
  StringVector content_names;

  // vectors
  vector<double> coords;
  vector<double> vels;
  vector<double> box;

  // data
  string title;
  int natoms = NA_INTEGER;
  double simtime = NA_REAL;
  
  // temp variables
  int i = 0;
  string tmps;
  double tmpd;
  vector<string> tmp_vec;
    
  // for reading
  string line;
  ifstream myfile;
  
  // open file and iterate over each line
  myfile.open(filename.c_str());
  
  if (myfile.is_open())  {
    while ( getline (myfile,line) ) {
      i++;
      
      // title line
      if(i==1) {
	title=trim(line);
	continue;
      }
      
      // natoms and simtime
      if(i==2) {
	tmps = trim(line);
	tmp_vec = split(tmps, ' ');
	
	if(tmp_vec.size()>0) {
	  natoms = stringToInt(tmp_vec[0]);
	}
	
	if(tmp_vec.size()>1) {
	  simtime = stringToDouble(tmp_vec[1]);
	}

	continue;
      }
      
      // coords
      if(natoms*3 > (int)coords.size()) {
	stringstream ss(line);
	
	while (ss >> std::setw(12) >> tmpd) {
	  coords.push_back(tmpd);
	}
	
	continue;
      }

      // velocities
      if((int)coords.size() == natoms*3 && (int)vels.size() < natoms*3) {
	stringstream ss(line);
	
	while (ss >> std::setw(12) >> tmpd) {
	  vels.push_back(tmpd);
	}
	
	continue;
      }

      // box size
      if((int)coords.size() == natoms*3) {
	stringstream ss(line);
	
	while (ss >> std::setw(12) >> tmpd) {
	  box.push_back(tmpd);
	}
	
	continue;
      }
    }
    myfile.close();

    
  }
  else {
    out = Rcpp::List::create(Rcpp::Named("error")="error");
    return(out);
  }
  
  // if velocities were not present, box info might be stored there
  if(vels.size()==6 && box.empty()) {
    box=vels;
    vels.clear();
  }
  
  out.push_back(title);
  out.push_back(natoms);
  out.push_back(simtime);
  out.push_back(coords);
  out.push_back(vels);
  out.push_back(box);
  
  content_names.push_back("title");
  content_names.push_back("natoms");
  content_names.push_back("time");
  content_names.push_back("xyz");
  content_names.push_back("velocities");
  content_names.push_back("box");
   
  out.attr("names") = content_names;
  return(out);
}
