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

struct dataformat {
  string fmt;
  string type;
  int length;
  int width;
};
  
dataformat getFormatFromLine(std::string line) {
  dataformat out;
  string tmp;
  string fmt;
  string type;
  int length;
  int width;
  size_t type_pos, tmplast;

  // format e.g. (5E16.8)
  fmt = trim(line.substr(7));

  // position of type definer
  type_pos = fmt.find_first_of("aEI");
  
  // fetch the type character
  type = fmt[type_pos];

  // length is imediately before the type char
  length = stringToInt(fmt.substr(1, type_pos-1));
  
  // width is after type char
  tmp = fmt.substr(type_pos+1);
  tmplast = tmp.find_first_of(".)");
  width = stringToInt(tmp.substr(0, tmplast));
  
  out.fmt=fmt;
  out.type=type;
  out.length=length;
  out.width=width;
  
  return(out);
}

// [[Rcpp::export('.read_prmtop')]]
List read_prmtop(std::string filename) {

  // out is a List object
  Rcpp::List out = Rcpp::List::create();
  StringVector content_names;

  // temporary vectors
  vector< vector<string> > str_data;
  vector< vector<int> > int_data;
  vector< vector<double> > dbl_data;
  
  // temp vectors
  vector<string> tmp_svector;
  vector<int> tmp_ivector;
  vector<double> tmp_dvector;
  
  // temp variables
  int tmpi;
  string tmps;
  double tmpd;
  
  // keep track of current fields
  dataformat fmt;
  string current_flag;
  string previous_flag = "";
  string current_format;
  string current_type;
  int current_width = 0;
  
  // for reading
  string line;
  ifstream myfile;
  
  // open file and iterate over each line
  myfile.open(filename.c_str());
  
  if (myfile.is_open())  {
    while ( getline (myfile,line) ) {
      line = trim(line);
      
      if(line.substr(0, 5)=="%FLAG") {
	current_flag=trim(line.substr(6));
	content_names.push_back(current_flag);
	
	// store the data
	if(tmp_svector.size()>0) {
	  out.push_back(tmp_svector);
	  tmp_svector.clear();
	}

	if(tmp_dvector.size()>0) {
	  out.push_back(tmp_dvector);
	  tmp_dvector.clear();
	}

	if(tmp_ivector.size()>0) {
	  out.push_back(tmp_ivector);
	  tmp_ivector.clear();
	}
	  
	continue;
      }

      else if(line.substr(0,7)=="%FORMAT"){
	fmt = getFormatFromLine(line);
	current_type=fmt.type;
	current_width=fmt.width;
	continue;	
      }
      
      else {
	stringstream ss(line);
	
	// string
	if(current_type=="a"){
	  while (ss >> std::setw(current_width) >> tmps) {
	    tmp_svector.push_back(tmps);
	  }

	  if(tmp_svector.size()==0)
	    tmp_svector.push_back("");
	  continue;
	}
	
	// double
	if(current_type=="E"){
	  while (ss >> std::setw(current_width) >> tmpd) {
	    tmp_dvector.push_back(tmpd);
	  }

	  if(tmp_dvector.size()==0)
	    tmp_dvector.push_back(NA_REAL);
	  continue;
	}
	
	// integers
	if(current_type=="I"){
	  while (ss >> std::setw(current_width) >> tmpi) {
	    tmp_ivector.push_back(tmpi);
	  }

	  if(tmp_ivector.size()==0)
	    tmp_ivector.push_back(NA_INTEGER);
	  continue;
	}
      }
    }
    myfile.close();

    // store last round of data
    if(tmp_svector.size()>0) {
      out.push_back(tmp_svector);
      tmp_svector.clear();
    }
    
    if(tmp_dvector.size()>0) {
      out.push_back(tmp_dvector);
      tmp_dvector.clear();
    }
    
    if(tmp_ivector.size()>0) {
      out.push_back(tmp_ivector);
      tmp_ivector.clear();
    }
  }
  else {
    out = Rcpp::List::create(Rcpp::Named("error")="error");
    return(out);
  }
  
  out.attr("names") = content_names;
  return(out);
}
