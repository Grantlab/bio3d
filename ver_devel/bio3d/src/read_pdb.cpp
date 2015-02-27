// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
using namespace Rcpp;


// trim from start
string ltrim(std::string s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}

// trim from end
string rtrim(std::string s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}

// trim from both ends
string trim(std::string s) {
  return ltrim(rtrim(s));
}

// [[Rcpp::export('.read_pdb')]]
List read_pdb(std::string filename) {
  Rcpp::List out;
  
  vector<string> type;
  vector<int> eleno;
  vector<string> elety;
  vector<string> alt;
  vector<string> resid;
  vector<string> chain;
  vector<int> resno;
  vector<string> insert;
  vector<double> x;
  vector<double> y;
  vector<double> z;
  vector<double> o;
  vector<double> b;
  vector<string> segid;
  vector<string> elesy;
  vector<double> charge;
  vector<double> xyz;

  int models = 0;
  string tmp;
  string line;
  ifstream myfile;
  myfile.open(filename);
  
  if (myfile.is_open())  {
    while ( getline (myfile,line) ) {
      if(line.substr(0,5)=="MODEL") {
	models+=1;
      }
      
      if(line.substr(0,4)=="ATOM" || line.substr(0,5)=="HETATM") {
	// coordinates
	double tmpx = std::stod(line.substr(30,8));
	double tmpy = std::stod(line.substr(38,8));
	double tmpz = std::stod(line.substr(46,8));

	xyz.push_back(tmpx);
	xyz.push_back(tmpy);
	xyz.push_back(tmpz);

	
	if(models<2) {
	  x.push_back(tmpx);
	  y.push_back(tmpy);
	  z.push_back(tmpz);
	  
	  type.push_back(rtrim(line.substr(0,5)));
	  eleno.push_back(std::stoi(line.substr(6,5)));
	  elety.push_back(trim(line.substr(12,4)));
	  alt.push_back(trim(line.substr(16,1)));
	  resid.push_back(trim(line.substr(17,4)));
	  chain.push_back(trim(line.substr(21,1)));
	  resno.push_back(std::stoi(line.substr(22,4)));
	  insert.push_back(line.substr(26,1));
	
	  // 4 last entries
	  o.push_back(std::stod(line.substr(54,6)));
	  b.push_back(std::stod(line.substr(60,6)));
	  
	  // elesy
	  tmp = line.substr(76,2);
	  if(tmp=="  ")
	    elesy.push_back("");
	  else
	    elesy.push_back(tmp);
	  
	  // charge
	  tmp = line.substr(78,2);
	  if(tmp=="  ")
	    charge.push_back(NA_REAL);
	  else
	    charge.push_back(std::stod(tmp));
	}
      }
    }
    myfile.close();
  }
  else {
    out = Rcpp::List::create(Rcpp::Named("error")="error");
    return(out);
  }
  
  out = Rcpp::List::create(Rcpp::Named("atom")=
			   Rcpp::DataFrame::create(
						   Rcpp::Named("type")=type,
						   Rcpp::Named("eleno")=eleno,
						   Rcpp::Named("elety")=elety,
						   Rcpp::Named("alt")=alt,
						   Rcpp::Named("resid")=resid,
						   Rcpp::Named("chain")=chain,
						   Rcpp::Named("resno")=resno,
						   Rcpp::Named("insert")=insert,
						   Rcpp::Named("x")=x,
						   Rcpp::Named("y")=y,
						   Rcpp::Named("z")=z,
						   Rcpp::Named("o")=o,
						   Rcpp::Named("b")=b,
						   Rcpp::Named("elesy")=elesy,
						   Rcpp::Named("charge")=charge,
						   Rcpp::Named("stringsAsFactors")=false
						   ),
			   Rcpp::Named("xyz")=xyz,
			   Rcpp::Named("models")=models
			   );
				
				
  return(out);
}
