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

// function get hexadecimal
int getHex(string hexstr) {
  return (int)strtol(hexstr.c_str(), 0, 16);
}

// [[Rcpp::export('.read_pdb')]]
List read_pdb(std::string filename, bool multi, bool hex) {
  // out is a List object
  Rcpp::List out;

  // keep track of number of atoms and models in PDB
  int natoms = 0;
  int models = 0;
  
  // assign vectors for building final 'atom' object
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
  vector<string> charge;

  // xyz object
  vector<double> xyz;
  
  // store HELIX / SHEET records
  vector<string> helix_chain;
  vector<int> helix_resno_start;
  vector<int> helix_resno_end;
  vector<string> helix_type;

  vector<string> sheet_chain;
  vector<int> sheet_resno_start;
  vector<int> sheet_resno_end;
  vector<string> sheet_sense;

  // store SEQRES
  vector<string> seqres;
  vector<string> seqres_chain;

  // temp variables
  string tmp;
  int tmp_eleno;

  // for reading
  string line;
  ifstream myfile;

  // open file and iterate over each line
  myfile.open(filename);
  
  if (myfile.is_open())  {
    while ( getline (myfile,line) ) {
      
      // keep of track of number of models in PDB file
      if(line.substr(0,5)=="MODEL") {
	models+=1;
	
	// break out of loop if we dont want multi-model
	if(!multi && models > 1) {
	  models=1;
	  break;
	}
      }
      
      // store helix info
      else if(line.substr(0,5)=="HELIX") {
	helix_chain.push_back(trim(line.substr(19,1)));
	helix_resno_start.push_back(std::stoi(line.substr(21,4)));
	helix_resno_end.push_back(std::stoi(line.substr(33,4)));
	helix_type.push_back(trim(line.substr(38,2)));
      }

      // store sheet info
      else if(line.substr(0,5)=="SHEET") {
	sheet_chain.push_back(trim(line.substr(21,1)));
	sheet_resno_start.push_back(std::stoi(line.substr(22,4)));
	sheet_resno_end.push_back(std::stoi(line.substr(33,4)));
	sheet_sense.push_back(trim(line.substr(38,2)));
      }
      
      // store SEQRES info
      else if(line.substr(0,6)=="SEQRES") {
	for(int i=0; i<13; i++) { 
	  tmp = trim(line.substr(19+(i*4),3));
	  if(tmp != "") {
	    seqres.push_back(tmp);
	    seqres_chain.push_back(trim(line.substr(11,1)));
	  }
	}
      }
      
      // store ATOM/HETATM records
      else if(line.substr(0,4)=="ATOM" || line.substr(0,5)=="HETATM") {
	// read coordinates
	double tmpx = std::stod(line.substr(30,8));
	double tmpy = std::stod(line.substr(38,8));
	double tmpz = std::stod(line.substr(46,8));
	
	// always store coords in xyz object
	xyz.push_back(tmpx);
	xyz.push_back(tmpy);
	xyz.push_back(tmpz);
	
	// only store other items for first MODEL
	if(models < 2) {
	  natoms++;
	  
	  // x, y, z for 'atom'
	  x.push_back(tmpx);
	  y.push_back(tmpy);
	  z.push_back(tmpz);

	  // eleno can be hexadecimal (e.g. from VMD)
	  if(hex && natoms > 99999) {
	    tmp_eleno = getHex(trim(line.substr(6,5)));
	  }
	  else {
	    tmp_eleno = std::stoi(line.substr(6,5));
	  }
	  eleno.push_back(tmp_eleno);
	  
	  // read all others items as they are
	  resno.push_back(std::stoi(line.substr(22,4)));
	  type.push_back(rtrim(line.substr(0,5)));
	  elety.push_back(trim(line.substr(12,4)));
	  alt.push_back(trim(line.substr(16,1)));
	  resid.push_back(trim(line.substr(17,4)));
	  chain.push_back(trim(line.substr(21,1)));
	  insert.push_back(line.substr(26,1));
	  o.push_back(std::stod(line.substr(54,6)));
	  b.push_back(std::stod(line.substr(60,6)));
	  segid.push_back(trim(line.substr(72,4)));
	  elesy.push_back(trim(line.substr(76,2)));
	  charge.push_back(trim(line.substr(78,2)));
	}
      }
    }
    myfile.close();
  }
  else {
    out = Rcpp::List::create(Rcpp::Named("error")="error");
    return(out);
  }
  
  // build output List
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
						   Rcpp::Named("segid")=segid,
						   Rcpp::Named("elesy")=elesy,
						   Rcpp::Named("charge")=charge,
						   Rcpp::Named("stringsAsFactors")=false
						   ),
			   Rcpp::Named("xyz")=xyz,
			   Rcpp::Named("models")=models,
			   Rcpp::Named("seqres")=seqres,
			   Rcpp::Named("seqres_chain")=seqres_chain,

			   Rcpp::Named("helix")=
			   Rcpp::List::create(
					      Rcpp::Named("start")=helix_resno_start,
					      Rcpp::Named("end")=helix_resno_end,
					      Rcpp::Named("chain")=helix_chain, 
					      Rcpp::Named("chain")=helix_type
					      ),
			   
			   Rcpp::Named("sheet")=
			   Rcpp::List::create(
					      Rcpp::Named("start")=sheet_resno_start,
					      Rcpp::Named("end")=sheet_resno_end,
					      Rcpp::Named("chain")=sheet_chain,
					      Rcpp::Named("sense")=sheet_sense
					      )
			   
			   );
				
				
  return(out);
}
