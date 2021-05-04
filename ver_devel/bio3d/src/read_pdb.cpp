#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>

#include "gzstream.h"
#include "utils.h"
#include "convert.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::export('.read_pdb')]]
List read_pdb(std::string filename, bool multi=false, bool hex=false, int maxlines=-1,
	      bool atoms_only=false) {
  // out is a List object
  Rcpp::List out;

  // keep track of number of atoms and models in PDB
  int natoms = 0;
  int models = 0;
  int nresi = 0;

  // reset this when reading MODEL, END, ENDMDL
  // bump models when this is 0
  int natoms_model = 0;
  
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
  vector<string> helix_inserti;
  vector<string> helix_inserte;

  vector<string> sheet_chain;
  vector<int> sheet_resno_start;
  vector<int> sheet_resno_end;
  vector<string> sheet_sense;
  vector<string> sheet_inserti;
  vector<string> sheet_inserte;

  // store SEQRES
  vector<string> seqres;
  vector<string> seqres_chain;

  // store REMARK 350 lines
  vector<string> remark350;

  // store HEADER
  string header;
  
  // temp variables
  string tmp;
  int tmp_eleno;
  int tmp_resno;

  // counter for seqres while loop
  int i;

  // values from prev 
  string curr_record_str;
  string prev_eleno_str;
  string curr_eleno_str;
  string prev_resno_str;
  string curr_resno_str;
  string prev_chain;
  string curr_chain;
  int counter = 0;

  bool hex_eleno = false;
  bool hex_resno = false;

  // for reading
  vector<string> raw_lines;
  string line;
  
  igzstream mystream;
  mystream.open(filename.c_str());
  
  
  if (mystream.is_open())  {

    // read line by line
    while ( getline (mystream, line) ) {
  
      // break if user has provided maxlines argument
      counter++;
      if(maxlines != -1 && counter > maxlines) {
	break;
      }

      curr_record_str = trim(line.substr(0,6));
      //cout << curr_record_str;

      // reset 'natoms_model' when reading END, ENDMDL or MODEL
      if(curr_record_str == "END" ||
	 curr_record_str == "ENDMDL" ||
	 curr_record_str == "MODEL") {

	natoms_model = 0;
      }
      
      // output header only if verbose=true
      else if(curr_record_str=="HEADER") {
	header = line;
      }

      // store helix info
      else if(curr_record_str=="HELIX" && !atoms_only) {
	helix_chain.push_back(trim(line.substr(19,1)));
	helix_resno_start.push_back(stringToInt(line.substr(21,4)));
	helix_resno_end.push_back(stringToInt(line.substr(33,4)));
	//helix_resno_start.push_back(std::stoi(line.substr(21,4)));
	//helix_resno_end.push_back(std::stoi(line.substr(33,4)));
	helix_type.push_back(trim(line.substr(38,2)));
	helix_inserti.push_back(trim(line.substr(25,1)));
	helix_inserte.push_back(trim(line.substr(37,1)));
      }
      
      // store sheet info
      else if(curr_record_str=="SHEET" && !atoms_only) {
	sheet_chain.push_back(trim(line.substr(21,1)));
	sheet_resno_start.push_back(stringToInt(line.substr(22,4)));
	sheet_resno_end.push_back(stringToInt(line.substr(33,4)));
	//sheet_resno_start.push_back(std::stoi(line.substr(22,4)));
	//sheet_resno_end.push_back(std::stoi(line.substr(33,4)));
	sheet_sense.push_back(trim(line.substr(38,2)));
	sheet_inserti.push_back(trim(line.substr(26,1)));
	sheet_inserte.push_back(trim(line.substr(37,1)));
      }
      
      // store SEQRES info
      else if(curr_record_str=="SEQRES" && !atoms_only) {

	i=0;
	while(19+(i*4) < (int) line.length()) {
	  tmp = trim(line.substr(19+(i*4), 3));

	  if(tmp != "") {
	    seqres.push_back(tmp);
	    seqres_chain.push_back(trim(line.substr(11,1)));
	  }
	  i++;
	}
      }

      // store REMARK 350 records (raw lines only)
      else if(line.substr(0,10)=="REMARK 350" && !atoms_only) {
	remark350.push_back(trim(line));
      }
      
      // store ATOM/HETATM records
      else if(curr_record_str=="ATOM" || curr_record_str=="HETATM") {
	// number of atoms in current model
	natoms_model++;

	// bump 'models' when reading first atom in model
	if(natoms_model == 1) {
	  models++;
	}

	// break out of loop when visiting 2nd model
	if(multi == false && models > 1) {
	  models = 1;
	  break;
	}
	
	// read coordinates
	double tmpx = stringToDouble(line.substr(30,8));
	double tmpy = stringToDouble(line.substr(38,8));
	double tmpz = stringToDouble(line.substr(46,8));
	
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

	  // read atom number
	  curr_eleno_str = trim(line.substr(6,5));

	  // enable parsing hexadecimals (e.g. from VMD)
	  // issue warning if hexdec found
	  if(natoms > 99999 && !hex_eleno && 
	     curr_eleno_str == "186a0" && prev_eleno_str == "99999") {
	    
	    if(!hex) {
	      Rf_warning("possible hexadecimal atom numbering. re-run with 'hex=TRUE'");
	    }
	    else {
	      Rf_warning("hexadecimal atom numbering used for eleno > 99.999");
	    }
	    hex_eleno = true;
	  }
	    
	  if(hex && hex_eleno) {
	    tmp_eleno = getHex(curr_eleno_str);
	  }
	  else {
	    try {
	      tmp_eleno = stringToInt(curr_eleno_str);
	    }
	    catch (...) {
	      tmp_eleno = NA_INTEGER;
	    }
	  }
	  prev_eleno_str = curr_eleno_str;


	  // read residue number and chain identifier
	  curr_resno_str = trim(line.substr(22,4));
	  curr_chain = trim(line.substr(21,1));

	  // count number of residues
	  if(curr_resno_str != prev_resno_str) {
	    // reset nresi if chain change
	    if(curr_chain != prev_chain) {
	      nresi = 0;
	      hex_resno = false;
	    }
	    prev_chain = curr_chain;
	    nresi++;
	  }

	  // enable parsing hexadecimals (e.g. from VMD)
	  // issue warning if hexdec found
	  if(nresi > 9999 && !hex_resno &&
	     curr_resno_str == "2710" && prev_resno_str == "9999") {

	    if(!hex) {
	      Rf_warning("possible hexadecimal residue numbering. re-run with 'hex=TRUE'");
	    }
	    else {
	      Rf_warning("hexadecimal residue numbering used for resno > 9.999");
	    }
	    hex_resno = true;
	  }
	  
	  if(hex && hex_resno) {
	    tmp_resno = getHex(curr_resno_str);
	  }
	  else {
	    try {
	      tmp_resno = stringToInt(curr_resno_str);
	    }
	    catch (...) {
	      tmp_resno = NA_INTEGER;
	    }
	  }
	  prev_resno_str = curr_resno_str;
	  

	  // store eleno resno chain to arrays
	  eleno.push_back(tmp_eleno);
	  resno.push_back(tmp_resno);
	  chain.push_back(curr_chain);

	  // read all others items as they are
	  type.push_back(rtrim(line.substr(0,6)));
	  elety.push_back(trim(line.substr(12,4)));
	  alt.push_back(trim(line.substr(16,1)));
	  resid.push_back(trim(line.substr(17,4)));
	  insert.push_back(trim(line.substr(26,1)));

	  try{
	     o.push_back(stringToDouble(trim(line.substr(54,6))));
	  } catch(...) {
	     o.push_back(NA_REAL);
	  }

	  try{
	     b.push_back(stringToDouble(trim(line.substr(60,6))));
	  } catch(...) {
	     b.push_back(NA_REAL);
	  }
	  
	  try{
	    segid.push_back(trim(trim(line.substr(72,4))));
	  } catch(...) {
	    segid.push_back("");
	  }
	  
	  try{
	    elesy.push_back(trim(line.substr(76,2)));
	  } catch(...) {
	    elesy.push_back("");
	  }
	  
	  try{
	    charge.push_back(trim(line.substr(78,2)));
	  } catch(...) {
	    charge.push_back("");
	  }
	}
      }
    } // while ( getline (myfile,line) ) {
    mystream.close();
  } // end   if (myfile.is_open())  {

  else {
    out = Rcpp::List::create(Rcpp::Named("error")="Error when reading file");
    return(out);
  }

  // add names to helix / sheet vector
  NumericVector helix_start_out = wrap(helix_resno_start);
  helix_start_out.names() = helix_inserti;
  NumericVector helix_end_out = wrap(helix_resno_end);
  helix_end_out.names() = helix_inserte;

  NumericVector sheet_start_out = wrap(sheet_resno_start);
  sheet_start_out.names() = sheet_inserti;
  NumericVector sheet_end_out = wrap(sheet_resno_end);
  sheet_end_out.names() = sheet_inserte;
  
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
					      Rcpp::Named("start")=helix_start_out,
					      Rcpp::Named("end")=helix_end_out,
					      Rcpp::Named("chain")=helix_chain, 
					      Rcpp::Named("type")=helix_type
					      ),
			   
			   Rcpp::Named("sheet")=
			   Rcpp::List::create(
					      Rcpp::Named("start")=sheet_start_out,
					      Rcpp::Named("end")=sheet_end_out,
					      Rcpp::Named("chain")=sheet_chain,
					      Rcpp::Named("sense")=sheet_sense,
					      Rcpp::Named("inserti")=sheet_inserti
					      ),

			   Rcpp::Named("remark350")=remark350,

			   Rcpp::Named("header")=header
			   
			   );
  
  return(out);
}
