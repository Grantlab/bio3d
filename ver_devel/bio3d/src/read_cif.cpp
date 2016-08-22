#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>

#include "gzstream.h"
#include "utils.h"
#include "convert.h"
using namespace std;
using namespace Rcpp;

// [[Rcpp::export('.read_cif')]]
List read_cif(std::string filename, bool multi=false, bool hex=false, int maxlines=-1, bool atoms_only=false) {
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
  vector<int> entid;
  vector<string> insert;
  vector<double> x;
  vector<double> y;
  vector<double> z;
  vector<double> o;
  vector<double> b;
  vector<string> segid;
  vector<string> elesy;
  vector<string> charge;
  vector<int> model;

  // xyz object
  vector<double> xyz;
  
  // store HELIX / SHEET records
  /**
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
  */

  // store SEQRES
  vector<string> seqres;
  vector<string> seqres_chain;

  // store REMARK 350 lines
  vector<string> remark350;

  // temp variables
  string tmp;
  int counter = 0;

  string tmps;
  vector<string> tmp_vec;

  // for reading
  vector<string> raw_lines;
  string line;

  igzstream mystream;
  mystream.open(filename.c_str());
  
  if (mystream.is_open())  {
    while ( getline (mystream, line) ) {
  
      // break if user has provided maxlines argument
      counter++;
      
      if(maxlines != -1 && counter > maxlines) {
	break;
      }

      // should be done reading _loop lines instead I guess
      /**
      else if(line.substr(0,6)=="HELX_P") {
	
	tmps = trim(line);
	tmp_vec = split(tmps, ' ');
	
	helix_chain.push_back(trim(tmp_vec[4]));
	helix_resno_start.push_back(stringToInt(tmp_vec[5]));
	helix_resno_end.push_back(stringToInt(tmp_vec[9]));
	helix_type.push_back(trim(tmp_vec[17]));
	helix_inserti.push_back(trim(tmp_vec[6]));
	helix_inserte.push_back(trim(tmp_vec[10]));
	
      }
      */

      // store ATOM/HETATM records
      else if(line.substr(0,4)=="ATOM" || line.substr(0,6)=="HETATM") {
	
	tmps = trim(line);
	tmp_vec = split(tmps, ' ');
	//cout << tmps;
	//cout << "\n";

	// read coordinates
	double tmpx = stringToDouble(tmp_vec[10]);
	double tmpy = stringToDouble(tmp_vec[11]);
	double tmpz = stringToDouble(tmp_vec[12]);
	
	// always store coords in xyz object
	xyz.push_back(tmpx);
	xyz.push_back(tmpy);
	xyz.push_back(tmpz);

	// include here check for model number as in read_pdb.cpp
	// add checks for "?"
	models = stringToInt(tmp_vec[25]);
	  
	// only store other items for first MODEL
	if(models == 1) {
	  natoms++;
	
	  // x, y, z for 'atom'
	  x.push_back(tmpx);
	  y.push_back(tmpy);
	  z.push_back(tmpz);

	  //_atom_site.group_PDB
	  type.push_back(tmp_vec[0]);
	  
	  //_atom_site.id
	  eleno.push_back(stringToInt(tmp_vec[1]));
	
	  //_atom_site.type_symbol
	  elesy.push_back(tmp_vec[2]);
	  
	  //_atom_site.label_atom_id
	  //_atom_site.label_alt_id
	  alt.push_back(tmp_vec[4]);
	  
	  //_atom_site.label_comp_id
	  //_atom_site.label_asym_id
	  //_atom_site.label_entity_id
	  entid.push_back(stringToInt(tmp_vec[7]));
	  
	  //_atom_site.label_seq_id
	  //_atom_site.pdbx_PDB_ins_code
	  insert.push_back(tmp_vec[9]);
	  
	  //_atom_site.Cartn_x
	  //_atom_site.Cartn_y
	  //_atom_site.Cartn_z
	  //_atom_site.occupancy
	  o.push_back(stringToDouble(tmp_vec[13]));
	  //_atom_site.B_iso_or_equiv
	  b.push_back(stringToDouble(tmp_vec[14]));
	  
	  //_atom_site.Cartn_x_esd
	  //_atom_site.Cartn_y_esd
	  //_atom_site.Cartn_z_esd
	  //_atom_site.occupancy_esd
	  //_atom_site.B_iso_or_equiv_esd
	  //_atom_site.pdbx_formal_charge
	  charge.push_back(tmp_vec[20]);
	  
	  //_atom_site.auth_seq_id
	  resno.push_back(stringToInt(tmp_vec[21]));
	  //_atom_site.auth_comp_id
	  resid.push_back(tmp_vec[22]);
	  //_atom_site.auth_asym_id
	  chain.push_back(tmp_vec[23]);
	  //_atom_site.auth_atom_id
	  elety.push_back(tmp_vec[24]);
	  //_atom_site.pdbx_PDB_model_num
	  model.push_back(stringToInt(tmp_vec[25]));
	  
	} // if(models < 2) {
      } // if(line.substr(0,4)=="ATOM" 
    } // while ( getline (myfile,line) ) {
    
    mystream.close();
  } // end   if (myfile.is_open())  {
  
  else {
    out = Rcpp::List::create(Rcpp::Named("error")="Error when reading file");
    return(out);
  }

  // add names to helix / sheet vector
  /**
  NumericVector helix_start_out = wrap(helix_resno_start);
  helix_start_out.names() = helix_inserti;
  NumericVector helix_end_out = wrap(helix_resno_end);
  helix_end_out.names() = helix_inserte;
  */  

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
						   Rcpp::Named("entid")=entid,
						   Rcpp::Named("elesy")=elesy,
						   Rcpp::Named("charge")=charge,
						   Rcpp::Named("stringsAsFactors")=false
						   ),
			   /**
			   Rcpp::Named("helix")=
			   Rcpp::List::create(
			   		      Rcpp::Named("start")=helix_start_out,
			   		      Rcpp::Named("end")=helix_end_out,
			   		      Rcpp::Named("chain")=helix_chain, 
			   		      Rcpp::Named("type")=helix_type
			   		      ),
			   */
			   
			   Rcpp::Named("xyz")=xyz,
			   Rcpp::Named("models")=models
			   
			   );
  
  return(out);
}
