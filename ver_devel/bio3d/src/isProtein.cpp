// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// [[Rcpp::export('.isProteinCpp')]]
LogicalVector isProtein(NumericVector resno, 
			CharacterVector chain, 
			CharacterVector insert, 
			CharacterVector elety) {
  
  // atoms in input
  int n = resno.size();

  // output (false by default)
  LogicalVector out(n);
  
  // atoms seen in current residue
  std::unordered_set<int> seen;

  // indices for the atoms in current residue
  std::vector<int> inds;

  // track when to jump to next residue 
  std::string prev_chain = "null";
  std::string prev_insert = "null";
  int prev_res = -99999;

  // iterate over the atoms
  for( int i = 0; i < n; i++ ) {
    
    if(resno[i] != prev_res || 
       Rcpp::as<std::string>(chain[i]) != prev_chain || 
       Rcpp::as<std::string>(insert[i]) != prev_insert ) {
      
      if(i>0) {
	int m = inds.size();
	int l = seen.size();
	if(l==4) {
	  for( int j=0; j<m; j++) {
	    out[inds[j]]=true;
	  }
	}
	
	inds.erase(inds.begin(), inds.end());
	seen.erase(seen.begin(), seen.end());
      }
    }

    if(elety[i]=="CA")
      seen.insert(1);
    if(elety[i]=="C")
      seen.insert(2);
    if(elety[i]=="N")
      seen.insert(3);
    if(elety[i]=="O")
      seen.insert(4);

    inds.push_back(i);    
    prev_res = resno[i];
    prev_chain = chain[i];
    prev_insert = insert[i];
  }
  
  return out;
}
