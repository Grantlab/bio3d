// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// [[Rcpp::export('.isProteinCpp')]]
LogicalVector isProtein(NumericVector resno, CharacterVector chain, 
			CharacterVector insert, CharacterVector elety) {
  
  // atoms in input
  int n = resno.size();

  // output (false by default)
  LogicalVector out(n);
  
  // atoms seen in current residue
  std::unordered_set<int> seen;

  // indices for the atoms in current residue
  std::vector<int> inds;

  // true when last atom in residue
  bool last_resatom = false;

  // iterate over the atoms
  for( int i = 0; i < n; i++ ) {

    // if new residue
    if(last_resatom) {
      // empty vectors for tracking
      inds.erase(inds.begin(), inds.end());
      seen.erase(seen.begin(), seen.end());
    }
    
     // add atom to this residue
    inds.push_back(i);  
    
    // add backbone atoms 
    if(elety[i]=="CA")
      seen.insert(1);
    if(elety[i]=="C")
      seen.insert(2);
    if(elety[i]=="N")
      seen.insert(3);
    if(elety[i]=="O")
      seen.insert(4);
    
    // check if current atom is the last in current residue
    last_resatom = true;
    if(i<(n-1)) {
      if( resno[i] != resno[i+1] || chain[i] != chain[i+1] || insert[i] != insert[i+1] ) {
	last_resatom = true;
      }
      else {
        last_resatom= false;
      }
    }
    
    // mark atoms in current residue as protein
    int l = seen.size();
    if(l==4 && last_resatom) {
      int m = inds.size();
      for( int j=0; j<m; j++) {
	out[inds[j]]=true;
      }
    }
  }
   
  
  return out;
}
