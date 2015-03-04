#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include "utils.h"
using namespace std;

// trim from start
string ltrim(string s) {
  s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
  return s;
}

// trim from end
string rtrim(string s) {
  s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
  return s;
}

// trim from both ends
string trim(string s) {
  return ltrim(rtrim(s));
}

// function get hexadecimal
int getHex(string hexstr) {
  return (int)strtol(hexstr.c_str(), 0, 16);
}
