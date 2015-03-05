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

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    if(!item.empty())
      elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}
