//=================================
// include guard
#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__

//=================================
// included dependencies
#include <string>
using namespace std;

//=================================
// the actual functions
// trim from start
string ltrim(string s);

// trim from end
string rtrim(string s);

// trim from both ends
string trim(string s);

// function get hexadecimal
int getHex(string hexstr);

// split a string
vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);


#endif
