#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

class BadConversion : public std::runtime_error {
 public:
 BadConversion(const std::string& s)
   : std::runtime_error(s)
    { }
};

inline double stringToDouble(const std::string& s)
{
  if(s=="") {
    return(NA_REAL);
  }
    
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    throw BadConversion("stringToDouble(\"" + s + "\")");
  return x;
}

inline int stringToInt(const std::string& s)
{
  if(s=="") {
    return(NA_INTEGER);
  }
  
  std::istringstream i(s);
  int x;
  if (!(i >> x))
    throw BadConversion("stringToInt(\"" + s + "\")");
  return x;
}
