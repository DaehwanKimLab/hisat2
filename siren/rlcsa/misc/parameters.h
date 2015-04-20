#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstdio>
#include <fstream>
#include <map>

#include "definitions.h"


namespace CSA
{


typedef std::pair<std::string, usint> parameter_type;


class Parameters
{
  public:
    Parameters() {};
    ~Parameters() {};

    bool contains(const std::string& key) const;
    usint get(const std::string& key) const;
    usint get(const parameter_type& param) const;
    void set(const std::string& key, usint value);
    void set(const parameter_type& param);

    void read(std::ifstream& file);
    void read(FILE* file);
    void read(const std::string& file_name);
    void print() const;
    void write(std::ostream& stream) const;
    void write(FILE* file) const;
    void write(const std::string& file_name) const;

  private:
    std::map<std::string, usint> parameters;
};


} // namespace CSA


#endif
