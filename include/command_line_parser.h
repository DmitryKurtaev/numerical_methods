// Copyright 2015 Dmitry Kurtaev

#ifndef INCLUDE_COMMAND_LINE_PARSER_H_
#define INCLUDE_COMMAND_LINE_PARSER_H_

#include <stdlib.h>
#include <string>
#include <typeinfo>
#include <iostream>

class CommandLineParser {
 public:
  CommandLineParser(int argc, char** argv);

  template<typename T>
  T Get(const std::string& key);

  bool Exists(const std::string& key);

 private:
  int argc_;
  char** argv_;
};

template<typename T>
T CommandLineParser::Get(const std::string& key) {
  int idx = -1;
  for (int i = 1; i < argc_; ++i) {
    if (argv_[i] == "-" + key ||
        argv_[i] == "--" + key) {
      idx = i;
      break;
    }
  }

  if (idx != -1) {
    int int_value;
    float float_value;
    double double_value;
    if (typeid(T) == typeid(int_value))
      return atoi(argv_[idx + 1]);
    if (typeid(T) == typeid(float_value) ||
        typeid(T) == typeid(double_value))
      return atof(argv_[idx + 1]);
  }
  return 0;
}

#endif  // INCLUDE_COMMAND_LINE_PARSER_H_
