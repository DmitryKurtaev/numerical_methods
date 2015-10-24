// Copyright 2015 Dmitry Kurtaev

#include "include/command_line_parser.h"

CommandLineParser::CommandLineParser(int argc, char** argv)
  : argc_(argc),
    argv_(argv) {
}

bool CommandLineParser::Exists(const std::string& key) {
  for (int i = 0; i < argc_; ++i) {
    if (argv_[i] == "-" + key ||
        argv_[i] == "--" + key) {
      return true;
    }
  }
  return false;
}


