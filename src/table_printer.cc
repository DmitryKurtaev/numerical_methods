#include "include/table_printer.h"

#include <iostream>
#include <iomanip>
#include <algorithm>

void TablePrinter::Print(const std::vector<std::vector<std::string> >& data,
                         int max_n_rows, int max_n_cols, int col_width) {
  bool use_cols_separator = false;

  PrintRows(data, 0, 0, max_n_cols, col_width, use_cols_separator);
  DrawRow('-', data[0].size(), col_width, max_n_cols);
  

  if (data.size() > max_n_rows) {
    PrintRows(data, 1, max_n_rows / 2 - 1, max_n_cols, col_width,
              use_cols_separator);
    DrawRow('.', data[max_n_rows / 2 - 1].size(), col_width, max_n_cols);
    PrintRows(data, data.size() - max_n_rows / 2, data.size() - 1,
              max_n_cols, col_width, use_cols_separator);
  } else {
    PrintRows(data, 1, data.size() - 1, max_n_cols, col_width,
              use_cols_separator);
  }
}

void TablePrinter::PrintRows(const std::vector<std::vector<std::string> >& data,
                             int from, int to, int max_n_cols, int col_width,
                             bool& use_cols_separator) {
  std::cout.fill(' ');
  for (int i = from; i <= to; ++i) {
    int size = data[i].size();

    use_cols_separator = use_cols_separator || size > max_n_cols;

    std::cout << "|";
    if (size > max_n_cols) {
      for (int j = 0; j < max_n_cols / 2; ++j) {
        std::cout << std::setw(col_width) << data[i][j] << "|";
      }
      if (use_cols_separator) {
        std::cout.fill('.');
        std::cout << std::setw(4) << "|";
        std::cout.fill(' ');
      }
      for (int j = size - max_n_cols / 2; j < size; ++j) {
        std::cout << std::setw(col_width) << data[i][j] << "|";
      }
    } else {
      for (int j = 0; j < size; ++j) {
        std::cout << std::setw(col_width) << data[i][j] << "|";
      }
    }
    std::cout << std::endl; 
  }
}

void TablePrinter::DrawRow(char filler, int n_cols, int col_width,
                           int max_n_cols) {
  std::cout << "|";
  std::cout.fill(filler);

  if (n_cols > max_n_cols) {
    for (int i = 0; i < std::min(max_n_cols / 2, n_cols - 1); ++i) {
      std::cout << std::setw(col_width + 1) << "|";
    }
    if (n_cols > max_n_cols) {
      std::cout.fill('.');
      std::cout << std::setw(4) << "|";
      std::cout.fill(filler);
    }
    for (int i = n_cols - max_n_cols / 2; i < n_cols; ++i) {
      std::cout << std::setw(col_width + 1) << "|";
    }
  } else {
    for (int i = 0; i < n_cols; ++i) {
      std::cout << std::setw(col_width + 1) << "|";
    }
  }


  std::cout << std::endl;
}