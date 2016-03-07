#ifndef INCLUDE_TABLE_PRINTER_H_
#define INCLUDE_TABLE_PRINTER_H_

#include <string>
#include <vector>

class TablePrinter {
 public:
  static void Print(const std::vector<std::vector<std::string> >& data,
                    int max_n_rows = 25, int max_n_cols = 6,
                    int col_width = 13);

 private:
  static void PrintRows(const std::vector<std::vector<std::string> >& data,
                        int from, int to, int max_n_cols, int col_width,
                        bool& use_cols_separator);

  static void DrawRow(char filler, int n_cols, int col_width, int max_n_cols);
};

#endif  // INCLUDE_TABLE_PRINTER_H_