#pragma once

#include <string>
#include "RooWorkspace.h"
#include "memory" // std::unique_ptr

class JpsiFitter {
public:
  JpsiFitter(RooWorkspace &ws);
  ~JpsiFitter();

  void processTree(const std::string &filePath, const std::string &treeName, const std::string &obs, const std::string &cut);
  void print(); // show the contents of workspace

private:
  RooWorkspace &m_ws;
};