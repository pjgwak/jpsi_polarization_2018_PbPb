#pragma once

#include "RooWorkspace.h"
#include <string>
#include <vector>
#include <map>

class ModelBuilder {
public:
  ModelBuilder(RooWorkspace &ws);
  
  // set fitting models
  void buildMassSignal();
  void buildMassBkg();
  void buildMassModel();
  // void combineModels();

private:
  RooWorkspace &m_ws;
};