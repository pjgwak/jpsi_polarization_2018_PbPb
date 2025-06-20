#include "Algorithm.h"

Algorithm::Algorithm() {
  // skip
}

void Algorithm::initialize(Data* data) {
  m_data = data;
  
}

void Algorithm::execute() {
  // fill the histogram
  if (h_runNb_m)
    h_runNb_m->Fill(m_data->runNb);
}

