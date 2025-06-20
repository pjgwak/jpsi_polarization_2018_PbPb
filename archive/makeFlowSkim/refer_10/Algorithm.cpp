#include "Algorithm.h"

Algorithm::Algorithm() {
  // skip
}

void Algorithm::initialize(std::shared_ptr<Data> data)
{
  m_data = data;
}

void Algorithm::execute() {
  // apply selection
  if (!passedSelection())
    return;

  // calculate the event weight
  calculateWeight();

  // fill the plots
  fillPlots();
}

bool Algorithm::passedSelection() {
  // select specific flavor cominations
  if (cut_selectSF && m_data->tau_0 != m_data->tau_1)
    return false;
  if (cut_selectDF && m_data->tau_0 == m_data->tau_1)
    return false;
  
  //  select invariant mass of the tau letpon visible decay products
  if(cut_maxVisMass.enabled() && cut_maxVisMass < (*(m_data->tau_0_p4) + *(m_data->tau_1_p4)).M())
    return false;

  // passed all cuts
  return true;
}

void Algorithm::calculateWeight() {
  // apply the sample weight
  m_eventWeight = p_sampleWeight;

  // apply the MC weight, if requested
  if (p_isMC)
    m_eventWeight *= m_data->weight_mc;
}
void Algorithm::fillPlots(){
  // fill the histograms
  if (h_ditau_m)
    h_ditau_m->Fill(m_data->ditau_mmc_mlm_m, m_eventWeight);
  if (h_ditau_visM)
    h_ditau_visM->Fill((*(m_data->tau_0_p4) + *(m_data->tau_1_p4)).M(), m_eventWeight);
}