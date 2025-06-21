config = {
  # dictionary
  'io': { # input, output setting
    'input_file': '/work/pjgwak/pol24/files_roodata/RooDataSet_miniAOD_isMC1_NP_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root',
    'dataset_name': 'dataset', # name of dataset inside input file
    'observable': 'mass',
    'selection_cut': 'mass > 2.6 && mass < 3.5', # for later
    'output_plot': 'fit_plot.png'
  },

  'fit_config': {
    'is_mc': True,
    'fit_range': [2.6, 3.5], # mass fit range
    'sig_pdf': 'doubleCb',
    'bkg_pdf': 'cheby3'
  },

  # sig parameters
  'mass_sig_params': {
    'doubleCb': {
      "mean_mass":   [3.096, 3.08, 3.11],
      "sigma1_mass": [0.04, 0.01, 0.1],
      "alpha1_mass": [1.2, 0.5, 5.0],
      "n1_mass":     [1.5, 0.5, 5.0],
      "x_mass":      [1.5, 0.1, 5.0], # sigma2 = sigma1 * x
      "frac_mass":   [0.5, 0.0, 1.0]
    }
  },

  # bkg parameters
  'mass_bkg_params': {
    'cheby3': {
      "s1_mass": [0.1, -1, 1],
      "s2_mass": [-0.1, -1, 1],
      "s3_mass": [0.05, -1, 1] 
    }
  },

  # sig, bkg particles numbers
  'yields': {
    'nSig': [1000, 1, 50000],
    'nBkg': [50000, 1, 200000],
    # for MC
    'nSig_MC': [100000, 1, 8000000]
  }
}