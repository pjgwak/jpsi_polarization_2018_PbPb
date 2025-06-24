config = {
  # dictionary
  'io': { # input, output setting
    'input_file': '/work/pjgwak/pol24/files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250302.root',
    'dataset_name': 'ds_mass', # name of dataset inside input file
    'observable': 'mass',
    'output_plot': 'figs/mass.png'
  },

  'fit_config': {
    'is_mc': False,
    'fit_range': [2.6, 3.5], # mass fit range
    'sig_pdf': 'doubleCb',
    'bkg_pdf': 'cheby2'
  },

  # sig parameters
  'mass_sig_params': {
    'doubleCb': {
      "mean_mass":   [3.096, 3.08, 3.11],
      "sigma1_mass": [0.02, 0.001, 0.1],
      "alpha1_mass": [1.2, 0.5, 5.0], # dummy. will be fixed
      "n1_mass":     [1.5, 0.5, 5.0], # dummy
      "x_mass":      [1.5, 0.1, 5.0], # dummy
      "frac_mass":   [0.5, 0.0, 1.0]  # dummy
    }
  },

  # bkg parameters
  'mass_bkg_params': {
    'cheby2': {
      "s1_mass": [0.01, -1, 1],
      "s2_mass": [0.01, -1, 1],
    },
    'cheby3': {
      "s1_mass": [0.01, -1, 1],
      "s2_mass": [0.01, -1, 1],
      "s3_mass": [0.05, -1, 1] 
    }
  },

  # sig, bkg particles numbers
  'yields': {
    'nSig': [10000, 5000, 100000],
    'nBkg': [50000, 30000, 1000000],
  },

  # basic cuts
  'cuts': {
      'kinematic': {
          'ptLow': 9, 'ptHigh': 12.0,
          'yLow': 0.0, 'yHigh': 1.6,
          'massLow': 2.6, 'massHigh': 3.5,
          'cLow': 0, 'cHigh': 180
      },
      # fixed cuts - dont' touch
      'acceptance': "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2) && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )",
      'opposite_sign': "recoQQsign==0"
  },
  
  'output': {
      'plot': 'figs/mass.png',
      'workspace_root_file': 'roots/mass.root'
  }
}