config = {
    # later: use published mass config according to kinematic bin - but it's a test
    'mass_fit_rootfile': 'roots/mass.root', # use mass fit results

    # input, output
    'io': {
        'input_file': '/work/pjgwak/pol24/files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250302.root',
        'dataset_name': 'splot_dataset',
        'observable': 'ctau3DErr'
    },

    # cuts
    'cuts': {
        'kinematic': {
            'ptLow': 3.0, 'ptHigh': 6.5,
            'yLow': 1.6, 'yHigh': 2.4,
            'massLow': 2.6, 'massHigh': 3.5,
            'cLow': 60, 'cHigh': 180
        },
        'angle_cut': {
            'cos_low': 0.0, 'cos_high': 0.1
        },
        'acceptance': "( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2) && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )",
        'opposite_sign': "recoQQsign==0"
    },

    # histogram setting
    'hist_config': {
        'initial_range': [0.0, 0.25],
        'initial_bins': 100,
        'auto_range_finding': True
    },

    # output path
    'output': {
        'plot': 'figs/splot_ctauErr.png',
        'root_file': 'roots/splot_ctauErr.root'
    }
}