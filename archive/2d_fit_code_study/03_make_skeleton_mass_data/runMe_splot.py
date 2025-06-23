from JpsiFitter import JpsiFitter

analysis = JpsiFitter(config_module_name='splot_ctauErr')
# analysis.fitter.loadMassResult('roots/mass.root')
# analysis.fitter.doSplot('/work/pjgwak/pol24/files_roodata/RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250302.root', 'pt>9&&pt<12&&( ((abs(eta1) <= 1.2) && (pt1 >=3.5)) || ((abs(eta2) <= 1.2) && (pt2 >=3.5)) || ((abs(eta1) > 1.2) && (abs(eta1) <= 2.1) && (pt1 >= 5.47-1.89*(abs(eta1)))) || ((abs(eta2) > 1.2) && (abs(eta2) <= 2.1) && (pt2 >= 5.47-1.89*(abs(eta2)))) || ((abs(eta1) > 2.1) && (abs(eta1) <= 2.4) && (pt1 >= 1.5)) || ((abs(eta2) > 2.1) && (abs(eta2) <= 2.4) && (pt2 >= 1.5)) )')
# analysis.fitter.makeSplotPdfs(False, 0.04)
# analysis.fitter.drawSplot('figs/splot_ctauErr.png')
analysis.run_splot_analysis()


print("===== Finish runMe_splot.py =====")