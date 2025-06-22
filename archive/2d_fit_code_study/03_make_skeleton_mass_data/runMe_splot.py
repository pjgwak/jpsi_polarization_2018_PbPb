from JpsiFitter import JpsiFitter

analysis = JpsiFitter(config_module_name='splot_ctauErr')
analysis.fitter.loadMassResult()
analysis.fitter.doSplot()
analysis.fitter.makeSplotPdfs()
analysis.fitter.drawSplot()


print("===== Finish runMe_splot.py =====")