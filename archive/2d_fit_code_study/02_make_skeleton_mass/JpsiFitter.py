# for quick test -> Don't need it later because a main script, runMe.py, will call it
from ROOT import gSystem
gSystem.Load('Analysis.so')

from ROOT import JpsiFitter as JpsiFitterCpp, RooWorkspace, ModelBuilder, RooFit, TCanvas
import importlib

class JpsiFitter(object):
  ''' wrapper around cpp class
  '''
  def __init__(self, config_module_name='mc_mass'):
    # get the fit config
    print(f"bring a 'configs/{config_module_name}'")
    config_moudle = importlib.import_module(f'configs.{config_module_name}')
    self.config = config_moudle.config

    # create cpp instances
    print('create cpp instacnes')
    self.ws = RooWorkspace('ws', '')
    self.model = ModelBuilder(self.ws)
    self.fitter = JpsiFitterCpp(self.ws)

    # fit results
    self.fit_result = None
    print('Finish initiation: JpsiFitter')
  
  def _load_data(self):
    '''
    read config's 'io' and call the data
    '''
    print('load input sample')
    io_config = self.config['io']
    self.fitter.processTree(
      io_config['input_file'],
      io_config['dataset_name'],
      io_config['observable'],
      io_config['selection_cut']
    )
  
  def _setup_model(self):
    '''
    build fit models
    '''
    print('set fit models up ')
    fit_config = self.config['fit_config']
    sig_pdf_type = fit_config['sig_pdf']
    bkg_pdf_type = fit_config['bkg_pdf']
    is_mc = fit_config['is_mc']

    # pdf params
    sig_params = self.config['mass_sig_params'].get(sig_pdf_type, {})
    bkg_params = self.config['mass_bkg_params'].get(bkg_pdf_type, {})
    yield_params = self.config.get('yields', {})
    all_params = {**sig_params, **bkg_params, **yield_params} # combine sig and bkg params

    # build pdfs in cpp
    self.model.setParameters(all_params)
    self.model.buildMassSig(sig_pdf_type)
    if not is_mc:
      self.model.buildMassBkg(bkg_pdf_type)
    self.model.buildMassModel(is_mc)
  
  def _perform_fit(self):
    print('start a fit')
    pdf = self.ws.pdf('massModel')
    dataset = self.ws.data(self.config['io']['dataset_name'])
    mass = self.ws.var(self.config['io']['observable'])

    fit_range = self.config['fit_config']['fit_range']
    mass.setRange('fit_range', fit_range[0], fit_range[1])

    self.fit_result = pdf.fitTo(dataset, RooFit.Save(), RooFit.Range('fit_range'), RooFit.Extended())
    self.fit_result.Print()
  
  def _plot_and_save(self):
    canvas = TCanvas('c', '', 800, 600)
    mass = self.ws.var(self.config['io']['observable'])
    dataset = self.ws.data(self.config['io']['dataset_name'])
    pdf = self.ws.pdf('massModel')

    frame = mass.frame(RooFit.Title("Fit Result"), RooFit.Range("fit_range"))

    dataset.plotOn(frame)
    pdf.plotOn(frame, RooFit.NormRange("fit_range"), RooFit.Range("fit_range"))

    # pdf.plotOn(frame, RooFit.Components("mass_sig"), RooFit.LineStyle(2), RooFit.LineColor(2))
    # pdf.plotOn(frame, RooFit.Components("mass_bkg"), RooFit.LineStyle(3), RooFit.LineColor(4))

    frame.Draw()
    output_file = self.config['io']['output_plot']
    canvas.SaveAs(output_file)

  
  def run(self):
    print('start analysis')
    self._load_data()
    self._setup_model()
    self._perform_fit()
    self._plot_and_save()
    print('analysis finished')


print('\n===== Test announcement: Finish JpsiFitter.py =====')