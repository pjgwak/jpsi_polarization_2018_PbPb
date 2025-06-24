# for quick test -> Don't need it later because a main script, runMe.py, will call it
from ROOT import gSystem
gSystem.Load('Analysis.so')

from ROOT import JpsiFitter as JpsiFitterCpp, RooWorkspace, ModelBuilder, RooFit, TCanvas, TFile
import importlib, pprint
from datetime import datetime

class JpsiFitter(object):
  ''' wrapper around cpp class
  '''
  def __init__(self, config_module_name, buffer_path=None):
    # get the fit config
    self.config_name = config_module_name
    print(f"bring a 'configs/{config_module_name}'")
    config_module= importlib.import_module(f'configs.{config_module_name}')
    self.config = config_module.config

    self.buffer_data = None
    if buffer_path:
      print(f'refer {buffer_path}')
      try:
        buffer_module_path = buffer_path.replace('/', '.').replace('.py', '')
        buffer_module = importlib.import_module(buffer_module_path)
        self.buffer_data = buffer_module.snapshot
      except ImportError:
        print(f'there is no {buffer_path}. Continue with default setting')

    # create cpp instances
    print('create cpp instances')
    self.ws = RooWorkspace('wsMy', '')
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
    cuts_config = self.config.get('cuts', {})
    
    if not cuts_config:
      print('Warning: no kinematic cuts are applied')
      selection_cut = io_config.get('selection_cut', '') # no cut -> true
    else:
      kine_vars = cuts_config.get('kinematic', {})
      kine_cut = ("pt>{ptLow} && pt<{ptHigh} && abs(y)>{yLow} && abs(y)<{yHigh} && "
                  "mass>{massLow} && mass<{massHigh} && cBin>={cLow} && cBin<{cHigh}").format(**kine_vars)
      acc_cut = cuts_config.get('acceptance', '1')
      os_cut = cuts_config.get('opposite_sign', '1')

      all_cuts = [cut for cut in [os_cut, acc_cut, kine_cut] if cut !='1']
      selection_cut = ' && '.join(all_cuts)
    
    self.fitter.processTree(
      io_config['input_file'],
      io_config['dataset_name'],
      io_config['observable'],
      selection_cut
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

    # use buffer parameters if possible
    if self.buffer_data:
      final_mc_params = self.buffer_data['final_parameters']

      # free parameter list
      free_params = ['mean_mass', 'sigma1_mass']

      for param_name in sig_params:
        if param_name in final_mc_params:
          mc_val = final_mc_params[param_name]['value']

          if param_name in free_params:
            # use mc results as initial value but it's free
            sig_params[param_name][0] = mc_val
          else:
            sig_params[param_name] = [mc_val]


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

    # is_extended = not self.config['fit_config']['is_mc']
    is_extended = True
    self.fit_result = pdf.fitTo(dataset, RooFit.Save(), RooFit.Range('fit_range'), RooFit.Extended(is_extended))
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

  def _save_to_buffer(self):
    ''' save the fit result as py
    '''
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    buffer_filename = f'buffers/{self.config_name}_{timestamp}.py'

    final_params ={}
    fit_params = self.fit_result.floatParsFinal()
    for i in range(fit_params.getSize()):
      param = fit_params.at(i)
      final_params[param.GetName()] = {
        'value': param.getVal(),
        'error': param.getError()
      }
    
    snapshot_data = {
      'source_config': self.config,
      'fit_status': {
        'status': self.fit_result.status(),
        'covQual': self.fit_result.covQual(),
        'minNll': self.fit_result.minNll()
      },
      'final_parameters': final_params
    }

    # transform a dictionart to string format
    buffer_content_string = pprint.pformat(snapshot_data)

    with open(buffer_filename, 'w', encoding='utf-8') as f:
      f.write(f'# fit snapshot generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n')
      f.write(f'snapshot = {buffer_content_string}')
    
    print(f'fit results saved to {buffer_filename}')
  
  def _save_workspace_to_rootfile(self):
    output_path = self.config.get('output', {}).get('workspace_root_file')

    root_file = TFile(output_path, 'RECREATE')
    if not root_file.IsOpen():
      print(f"Error: Can't open file {output_path} for writing.")
      return
    
    self.ws.Write()
    root_file.Close()

  def run(self):
    print('start analysis')
    self._load_data()
    self._setup_model()
    self._perform_fit()
    self._plot_and_save()
    self._save_to_buffer()
    self._save_workspace_to_rootfile()
    print('analysis finished')

  def build_cut_string(self):
    cuts_config = self.config.get('cuts', {})
    if not cuts_config:
      return "true" # pass all events
    
    # kinematic cuts
    kine_vars = cuts_config.get('kinematic', {})
    kine_cut = ("pt>{ptLow} && pt<{ptHigh} && abs(y)>{yLow} && abs(y)<{yHigh} && "
                "mass>{massLow} && mass<{massHigh} && cBin>={cLow} && cBin<{cHigh}").format(**kine_vars)
                
    acc_cut = cuts_config.get('acceptance', '1')
    os_cut = cuts_config.get('opposite_sign', '1')

    all_cuts = [cut for cut in [os_cut, acc_cut, kine_cut] if cut !='1']
    return ' && '.join(all_cuts)
  
  def run_splot_analysis(self):
    print('===== start SPlot process =====')
    # get mass fit results
    mass_fit_file = self.config.get('mass_fit_rootfile')
    if not mass_fit_file:
      print('Error: mass_fit_rootfile is not found in config')
    self.fitter.loadMassResult(mass_fit_file)

    # get original dataset path and cuts
    splot_inputs = self.config.get('inputs', {})
    input_file = splot_inputs.get('input_file')
    cut_string = self.build_cut_string()
    if not input_file:
      print('Error: inputs:input_file is not found in config')
      return
    self.fitter.doSplot(input_file, cut_string)

    # make PDFs
    hist_config = self.config.get('hist_config', {})
    use_forced_max = hist_config.get('use_forced_range', False)
    forced_max_val = hist_config.get('forced_ctauErrMax', 10)
    self.fitter.makeSplotPdfs(use_forced_max, forced_max_val)

    # draw plot
    output_config = self.config.get('output', {})
    plot_path = output_config.get('plot', 'figs/splot_default.png')
    self.fitter.drawSplot(plot_path)

    print('===== finish SPlot process =====')

