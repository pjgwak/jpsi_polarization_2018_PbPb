# set common fit setting, such as
# input paths
# signal pdf
# bkg pdf
# fixed parameters etc

from pathlib import Path


def get_common_config():
  this_dir = Path(__file__).resolve().parent # directory including the local_config file
  
  files_root = (this_dir / '..' / '..' / 'files_roodata').resolve() # RooDataset folder path
  lib_path = (this_dir / '..' / 'Analysis.so').resolve() # Analysis.so path

  config = {
    # --- input ---
    'mc_input_path': str(files_root / 'RooDataSet_miniAOD_isMC1_PR_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root'),
    'data_input_path': str(files_root / 'RooDataSet_miniAOD_isMC0_Jpsi_cent0_180_Effw1_Accw1_PtW1_TnP1_HFNom_250221.root'),
    'shared_lib_path': str(lib_path),

    # --- output base dir ---
    'output_base_dir': str(this_dir),

    # --- user custom directory label ---
    'date_tag': 'No_Weight_2',

    # --- fit models ---
    # mass fit
    'default_pdf_mass_sig': 'doubleCB',
    'default_pdf_mass_bkg': 'cheby2',

    # ctau fit
    

    # sample information flags
    'default_PR': 2,
    'default_PRw': 1,
    'default_fEffW': False,
    'default_fAccW': False,
    'default_isPtW': False,
    'default_isTnP': False
  }

  return config

# --- check configs ---
# config = get_common_config()

# for key in sorted(config):
#   print(f'{key:<20} : {config[key]}')
