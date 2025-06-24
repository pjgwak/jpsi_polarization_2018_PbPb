from JpsiFitter import JpsiFitter
import os, glob

# find latest MC buffer file
buffer_dir = 'buffers'
mc_fit_results = glob.glob(os.path.join(buffer_dir, 'mc_mass_*.py'))

if not mc_fit_results:
    print(f"Error: Can't find mc results")
    print('Please fitthe mc mass first')
else:
    latest_buffer_file = max(mc_fit_results, key=os.path.getctime)
    print(f'bring: {latest_buffer_file}')

    analysis = JpsiFitter(config_module_name='data_mass', buffer_path=latest_buffer_file)
    analysis.run()
    