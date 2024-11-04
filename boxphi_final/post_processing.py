import numpy as np
import sys
import os
import shutil
from cobaya.yaml import yaml_load_file, yaml_dump_file

fileroot = sys.argv[1]
filename = sys.argv[2]

idx = 1
while True:
  ipath = os.path.join(fileroot, filename)+'.%d.txt'%idx
  print('Looking for %s'%ipath)
  if not os.path.isfile(ipath):
    print('File %s does not exist. Stop!'%ipath)
    break
  opath = os.path.join(fileroot, 'final/%s'%filename)+'.%d.txt'%idx
  tmp = open(ipath)
  old_names = tmp.readline().replace('#','').split()
  tmp.close()

  old_chain = np.loadtxt(ipath)
  if len(old_names) != old_chain.shape[1]:
    print('name not match, sth wrong')
    print(old_names)
    print(old_chain[0])
    break
  new_chain = np.empty([old_chain.shape[0], old_chain.shape[1]+2])
  new_names = old_names.copy()
  new_names.insert(old_names.index('chi2__CMB')+1, 'chi2__CMBnL')
  new_names.insert(old_names.index('chi2__CMB')+2, 'S8')
  for par in old_names:
    new_chain[:, new_names.index(par)] = old_chain[:, old_names.index(par)]
  new_chain[:, new_names.index('chi2__CMBnL')] = new_chain[:, new_names.index('chi2__planck_2020_hillipop.TTTEEE')] + new_chain[:, new_names.index('chi2__planck_2018_lowl.TT')] + new_chain[:, new_names.index('chi2__planck_2020_lollipop.lowlE')]
  new_chain[:, new_names.index('chi2__CMB')] = new_chain[:, new_names.index('chi2__CMBnL')] + new_chain[:, new_names.index('chi2__planckpr4lensing')]
  new_chain[:, new_names.index('S8')] = new_chain[:, new_names.index('sigma8')]*np.sqrt(new_chain[:, new_names.index('omegam')]/0.3)
  
  with open(opath, 'w') as o:
    t = '#\t'
    t += '\t'.join(new_names)
    t += '\n'
    o.write(t)
    for line in new_chain:
      dat = ''
      for i in range(len(line)):
        if new_names[i] == 'weight':
          dat += '\t%d'%line[i]
        else:
          dat += '\t%.8e'%line[i]
      dat += '\n'
      o.write(dat)  
  
  print('File %s done and outputed to %s'%(ipath, opath))
  idx += 1

shutil.copyfile(os.path.join(fileroot, filename)+'.covmat', os.path.join(fileroot, 'final/%s'%filename)+'.covmat')
old_yaml = yaml_load_file(os.path.join(fileroot, filename)+'.updated.yaml')
new_yaml = old_yaml.copy()
new_yaml['params']['chi2__CMBnL'] = {'derived': True, 'latex': '\\chi^2_{CMBnL}'}
new_yaml['params']['S8'] = {'value': 'lambda sigma8,omegam: sigma8*np.sqrt(omegam/0.3)', 'derived': True, 'latex': 'S_8'}
yaml_dump_file(os.path.join(fileroot, 'final/%s'%filename)+'.updated.yaml', new_yaml,error_if_exists=False)
  