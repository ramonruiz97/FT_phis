# Flavour Combination of OS taggers
__author__ = ["Ramón Ángel Ruiz Fernández"],


import numpy as np
import ipanema
import complot
from uncertainties import unumpy as unp
import argparse
import yaml
import os
import uproot3 as uproot




def taggerCombination(qs, etas):
    pB = np.ones_like(qs[0])
    pBbar = np.ones_like(qs[0])
    
    #Events w/ omega > 0.5 < 0. do not contribute 
    #Warning Sevda try ! -> Comment following lines for Sevda try
    for t in range(len(etas)):
      condition = (etas[t]>0.5) | (etas[t]<0.)
      etas[t] = np.where(condition, 0.5, etas[t])
      qs[t] = np.where(condition, 0., qs[t])


    for i in range(len(qs)):
        pB    *= 0.5*(1.+qs[i]) - qs[i]*(1.-etas[i])
        pBbar *= 0.5*(1.-qs[i]) + qs[i]*(1.-etas[i])

    q = np.zeros_like(qs[0])
    eta = 0.5*np.ones_like(qs[0])

    PB = pB/(pB+pBbar)
    PBbar = 1. - PB

    q = np.where(PB > PBbar, -1., q)
    eta = np.where(PB > PBbar, 1.-PB, eta)
    q = np.where(PBbar > PB, 1., q)
    eta = np.where(PBbar > PB, 1.-PBbar, eta)
    return q, eta




if __name__ == '__main__':
  parser = argparse.ArgumentParser(
        description='Tagging Calibration for B+ and MC')
  parser.add_argument('--data', 
                      help='Input data.')

  parser.add_argument('--weight', 
                      help='Weight to apply to data')

  parser.add_argument('--mode', 
                      help='Mode')

  parser.add_argument('--tagger', 
                      help='OSCombination', default="OSCombination")

  parser.add_argument('--version', 
                      help='Version of tuples')

  parser.add_argument('--model', default='linear')

  parser.add_argument('--calibrations', default=False)

  parser.add_argument('--output-sample', 
						help='root file with Combined omega')

  args = vars(parser.parse_args())

  tp = args["data"]
  data = ipanema.Sample.from_root(tp)
  tagger = args["tagger"]
  mode = args["mode"]
  order = args["model"]

  if args['calibrations']:
    calibrations = args["calibrations"].split(",")
  else:
    calibrations = False
  
  #We could take from yaml w/ taggers information
  q_list = ['B_OSKaonLatest_TAGDEC',
						'B_OSMuonLatest_TAGDEC',
						'B_OSElectronLatest_TAGDEC',
						'B_OSCharm_TAGDEC',
						'B_OSVtxCh_TAGDEC']

  eta_list = ['B_OSKaonLatest_TAGETA',
							'B_OSMuonLatest_TAGETA',
							'B_OSElectronLatest_TAGETA',
							'B_OSCharm_TAGETA',
							'B_OSVtxCh_TAGETA']

  q_list = [np.float64(data.df[b]) for b in q_list]
  eta_list = [np.float64(data.df[b]) for b in eta_list]


  if calibrations:
    list_pars = [ipanema.Parameters.load(p) for p in calibrations]
    if order=="linear":
      cal = [[p['p0'].value, p['p1'].value, p['eta_bar'].value] for p in list_pars]
      eta_cal = [np.float64(cal[i][0] + cal[i][1] * (eta_list[i]-cal[i][2])) for i in range(len(eta_list))]
    else:
      cal = [[p['p0'].value, p['p1'].value, p['p2'].value, p['eta_bar'].value] for p in list_pars]
      eta_cal = [np.float64(cal[i][0] + cal[i][1] * (eta_list[i]-cal[i][3])+ cal[i][2] * (eta_list[i]-cal[i][3])**2) for i in range(len(eta_list))]
    q, eta = taggerCombination(q_list, eta_cal)
  else:
    q, eta = taggerCombination(q_list, eta_list)

  data.df["OSCombination_TAGDEC"] = q
  data.df["OSCombination_TAGETA"] = eta

  with uproot.recreate(args['output_sample']) as f:
    _branches = {}
    for k, v in data.df.items():
      if 'int' in v.dtype.name:
        _v = np.int32
      elif 'bool' in v.dtype.name:
        _v = np.int32
      else:
        _v = np.float64
      _branches[k] = _v
    mylist = list(dict.fromkeys(_branches.values()))
    f["DecayTree"] = uproot.newtree(_branches)
    f["DecayTree"].extend(data.df.to_dict(orient='list'))

	










