# linear_fit
#
#

__all__ = []
__author__ = ["name"]
__email__ = ["email"]

import os
import uproot3 as uproot
import ipanema
import numpy as np
import complot
import matplotlib.pyplot as plt
import argparse

ipanema.initialize('cuda', 1) #To change to also use clusterA -> env variable

def linear_fit(sigmat, sigmaeff, sigmaeff_err):
  def model(x, *args): #Model
	  return args[0] + args[1]*x

  def fcn(params, x, y): #Func to minimize by Minuit
    p = list(params.valuesdict().values()) #p0 and p1
    return (y - model(x, *p))**2/(sigmaeff_err**2)

  pars = ipanema.Parameters()
  pars.add(dict(name='p0', value=0, free=True))
  pars.add(dict(name='p1', value=1., free=True))
  res = ipanema.optimize(fcn, pars, fcn_args=(sigmat, sigmaeff),
												 method='minos', verbose=True)
  return res.params

def calibration(df, params_bin, sigma_range, 
								time_range, nBin, plots=False):

  tLL, tUL = time_range
  sLL, sUL = sigma_range[0], sigma_range[nBin-1]
	
  sigmat = []
  sigmat_errh = [] #For plotting
  sigmat_errl = [] #For plotting
  sigma_eff = []
  sigmaeff_err = []
  _n = [] # To norm
  
  for i, sLL, sUL in zip(range(nBin), sigma_range[:-1], sigma_range[1:]): #TODO: limit of Nbin
    list_of_cuts = [
            f"time>{tLL} & time<{tUL}",
            f"sigmat>{sLL} & sigmat<{sUL}"
        ]
    cut = "(" + ") & (".join(list_of_cuts) + ")"
    cdf = df.query(cut)
    sigmat.append(params_bin[i]['sigmaAverage'].value)
    sigma_eff.append(params_bin[i][f'sigmaeff'].value)
    sigmaeff_err.append(params_bin[i][f'sigmaeff'].stdev)
  	#For plotting:
    sigmat_errl.append(sigmat[-1] - sLL)
    sigmat_errh.append(sUL - sigmat[-1])
    _n.append(cdf.shape[0])

  sigmat = np.array(sigmat)
  sigmat_errl = np.array(sigmat_errl)
  sigmat_errh = np.array(sigmat_errh)
  sigma_eff = np.array(sigma_eff)
  sigmaeff_err = np.array(sigmaeff_err)

  print(sigmat)
  print(sigmat_errl)
  print(sigmat_errh)
  print(sigma_eff)
  print(sigmaeff_err)


  pars = linear_fit(sigmat, sigma_eff, sigmaeff_err)
  if plots:
    def model(x, p):
      return p[0] + p[1]*(x)
    sigma_proxy = np.linspace(np.min(sigma_range), np.max(sigma_range), 100)
    fig, axplot = complot.axes_plot()
    axplot.plot(sigma_proxy, model(sigma_proxy, pars.valuesarray()))
    axplot.errorbar(sigmat, sigma_eff, yerr=sigmaeff_err,
                    xerr=[sigmat_errl, sigmat_errh],
                    fmt='.', color=f'k')
    result = "\n".join((f"p0 = {np.round(pars['p0'].value,5)} +- {np.round(pars['p0'].stdev,5)}",
		                f"p1 = {np.round(pars['p1'].value,5)} +- {np.round(pars['p1'].stdev,5)}"))
    axplot.text(0.05, 0.95, result, transform=axplot.transAxes, fontsize=12,
		            verticalalignment='top')

    axplot.set_xlabel(r'$\sigma_t$ [ps]')
    axplot.set_ylabel(r'$\sigma_{eff}$ [ps]')
    fig.savefig(os.path.join(plots, "linear.pdf"))

  return pars

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
        description='Calibration of the decay time resolution.')
  parser.add_argument('--data', help='Input prompt data.')
  parser.add_argument('--json-bin', help='Calibrations of each bin')
  parser.add_argument('--maxbin', default=10,
											help= 'maximum value of the bin taken into account in the calibration')
  parser.add_argument('--output-json', help='Result of fit parameters')
  parser.add_argument('--output-plots', help='Linear plot')


  args = vars(parser.parse_args())
  branches = ["time", "sigmat"]
  df = uproot.open(args['data'])["DecayTree"].pandas.df(branches=branches)
  # time_range = [-3., 7.]
  nBin = int(args["maxbin"])
  timeres_binning = [0., 0.023, 0.027, 0.030, 0.033, 0.036, 0.038, 0.042, 0.047, 0.062, 0.090]
  
  p = [ipanema.Parameters.load(p) for p in args['json_bin'].split(',')]
  time_range = [p[0]["tLL"].value, p[0]["tUL"].value]
  print(time_range)

  os.makedirs(args['output_plots'], exist_ok=True)

  pars = calibration(df, p, timeres_binning, time_range, 
										 nBin=nBin,  plots=args['output_plots'])

  pars.dump(args["output_json"])





# vim: fdm=marker ts=2 sw=2 sts=2 sr noet
