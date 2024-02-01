__all__ = []
__author__ = ["Ramón Ángel Ruiz Fernández"]
__email__ = ["rruizfer@cern.ch"]


import uproot3 as uproot
import os
# import sys
# import json
# import math as m
import argparse
# import math
import uncertainties
import ipanema
import numpy as np
import complot
import matplotlib.pyplot as plt


ipanema.initialize('cuda', 1)
prog = ipanema.compile(open('scripts/timeres.c').read())



def model1(time, fres, mu, sigma, fprompt, fll, fsl, fwpv, taul, taus, tau1,
           tau2, share, tLL, tUL, prob=False, norm=1):
  """"
  Model1: 2 resolution gaussians + long-live component + wrong PV
  complete me
  """
  prog.kernel_time_fit(prob, time, fres,
                       np.float64(mu), sigma, np.float64(fprompt),
                       np.float64(fll), np.float64(fsl), np.float64(fwpv),
                       np.float64(taul), np.float64(taus),
                       np.float64(tau1), np.float64(tau2), np.float64(share),
                       np.float64(tLL), np.float64(tUL),
                       global_size=(len(time),))
  return norm * prob


def fcn(params, time, prob, weight, norm=1):
  p = params.valuesdict()
  pbin = ''
  fres = ipanema.ristra.allocate(np.float64([1 - p[f'f{pbin}'], p[f'f{pbin}']])) #2 Gaussians
  # fres = ipanema.ristra.allocate(np.float64([p[f'f1'], p[f'f2'], p[f'f3']])) #3 Gauss
  mu = p[f'mu{pbin}']
  sigma = ipanema.ristra.allocate(np.float64(           # 2 Gaussians
      [p[f'sigma1{pbin}'], p[f'sigma2{pbin}']]))
  # sigma = ipanema.ristra.allocate(np.float64(
      # [p[f'sigma1'], p[f'sigma2'], p[f'sigma3']] ) ) # 3 Gaussians
  fprompt = p[f'fprompt{pbin}']
  fll = p[f'fll{pbin}']
  fsl = p[f'fsl{pbin}']
  fwpv = p[f'fwpv{pbin}']
  taul = p[f'taul{pbin}']
  taus = p[f'taus{pbin}']
  tau1 = p[f'tau1{pbin}']
  tau2 = p[f'tau2{pbin}']
  share = p[f'share{pbin}']
  tLL = p[f'tLL']
  tUL = p[f'tUL']
  model1(time=time, prob=prob, fres=fres, mu=mu, sigma=sigma,
         fprompt=fprompt, fll=fll, fsl=fsl, fwpv=fwpv, taul=taul, taus=taus,
         tau1=tau1, tau2=tau2, share=share, tLL=tLL, tUL=tUL, norm=norm)
  num = ipanema.ristra.get(prob)
  den = 1.
  return -2.0 * np.log(num / den)*ipanema.ristra.get(weight)


def time_resolution(df, time_range, sigma_range, wpv,  mode,
                    weight=False, cut=False, figs=False):
  """
  Master function that creates the fit to the binned in the time error bins
  sample.
  TODO: type_simult functionality
  <> -< >-
  """
  if weight:
    print("time resolution: using the weight ", weight)

  # diego suggests not to use this if

  tLL, tUL = time_range
  sLL, sUL = sigma_range
  time = 'time'
  sigmat = 'sigmat'

  list_of_cuts = [
      f"{time}>{tLL} & {time}<{tUL}",
      f"{sigmat}>{sLL} & {sigmat}<{sUL}"
  ]

  cut = "(" + ") & (".join(list_of_cuts) + ")"
  print(cut)
  cdf = df.query(cut)
  weight = False
  # weight = "sig_sw"

  DM = 17.74

  pars = ipanema.Parameters()
  #Constants
  pars.add(dict(name="DM", value=DM, free=False))
  pars.add(dict(name="tLL", value=tLL, free=False))
  pars.add(dict(name="tUL", value=tUL, free=False))

  #Gauss
  pars.add(dict(name="fprompt", value=0.99,
                min=0.1, max=1, free=True, 
                latex=r"f_{\mathrm{prompt}}"))
  pars.add(dict(name="f", value=0.73,
                min=0, max=1, free=True)) #2 Gaussians


  # pars.add(dict(name="f1", value=0.73,
                # min=0, max=1, free=True)) #3 Gaussians
  # pars.add(dict(name="fproxy2", value=0.1,
                # min=0, max=1, free=True)) #3 Gaussians
  # pars.add(dict(name="f2", formula="(1-f1)*fproxy2")) #3 Gaussians
  # pars.add(dict(name="f3", formula="(1-f1)*(1-fproxy2)")) #3 Gaussians

  pars.add(dict(name="mu", value=0,
                min=tLL, max=tUL, free=True))

  # pars.add(dict(name="sigmap", value=0.025, #2 Gaussians
  #               min=0.0001, max=20, free=True))
  # pars.add(dict(name="sigmapp", value=0.010,
  #               min=0.0001, max=20, free=True))
  # pars.add(dict(name="sigma1",
  #               formula="sigmap - sigmapp * sqrt((f)/(1-f))"))
  # pars.add(dict(name="sigma2",
  #               formula="sigmap + sigmapp * sqrt((1-f)/(f))"))
  pars.add(dict(name="sigma1", value=0.05, min=0.001, max=1., free=True))
  pars.add(dict(name="sigma2", value=0.10, min=0.001, max=1., free=True))

  # pars.add(dict(name="sigmap", value=0.025, #2 Gaussians
                # min=0.001, max=20))
  # pars.add(dict(name="sigmapp", value=0.010,
                # min=0.001, max=20, free=True))
  # pars.add(dict(name="sigma1",
                # formula="sigmap - sigmapp * sqrt((f2)/(1-f2))"))
  # pars.add(dict(name="sigma2",
                # formula="sigmap + sigmapp * sqrt((1-f2)/(f2))"))
  # pars.add(dict(name="sigma1", value=0.25, min=0.01, max=0.6, free=True))
  # pars.add(dict(name="sigma2", value=0.10, min=0.01, max=1., free=True))
  # pars.add(dict(name="sigma3", value=.1, min=.01, max=0.7, free=True))

  pars.add(dict(name="fphys2",
                value = 0.9, #For the moment no more components to dspi
                min = 0.7,
                max = 1.0,
                free = True))
  # pars.add(dict(name="fphys2",
                # value=0.99,
                # free = False))

  pars.add(dict(name="fsl", value=0.74, 
                min=0, max=1, free=True,  #For the moment no more components
                latex=r"f_{\mathrm{short-lived}}"))
  pars.add(dict(name="fll", formula="(1-fprompt)*fphys2",
                latex=r"f_{\mathrm{long-lived}}")) #Needed if we use fwpv
  pars.add(dict(name="taul", value=0.1,
                min=0, max=10, free=True)) #For the moment no more components
  pars.add(dict(name="taus", value=0.01,
                min=0., max=10., free=True)) #For the moment no more components


  pars.add(dict(name="fwpv",
                formula="(1-fprompt)*(1-fphys2)",
                latex=r"f_{\mathrm{WPV}}"))
  #               # formula="(1-fprompt)*(1-fphys2)",
  #               # latex=r"f_{\mathrm{WPV}}"))
  # pars.add(dict(name="tau1", value=0.6,
  #               min=0.1, max=1., free=False))
  # pars.add(dict(name="tau2", value=0.6,
  #               min=0.1, max=1, free=False))
  # pars.add(dict(name="share", value=0.4,
  #               min=0.1, max=0.8, free=False))
  wpv.lock()
  # print(wpv)
  # exit()
  pars = pars + wpv #Warning
  print(pars)

  pars.add( #2 Gaussians
      dict(name="part1", formula="(1-f) * exp(-(1/2.) * (sigma1*sigma1) * (DM*DM))"))
  pars.add(
      dict(name="part2", formula="f  * exp(-(1/2.) * (sigma2*sigma2) * (DM*DM))"))
  pars.add(dict(name="dilution", formula="part1 + part2"))
  pars.add(dict(name="sigmaeff", formula="sqrt(-2*log(part1+part2))/DM"))

  # pars.add(
      # dict(name="part1", formula="(1-f2-f3) * exp(-(1/2.) * (sigma1*sigma1) * (DM*DM))"))

  # pars.add(
      # dict(name="part2", formula="f2  * exp(-(1/2.) * (sigma2*sigma2) * (DM*DM))"))

  # pars.add(
      # dict(name="part3", formula="f3  * exp(-(1/2.) * (sigma3*sigma3) * (DM*DM))"))
  #

  # pars.add(dict(name="dilution", formula="part1 + part2 + part3"))
  # pars.add(dict(name="sigmaeff", formula="sqrt(-2*log(part1+part2 + part3))/DM"))
  print(pars)

  timed = ipanema.ristra.allocate(np.float64(cdf[time].values))
  sigmatd = ipanema.ristra.allocate(np.float64(cdf[sigmat].values))
  if weight:
    weightd = ipanema.ristra.allocate(np.float64(cdf.eval(weight).values))
  else:
    weightd = timed/timed
  prob = 0 * timed
  res = ipanema.optimize(fcn, pars, fcn_args=(timed, prob, weightd),
                         method='minuit', tol=0.1, verbose=True)
  print(res)
  fpars = ipanema.Parameters.clone(res.params)
  for k, v in fpars.items():
    v.min = -np.inf
    v.max = +np.inf
    v.set(value=res.params[k].value, min=-np.inf, max=np.inf)
  
  #Functions to propagate uncs
  wexp = uncertainties.wrap(np.exp)
  wlog = uncertainties.wrap(np.log)
  wsqrt = uncertainties.wrap(np.sqrt)
  pars = res.params
  f = pars['f'].uvalue
  # sigma1 = pars['sigmap'].uvalue - pars['sigmapp'].uvalue * ((f) / (1 - f))**0.5
  # sigma2 = pars['sigmap'].uvalue + pars['sigmapp'].uvalue * ((1 - f) / (f))**0.5
  sigma1 = pars['sigma1'].uvalue
  sigma2 = pars['sigma2'].uvalue
  exp1 = wexp(-(1 / 2.) * (sigma1 * pars['DM'].value)**2)
  exp2 = wexp(-(1 / 2.) * (sigma2 * pars['DM'].value)**2)
  part1 = (1 - pars['f'].uvalue) * exp1
  part2 = (0 + pars['f'].uvalue) * exp2
  dilution = part1 + part2
  sigmaeff = wsqrt(-2 * wlog(part1 + part2)) / pars['DM'].value
  pars['sigmaeff'].set(value=sigmaeff.n, stdev=sigmaeff.s)
  pars['dilution'].set(value=dilution.n, stdev=dilution.s)
  # sigma_ave = np.mean(ipanema.ristra.get(sigmatd))
  sigma_ave = np.sum(ipanema.ristra.get(weightd*sigmatd))/(np.sum(ipanema.ristra.get(weightd)))
  pars.add(dict(name='sigmaAverage', value=sigma_ave, stdev=0))
  nevts = len(timed)
  nevts = uncertainties.ufloat(nevts, np.sqrt(nevts))
  pars.add(dict(name='nevts', value=nevts.n, stdev=nevts.s))
  print(f"Dilution:          {dilution:.2uL}")
  print(f"Sigma:             {sigmaeff:.2uL}")
  print(f"Average of sigmat: {sigma_ave}")
  print(f"Number of events:  {nevts:.2uL}")
  print("New set of parameters")
  print(pars)
  

  ##Plot
  timeh = ipanema.ristra.get(timed)
  # sigmath = ipanema.ristra.get(sigmatd)
  weighth = ipanema.ristra.get(weightd)
  # sigmath = np.float64(cdf[sigmat].values)
  probh = 0 * timeh


  fig, axplot, axpull = complot.axes_plotpull()
  hdata = complot.hist(timeh, bins=300, weights=weighth, density=False)
  # hdata = complot.hist(timeh, bins=100,  density=False)

  axplot.errorbar(hdata.bins, hdata.counts, yerr=hdata.yerr, xerr=hdata.xerr,
                  fmt=".k")

  proxy_time = ipanema.ristra.linspace(min(timeh), max(timeh), 5000)
  proxy_prob = 0 * proxy_time

  def pdf(params, time, prob, norm):
    #Function already has the params w/ weights taken into account:
    weight = time/time #Only not to give an error
    lkhd = fcn(params, time, prob, weight, norm=1)
    return norm * np.exp(-lkhd / 2)

  # plot subcomponents
  species_to_plot = ['fprompt', 'fll', 'fwpv']
  # species_to_plot = ['fprompt', 'fll']
  for icolor, pspecie in enumerate(species_to_plot):
    _color = f"C{icolor+1}"
    _label = rf"${fpars[pspecie].latex.split('f_{')[-1][:-1]}$"
    # print(_label)
    _p = ipanema.Parameters.clone(fpars)
    # print(_p)
    for f in _p.keys():
      # print("fraction", f)
      # _p[f].set(value=fpars[pspecie].value, min=-np.inf, max=np.inf)
      if f.startswith('f') and f != pspecie:
        # print(f"looking at {f}")
        if len(f) == 1:
          0
          # print('f as is')
        elif 'fphys' in f:
          0
          # print('fphys as is')
        elif 'fsl' in f:
          0
        else:
          _p[f].set(value=0, min=-np.inf, max=np.inf)
    # total_frac = fpars['fprompt'].value + fpars['fsl'].value + fpars['fll'].value + fpars['fwpv'].value
    total_frac = fpars['fprompt'].value + fpars['fll'].value + fpars['fwpv'].value
    # print(_p)
    # print(total_frac)
    # print("-----------------------------------------------")
    _prob = pdf(_p, proxy_time, proxy_prob, norm=hdata.norm)
    _prob = np.nan_to_num(_prob)
    # print(_prob, ipanema.ristra.sum(prob))
    axplot.plot(ipanema.ristra.get(proxy_time), _p[pspecie].value * ipanema.ristra.get(_prob) / total_frac,
                color=_color, linestyle='--', label=_label)
  # exit()
  # plot fit with all components and data
  _p = ipanema.Parameters.clone(res.params)
  _prob = pdf(_p, proxy_time, proxy_prob, norm=hdata.norm)
  axplot.plot(ipanema.ristra.get(proxy_time), _prob, color="C0",
              label=rf"full fit $ {sLL} < \sigma_t < {sUL}$")
  pulls = complot.compute_pdfpulls(ipanema.ristra.get(proxy_time), ipanema.ristra.get(_prob),
                                   hdata.bins, hdata.counts, *hdata.yerr)
  axpull.fill_between(hdata.bins, pulls, 0, facecolor="C0", alpha=0.5)

  # label and save the plot
  axpull.set_xlabel(r"$t$ [ps]")
  axplot.set_ylabel(r"Candidates")
  axpull.set_ylim(-1, 0.1)
  axplot.set_xlim(tLL, tUL)
  axpull.set_yticks([-5, 0, 5])
  # axpull.set_yticks([-2, 0, 2])
  # axpull.hlines(3, mLL, mUL, linestyles='dotted', color='k', alpha=0.2)
  # axpull.hlines(-3, mLL, mUL, linestyles='dotted', color='k', alpha=0.2)

  axplot.legend(loc="upper right", prop={'size': 8})
  os.makedirs(figs, exist_ok=True)
  fig.savefig(os.path.join(figs, "fit.pdf"))
  axplot.set_yscale("log")

  try:
    axplot.set_ylim(1e0, 1.5 * np.max(_prob))
  except:
    print("axes not scaled")

  fig.savefig(os.path.join(figs, "log_fit.pdf"))
  plt.close()

  return pars

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calibration of the decay time resolution.')
    parser.add_argument('--data', help='Input prompt data.')

    parser.add_argument('--wpv-shape', default='None',
                        help='File to the description of wpv shape')

    parser.add_argument('--time-range', default='None',
                        help='time range')

    parser.add_argument('--nbin', help='Bin of sigmat')

    parser.add_argument('--output-json', help='Location to save fit parameters')
    parser.add_argument('--output-plots', help='Location to create plots')

    args = vars(parser.parse_args())
    print(args)

    mode = "DsPi"
    # branches = ["time", "sigmat", "sig_sw", "gbw"] #TODO: In the future our sw8
    branches = ["time", "sigmat"]
    current_bin = int(args['nbin'])
    data = args['data']
    wpv_shape = args['wpv_shape']
    wpv = ipanema.Parameters.load(args['wpv_shape'])
    #From Kechen
    bines = [0., 0.023, 0.027, 0.030, 0.033, 0.036, 0.038, 0.042, 0.047, 0.062, 0.090]
    # weight = "sig_sw"
    # weight = "sig_sw" #Warning da muy mal con gbw8 -> Checkear wpv (?) TODO
    weight = False
    

    #DATAFRAME cut in bins
    df = uproot.open(data)["DecayTree"].pandas.df(branches=branches )
    sigma_range = bines[current_bin:current_bin+2]
    df = df.query(f"sigmat>{sigma_range[0]} & sigmat<{sigma_range[1]}")
    time_range = [float(i) for i in args["time_range"].split(",")]
    # time_range = [-3, 7.]  #Warning
    cut_time = True
    if cut_time:
      df.query(f"time>={time_range[0]} & time<={time_range[1]} ")

    pars = time_resolution(df, time_range, sigma_range, wpv, weight=weight, mode="DsPi", figs=args['output_plots'])
    pars.dump(args['output_json'])
    print(pars)



# vim: fdm=marker ts=2 sw=2 sts=2 sr et
