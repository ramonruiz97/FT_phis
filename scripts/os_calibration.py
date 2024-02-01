# Flavour Specific Calibration: Valid for MC and B+
__author__ = ["Ramón Ángel Ruiz Fernández",
              "Marcos Romero Lamas",
              "Veronika Georgieva Chobanova"]


import numpy as np
import ipanema
import complot
from uncertainties import unumpy as unp
import argparse
import yaml
import os

ipanema.initialize('cuda', 1)

prog = ipanema.compile(open('scripts/kernel.cu').read())
pdf = prog.os_calibration

def model(x, prob, p0=0, dp0=0, p1=1, dp1=0, p2=0, dp2=0, eta_bar=0.5):
  pdf(x, prob,
			np.float64(p0), np.float64(dp0),
			np.float64(p1),np.float64(dp1),
			np.float64(p2), np.float64(dp2),
			np.float64(eta_bar), 
			np.int32(len(prob)), global_size=(len(prob),))
  return ipanema.ristra.get(prob)


def fcn(pars, data):
  p = pars.valuesdict()
  prob = model(data.x, data.prob, **p)
  return -2. * np.log(prob) * ipanema.ristra.get(data.weight)

def calib(eta, p0, p1, p2, dp0, dp1, dp2, eta_bar, tag=1):
    result = 0
    result += (p0 + tag * 0.5 * dp0)
    result += (p1 + tag * 0.5 * dp1) * (eta - eta_bar)
    result += (p2 + tag * 0.5 * dp2) * (eta - eta_bar) * (eta - eta_bar)
    return result

def omega(eta, b_id, omega, p0=0, dp0=0, p1=1, dp1=0, p2=0, dp2=0, eta_bar=0.5):
  prog.calibrated_mistag(eta, b_id,  omega, 
                         np.float64(p0), np.float64(p1), np.float64(p2),
                         np.float64(dp0), np.float64(dp1), np.float64(dp2), 
                         np.float64(eta_bar),
                         np.int32(len(omega)), global_size=(len(omega),))
  return ipanema.ristra.get(omega)



def taggerCombination(qs, etas):
    pB = np.ones_like(qs[0])
    pBbar = np.ones_like(qs[0])
    
    #Events w/ omega > 0.5 < 0. do not contribute 
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

def os_calibration(data, mode, tagger, order="linear", 
                weight=False, calibrations=False, figs=False):

    with open('config/tagger.yaml') as config:
      config = yaml.load(config, Loader=yaml.FullLoader)
	
    cut = config[tagger]['cut']
    q_b = config[tagger]['branch']['decision']
    eta_b = config[tagger]['branch']['eta']

	
    # DP = False if "Bs" in mode else True 
    DP = True
    SECOND = True if "parabola" in order else False

    pars = ipanema.Parameters()
    pars.add(dict(name='p0', value=0.4, min=0.0, max=2.0, free=True))
    pars.add(dict(name='p1', value=0.9, min=0., max=2.1, free=True))
    pars.add(dict(name='dp0', value=0., min=-1, max=1, free=DP))
    pars.add(dict(name='dp1', value=0., min=-1, max=1, free=DP))
    pars.add(dict(name='p2', value=0., min=0.0, max=1.0, free=SECOND))
    pars.add(dict(name='dp2', value=0., min=-1, max=1., free=SECOND))

    pars.add(dict(name='eta_bar', value=0.5, min=-1, max=1, free=False))


    if calibrations:
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
      
      list_pars = [ipanema.Parameters.load(p) for p in calibrations]
      cal = [[p['p0'].value, p['p1'].value, p['eta_bar'].value] for p in list_pars]

      q_list = [np.float64(data.df[b]) for b in q_list]
      eta_list = [np.float64(data.df[b]) for b in eta_list]
      eta_cal = [np.float64(cal[i][0] + cal[i][1] * (eta_list[i]-cal[i][2])) for i in range(len(eta_list))]
      q, eta = taggerCombination(q_list, eta_cal)
      # q = data.df['OS_Combination_DEC'] #To check this branch from epm
      # eta = data.df['OS_Combination_ETA']

    else:
      q = data.df[q_b]
      eta = data.df[eta_b]


    # weight = "sw*gbw"
    print("Settings")
    print("Calibrations: ", calibrations)
    print("Branch Decision: ", q_b)
    print("Branch Mistag: ", eta_b)
    print("Model: ", order)
    corr = np.sum(data.df.eval(weight))/np.sum(data.df.eval(f"{weight}*{weight}"))
    data.df['q'] = q
    data.df['eta'] = eta
    # data.df['weight'] = data.df.eval(weight)
    b_id = 'B_ID_GenLvl' if 'MC' in mode else 'B_ID'
    data.df['id'] = data.df[f'{b_id}']
    data.df['weight'] = data.df.eval(f"{weight}*{corr}")
    print("Weight used to calibrate:", weight)


    data.chop(cut) #TODO: B_BKGCAT cuts?
    print("Cut applied to sample:", cut)

    data.allocate(weight='weight')
    data.allocate(x=['q', 'eta', 'id'])
    data.allocate(prob="0*weight")

    pars['eta_bar'].set(value=data.df.eval(f'weight*eta').sum()/data.df.eval('weight').sum())
    res = ipanema.optimize(fcn, pars, fcn_args=(data,),
							method='minuit', tol=0.05, verbose=False)

    if figs:
      os.makedirs(figs, exist_ok=True)
      fig, axplot = complot.axes_plot()

      x = np.linspace(0, 0.5, 800)
      pars_plot = ipanema.Parameters.clone(res.params)
      # pars['dp0'].set(value=0.)
      # pars['dp1'].set(value=0.)
      # eta_d = ipanema.ristra.allocate(data.df['eta'].array[0:1000])
      # id_d = ipanema.ristra.allocate(data.df['id'].array[0:1000])
      # prob_d = 0*eta_d
      # y = omega(eta_d, id_d, prob_d, **res.params.valuesdict())
      # x = ipanema.ristra.get(eta_d)

      #Calibration line
      y = calib(x, **pars_plot.valuesdict())
      _, edges = np.histogram(data.df['eta'], 10, weights=data.df['weight'])

      #w = W/(W+R)
      #WRONG + RIGHT
      total = complot.hist(data.df['eta'], 
                        bins=edges,
                    weights=data.df['weight'])
      #WRONG
      wrong = complot.hist(data.df.query('id*q<0')['eta'], bins=edges,
                    weights=data.df.query('id*q<0')['weight'])

      num = unp.uarray(wrong.counts, (wrong.yerr[0] + wrong.yerr[1])/2)
      den = unp.uarray(total.counts, (total.yerr[0] + total.yerr[1])/2)
      ratio = num / den
      axplot.errorbar(total.bins, unp.nominal_values(ratio),
                yerr=unp.std_devs(ratio), xerr=total.xerr, color='k', fmt='.',
                      label='w definition')


      axplot.plot(x, y, label="calibration")
      # axplot.plot(x, y, 'b.', label="calibration")
      axplot.set_xlabel("$\eta^{OS}$")
      axplot.set_ylabel("$\omega(\eta^{OS})$")
      axplot.set_ylim(0., 0.5)
      # axplot.legend()
      fig.savefig(os.path.join(figs, "calibration.pdf"))
      
    return res.params


if __name__ == '__main__':
  parser = argparse.ArgumentParser(
        description='Tagging Calibration for B+ and MC')
  parser.add_argument('--data', 
                      help='Input data.')

  parser.add_argument('--weight', 
                      help='Weight to apply to data')

  parser.add_argument('--tagger', 
						help='Branch to calibrate')

  parser.add_argument('--mode', 
                      help='Mode')

  parser.add_argument('--version', 
                      help='Version of tuples')

  parser.add_argument('--model', default='linear')

  parser.add_argument('--calibrations', default=False)

  parser.add_argument('--output-json', 
						help='Location to save fit parameters')

  parser.add_argument('--output-plots', 
						help='Location to create plots')

  args = vars(parser.parse_args())

  # tp = "/scratch48/ramon.ruiz/Inclusive_Tagging/tagging-for-phi_s/tuples/Bu2JpsiKplus/2017/Bu2JpsiKplus_2017_v1r0_tagged_OS.root"
  # data = ipanema.Sample.from_root(tp)
  # tagger =  "OS"
  # mode = "Bu2JpsiKplus"
  # weight = 'sw*gbw'
  # figs = False
  tp = args["data"]
  data = ipanema.Sample.from_root(tp, flatten=False)
  weight = args["weight"]
  tagger = args["tagger"]
  mode = args["mode"]
  order = args["model"]
  figs = args["output_plots"]
  if args['calibrations']:
    calibrations = args["calibrations"].split(",")
  else:
    calibrations = False

  pars = os_calibration(data, mode, tagger, order, weight, calibrations, figs)
  print(pars)
  pars.dump(args['output_json'])
	










