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

def model_cal(x, prob, p0=0, dp0=0, p1=1, dp1=0, p2=0, dp2=0, eta_bar=0.5):
  pdf(x, prob,
			np.float64(p0), np.float64(dp0),
			np.float64(p1),np.float64(dp1),
			np.float64(p2), np.float64(dp2),
			np.float64(eta_bar), 
			np.int32(len(prob)), global_size=(len(prob),))
  return ipanema.ristra.get(prob)


def fcn(pars, data):
  p = pars.valuesdict()
  prob = model_cal(data.x, data.prob, **p)
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




def os_calibration(data, mode, correction, model, 
                   weight, figs=False):

    with open('config/tagger.yaml') as config:
      config = yaml.load(config, Loader=yaml.FullLoader)
	
    cut = config["IFT"]['cut']
    q_b = config["IFT"]['branch']['decision']
    eta_b = config["IFT"]['branch']['eta']
    cut = f"({cut})"

    DP = False
    if "full" in model:
      DP = True


    pars = ipanema.Parameters()
    pars.add(dict(name='p0', value=0.4, min=0.0, max=2.0, free=True))
    pars.add(dict(name='p1', value=0.9, min=0., max=2.1, free=True))
    pars.add(dict(name='dp0', value=0., min=-1, max=1, free=DP))
    pars.add(dict(name='dp1', value=0., min=-1, max=1, free=DP))
    pars.add(dict(name='p2', value=0., min=0.0, max=1.0, free=False))
    pars.add(dict(name='dp2', value=0., min=-1, max=1., free=False))

    pars.add(dict(name='eta_bar', value=0.402, min=-1, max=1, free=False))


    q = data.df[q_b]
    eta = data.df[eta_b]


    print("Settings")
    print("Branch Decision: ", q_b)
    print("Branch Mistag: ", eta_b)
    corr = np.sum(data.df.eval(weight))/np.sum(data.df.eval(f"{weight}*{weight}"))
    # corr = 1.
    data.df['q'] = q
    data.df['eta'] = eta
    b_id = "B_TRUEID"
    data.df['id'] = data.df[f'{b_id}']
    data.df['weight'] = data.df.eval(f"{weight}*{corr}")
    print("Weight used to calibrate:", weight)

    if ("cut" in correction) or ("Bs2DsPi" in mode):
      cut += " & (B_PT > 2e3)"

    cut += " & ( (B_BKGCAT == 0) | (B_BKGCAT == 50) | (B_BKGCAT==20) )"


    data.chop(cut) #TODO: B_BKGCAT cuts?
    print("Cut applied to sample:", cut)

    data.allocate(weight='weight')
    data.allocate(x=['q', 'eta', 'id'])
    data.allocate(prob="0*weight")


    res = ipanema.optimize(fcn, pars, fcn_args=(data,),
							method='minuit', tol=0.05, verbose=False)
    if figs:
      os.makedirs(figs, exist_ok=True)
      fig, axplot = complot.axes_plot()

      x = np.linspace(0, 0.5, 800)
      pars_plot = ipanema.Parameters.clone(res.params)

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
      axplot.set_xlabel("$\eta^{IFT}$")
      axplot.set_ylabel("$\omega(\eta^{IFT})$")
      axplot.set_ylim(0., 0.5)
      # axplot.legend()
      fig.savefig(os.path.join(figs, "calibration.pdf"))
      
    return res.params


if __name__ == '__main__':
  parser = argparse.ArgumentParser(
        description='Tagging Calibration for B+ and MC')
  parser.add_argument('--data', 
                      help='Input data.')

  parser.add_argument('--mode', 
                      help='Mode')

  parser.add_argument('--version', 
                      help='Version of tuples')

  parser.add_argument('--model')

  parser.add_argument('--corr')

  parser.add_argument('--output-json', 
						help='Location to save fit parameters')

  parser.add_argument('--output-plots', 
						help='Location to create plots')

  args = vars(parser.parse_args())

  tp = args["data"]
  data = ipanema.Sample.from_root(tp, flatten=False)
  mode = args["mode"]
  model = args["model"]
  correction = args["corr"]

  if ("Bs2DsPi" in mode) & ("gbw" in correction):
    weight = "gbw"
  else: 
    weight = "time/time"

  figs = args["output_plots"]
  

  pars = os_calibration(data, mode, correction, model, weight, figs)
  pars.dump(args['output_json'])
	










