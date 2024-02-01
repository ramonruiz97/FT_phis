# ss_calibration
__author__ = ["Ramón Ángel Ruiz Fernández",
              "Marcos Romero Lamas",
              "Veronika Georgieva Chobanova"]


import numpy as np
import ipanema
import complot
import argparse
import hjson
import yaml
import os
from uncertainties import unumpy as unp

from scipy.stats import beta
from scipy.special import erf


# import analysis.badjanak as badjanak

from ipanema import uncertainty_wrapper
from ipanema.confidence import get_confidence_bands
# from ipanema.confidence.propagator import fast_jac, fast_hesse
# from scipy.interpolate import interp1d
from ipanema import extrap1d

#Taking Kernels
ipanema.initialize('cuda', 1)

prog = ipanema.compile(open('scripts/kernel.cu').read())

pdf = prog.ss_calibration
pdf_plot = prog.plot_ss
#Old acceptance: TODO, create a keyword for using this or not
# pdf_china = prog.ss_calibration_china

#scq style
# badjanak.config['knots'] = [0.3, 0.91, 1.96, 9.0]
# knots = [0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 12.0]
# badjanak.config['knots'] = knots
# badjanak.get_kernels(True)

#Translation DsPi style
from scripts.spline import print_spline
# from scripts.DsPiSpline import print_spline

def calib(eta, p0, p1, eta_bar):
  return p0 + p1*(eta-eta_bar)

def calib_2(eta, eta_bar, p):
  return p[0] + p[1]*(eta-eta_bar)

def omega_plot(time, sigma, q, id, prob,
	        G=0.6, DG=0.007, DM=17.6, sigma_0 = 0., sigma_1 = 1.,
					tLL=0.3, tUL=15., omega=0.5,
	        **coeffs):

  knots = [0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 12.0, 15.0]
  coeffs_d = ipanema.ristra.allocate(np.float64(
              print_spline(knots[1:-1], list(coeffs.values()), interpolate=1)
							))
  pdf_plot(time, sigma, q, id, prob, coeffs_d,
		  np.float64(G), np.float64(DG), np.float64(DM),
			np.float64(omega), 
			np.float64(sigma_0), np.float64(sigma_1),
		  np.float64(tLL), np.float64(tUL),
			np.int32(len(prob)), global_size=(len(prob),))

  return ipanema.ristra.get(prob)

def beta_central_interval_efficiency(k, N):
    # Returns the efficiency and its 1 sigma confidence interval of k events out of N passing a selection
    e = k / N
    nsigma = 1

    conf = erf(nsigma / np.sqrt(2))
    aa = k + 1
    bb = N - k + 1

    upper = beta.ppf((1 + conf) / 2, aa, bb)
    lower = beta.ppf((1 - conf) / 2, aa, bb)

    return [k / N, e - lower, upper - e]


def model(x, prob,
					p0=0, dp0=0, p1=1, dp1=0, p2=0, dp2=0, eta_bar=0.5,
	        G=0.6, DG=0.007, DM=17.6, sigma_0 = 0., sigma_1 = 1.,
					tLL=0.3, tUL=15.,
	        **coeffs):
  
  # coeffs_d = ipanema.ristra.allocate(np.float64(
         # badjanak.coeffs_to_poly(list(coeffs.values()))
  # ))
  # print(ipanema.ristra.get(coeffs_d))
  knots = [0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 12.0, 15.0]
  coeffs_d = ipanema.ristra.allocate(np.float64(
              print_spline(knots[1:-1], list(coeffs.values()), interpolate=1)
							))
	
  pdf(x, prob, coeffs_d,
		  np.float64(G), np.float64(DG), np.float64(DM),
			np.float64(p0), np.float64(dp0),
			np.float64(p1),np.float64(dp1),
			np.float64(p2), np.float64(dp2),
			np.float64(eta_bar),
			np.float64(sigma_0), np.float64(sigma_1),
		  np.float64(tLL), np.float64(tUL),
			np.int32(len(prob)), global_size=(len(prob),))

  return ipanema.ristra.get(prob)

# def model_china(time, sigmat, id, q, eta,
# 					prob,
# 					p0=0, dp0=0, p1=1, dp1=0, p2=0, dp2=0, eta_bar=0.5,
# 	        G=0.6, DG=0.007, DM=17.6, sigma_0 = 0., sigma_1 = 1.,
# 					tLL=0.3, tUL=15., a = 1., b=1., n=1., p0_epm=1.):
#    
# 	
#   pdf_china(time, sigmat, 
# 			      id, q, eta,
# 		      	prob, 
# 			      np.float64(a), np.float64(b), np.float64(n),
# 		        np.float64(G), np.float64(DG), np.float64(DM),
# 			      np.float64(p0), np.float64(dp0),
# 			      np.float64(p1),np.float64(dp1),
# 			      np.float64(p2), np.float64(dp2),
# 			      np.float64(eta_bar),
# 			      np.float64(sigma_0), np.float64(sigma_1),
# 		        np.float64(tLL), np.float64(tUL),
# 			      np.int32(len(prob)), global_size=(len(prob),))
# 	
# 	
#
#   return ipanema.ristra.get(prob)



def calibration_ss(data, mode, tagger, time_res, order="first", weight=False, figs=False):
  
  with open('config/tagger.yaml') as config:
    config = yaml.load(config, Loader=yaml.FullLoader)
	
  cut = config[tagger]['cut']
  q_b = config[tagger]['branch']['decision']
  q_eta = config[tagger]['branch']['eta']
	
  SECOND = True if "second" in order else False
  
	#Parameters for the fit
  pars = ipanema.Parameters()
  pars.add(dict(name='p0', value=0.44, min=0.0, max=.6, free=True))
  pars.add(dict(name='p1', value=0.64, min=0.4, max=1.3, free=True))
  pars.add(dict(name='dp0', value=0., min=-1, max=1, free=False))
  pars.add(dict(name='dp1', value=0., min=-1, max=1, free=False))
  pars.add(dict(name='p2', value=0., min=0.0, max=1.0, free=SECOND))
  pars.add(dict(name='dp2', value=0., min=-1, max=1., free=SECOND))
  pars.add(dict(name='eta_bar', value=0.5, min=-1, max=1, free=False))
  

	#Time resolution calibration -> we should do from a json
  # pars.add(dict(name='sigma_0', value=time_res['q0'].value, min=0., max=0.7, free=True))
  # pars.add(dict(name='sigma_1', value=time_res['q1'].value, min=0.5, max=1.5, free=True))
	#warning
  pars.add(dict(name='sigma_0', value=0.0029, min=0., max=0.7, free=False))
  pars.add(dict(name='sigma_1', value=1.0847, min=0.5, max=1.5, free=False))

  
	#Coefficients of time acc
	#SCQ style -> Bs2JpsiPhi param
  # pars.add(dict(name=f'c0', value=1., min=0., max=100., free=True)) #Convention
  # pars.add(dict(name=f'c1', value=1., min=0., max=100., free=True))
  # pars.add(dict(name=f'c2', value=1., min=0., max=100., free=True))
  # pars.add(dict(name=f'c3', value=1., min=0., max=100., free=True))
  # pars.add(dict(name=f'c4', value=1., min=0., max=100., free=True))
  # pars.add(dict(name=f'c5', value=1., min=0., max=100., free=False))

	#DsPi style -> Bs2DsPi time acc param
  pars.add(dict(name=f'c1', value=0.325, min=0., max=100., free=True)) #Convention
  pars.add(dict(name=f'c2', value=0.468, min=0., max=100., free=True))
  pars.add(dict(name=f'c3', value=0.779, min=0., max=100., free=True))
  pars.add(dict(name=f'c4', value=0.958, min=0., max=100., free=True))
  pars.add(dict(name=f'c5', value=1.107, min=0., max=100., free=True))
  pars.add(dict(name=f'c6', value=1.33, min=0., max=100., free=True))
  pars.add(dict(name=f'c7', value=1., min=0., max=100., free=False))
  pars.add(dict(name=f'c8', formula = f'c7 + (c6-c7)*({data.df["time"].max()}-12)/(3-12)'))
  # pars.add(dict(name=f'c8', value=0.9, min=0., max=100., free=True))

	#China style -> Old timeacc param
  # pars.add(dict(name=f'a', value=1.0691, min=0.7, max=1.2, free=True)) #Convention
  # pars.add(dict(name=f'b', value=0.021582, min=0.01, max=1., free=True))
  # pars.add(dict(name=f'n', value=2.3792, min=2., max=3., free=True))


  pars.add(dict(name='G', value=0.6596, min=0.60, max=0.73, free=True))
  pars.add(dict(name='DG', value=0.085, min=0.079, max=0.091, free=True))
  pars.add(dict(name='DM', value=17.766, min=15, max=19., free=True))

  pars.add(dict(name='tLL', value=data.df['time'].min(), free=False))
  pars.add(dict(name='tUL', value=data.df['time'].max(), free=False))
  # pars.add(dict(name='tLL', value=0.4, free=False))
  # pars.add(dict(name='tUL', value=15, free=False))
  pars.add(dict(name="p0_epm", formula="p0-eta_bar", free=False))
  pars.add(dict(name="tau", formula="1/G", free=False))
	#TODO: WARNING
  # pars.add(dict(name='tLL', value=0.3, free=False))
  # pars.add(dict(name='tUL', value=15., free=False))

	
  q = data.df[q_b]
  eta = data.df[q_eta]
	#Yuehong factor
  corr = np.sum(data.df.eval(weight))/np.sum(data.df.eval(f"{weight}*{weight}"))
	
  data.df['q'] = q
  data.df['eta'] = eta
  data.df['weight'] = data.df.eval(weight)
  b_id = 'B_TRUEID' if 'MC' in mode else 'B_ID'
  data.df['id'] = data.df[f'{b_id}']
  data.df['weight'] = data.df.eval(f"{weight}*{corr}")

  data.chop(cut) #Warning -> I have changed cut of the tagger -> CHECK
  data.chop("B_PT > 2000")
  print("Number of events: ", data.df.shape[0])
  print("Number of weighted events: ", np.sum(data.df.eval(f'{weight}')))
  print("Number of weighted events corrected: ", np.sum(data.df['weight']))


  data.allocate(weight='weight')
  data.allocate(x=['q', 'eta', 'id', 'time', 'sigmat'])
  data.allocate(prob="0*weight")

  # pars['eta_bar'].set(value=data.df.eval(f'weight*eta').sum()/data.df.eval('weight').sum())
  pars['eta_bar'].set(value=0.4239)
  print("Parameters before the fit")
  print(pars)
  # exit()

  def fcn(pars, data):
    p = pars.valuesdict()
    prob = model(data.x, data.prob, **p) #already in cpu

    chi2 = -2.*np.log(prob)*ipanema.ristra.get(data.weight)

    chi2_gauss = 0
    #Time Parameters:
    chi2_gauss += (p['DG'] - data.time_res['DG'].value)**2/data.time_res['DG'].stdev**2
    chi2_gauss += (p['DM'] - data.time_res['DM'].value)**2/data.time_res['DM'].stdev**2
    # chi2_gauss += (1/p['G'] - 1.516)**2/0.006**2
    chi2_gauss += (1/p['G'] - 1.520)**2/0.005**2

		#Time resolution:
    cov = np.asmatrix(data.time_res.cov(['q0', 'q1']))
    covInv = cov.I
		
    diff = np.matrix([
			                p['sigma_0'] - data.time_res['q0'].value,
		                  p['sigma_1'] - data.time_res['q1'].value
		              ])
    chi2_gauss += np.dot(np.dot(diff, covInv), diff.T).item(0) #warning
    # print(chi2_gauss)

    #Warning check:
    chi2_gauss  /= len(chi2)
    return chi2 + chi2_gauss
	
	
  # def fcn_china(pars, time, sigmat, q, id, eta, prob, weight):
  #   p = pars.valuesdict()
  #   num = model_china(time, sigmat, id, q, eta, prob, **p)
  #   # den = np.trapz(num, ipanema.ristra.get(time))
  #   prob = num
  #   chi2 = -2.*np.log(prob)*ipanema.ristra.get(weight)
  #
  #   chi2_gauss = 0
  #   #Time Parameters:
  #   chi2_gauss += (p['DG'] - data.time_res['DG'].value)**2/data.time_res['DG'].stdev**2
  #   chi2_gauss += (p['DM'] - data.time_res['DM'].value)**2/data.time_res['DM'].stdev**2
  #   chi2_gauss += (1/p['G'] - 1.516)**2/0.006**2
  #
  #   cov = np.asmatrix(data.time_res.cov(['q0', 'q1']))
  #   covInv = cov.I
		# 
  #   diff = np.matrix([
		# 	                p['sigma_0'] - data.time_res['q0'].value,
		#                   p['sigma_1'] - data.time_res['q1'].value
		#               ])
  #
		# #TODO Check if adding corr is needed -> Covariance matrix
  #   chi2_gauss += np.dot(np.dot(diff, covInv), diff.T)[0][0] #warning
  #
  #   chi2_gauss  /= len(chi2)
  #   return chi2 + chi2_gauss
	
  def fcn_plot(pars, time, sigma, q, id, prob, weight):
    p = pars.valuesdict()
    prob = omega_plot(time, sigma, q, id, prob, **p)
    return -2.*np.log(prob)*ipanema.ristra.get(weight)

  #TODO
  res = ipanema.optimize(fcn, pars, fcn_args=(data,),
                         method='minuit', verbose=True, timeit=True, policy='filter')
	
  print(res)

  # timed = ipanema.ristra.allocate(np.float64(data.df["time"].values))
  # sigmatd = ipanema.ristra.allocate(np.float64(data.df["sigmat"].values))
  # qd = ipanema.ristra.allocate(np.float64(data.df["q"].values))
  # idd = ipanema.ristra.allocate(np.float64(data.df["id"].values))
  # etad = ipanema.ristra.allocate(np.float64(data.df["eta"].values))
  # prob = 0*timed
  # weightd = ipanema.ristra.allocate(np.float64(data.df["weight"].values))
  #
  # res = ipanema.optimize(fcn_china, pars, fcn_args=(timed, sigmatd, qd, idd, etad, prob, weightd),
  #                        method='minuit', verbose=True)

	
  print(res.params)
  res.params.dump("forKechen.json")
  exit()

  # figs = True	
  if figs:
    os.makedirs(figs, exist_ok=True)
    fig, axplot = complot.axes_plot()

    x = np.linspace(0, 0.6, 800)
    names = ['p0', 'p1', 'eta_bar']
    pars_plot = ipanema.Parameters.build(res.params, names)
    # pars_plot = ipanema.Parameters.build(pars, names)
    y = calib(x, **pars_plot.valuesdict())

		#FT calib approximation sign(A) to correct maybe not enough for Bs
  #   nbins = 10
  #   eta_arr = np.array(data.df['eta'].array)
  #   sorted_eta = sorted(eta_arr)
  #   splitting  = np.array_split(sorted_eta, nbins)
  #   binning    = np.array([[s[0], s[-1]] for s in splitting])
		# #first try wo acceptance and resolutions maybe for the plot is enough -> EPM, ftcalib style
  #   def mixing_asymetry(time, DM, DG):
  #     return np.cos(DM*time) / np.cosh(0.5*DG*time)
		#
  #   Amix = mixing_asymetry(data.df['time'], pars['DM'].value, pars['DG'].value)
  #   pollution = 0.5 * (1-np.abs(Amix))
		#
  #   def true_mistag(eta_lo, eta_hi):
  #     data.df['pollution'] = pollution
  #     data.df['Amix'] = Amix
  #     data.df['prod_flav'] = data.df.eval('id*Amix/abs(Amix*id)')
		#
  #     eta_cut = f"(eta >= {eta_lo} & (eta < {eta_hi}))"
  #     right_cut = f"q == prod_flav & {eta_cut}"
  #     wrong_cut = f"q != prod_flav & {eta_cut}"
  #     # print(data.df.query(right_cut))
  #     # print(data.df.query(wrong_cut))
		#
  #     binweights = data.df.query(eta_cut)['weight']
  #     binpollution = data.df.query(eta_cut)['pollution']
  #     sumW = np.sum(binweights)
		#
  #     oscDsq = np.sum(binweights * (1 - 2 * binpollution)**2) / sumW
  #     asym = np.sum(data.df.query(right_cut)['weight'].array * (1-2*data.df.query(right_cut)['pollution'].array))
  #     asym -= np.sum(data.df.query(wrong_cut)['weight'].array * (1-2*data.df.query(wrong_cut)['pollution'].array))
  #     # print(oscDsq)
  #     # print(sumW)
  #     # print(asym)
		#
  #     asym /= (oscDsq*sumW)
		#
  #     effnum = (sumW**2 / np.sum(binweights**2)) * oscDsq
  #     # print(asym)
  #     Nw = 0.5 * (1 - asym) * effnum
  #     Nr = 0.5 * (1 + asym) * effnum
  #     return np.array([np.sum(data.df.query(eta_cut)['eta'].array * binweights) / sumW] + beta_central_interval_efficiency(Nw, Nw + Nr))
		#
  #   omega_true = np.array([true_mistag(lo, hi) for lo, hi in binning])
  #   # print("Binning:", binning)
  #   # print("eta: ", omega_true[:, 0])
  #   # print("Omega_True: ", omega_true[:, 1])
  #   # print("Errors: ", omega_true[:, 2])
  #   # print("Errors: ", omega_true[:, 3])
  #   # exit()
		# # Eta and omega ranges
  #   eta_range = [binning[0][0], binning[-1][1]]
  #   eta_interv = eta_range[1] - eta_range[0]
  #   eta_range[0] -= eta_interv * 0.02
  #   eta_range[1] += eta_interv * 0.02
  #   eta_range = [0., 0.55]
  #   omega_range = [0.2, 0.55]

	  #This is like EPM ft-calib style-> Omega_true from sign(A_mix)*decay_flavour = prod flavor  
    # axplot.errorbar(x = omega_true[:, 0], #eta 
				# 					  y = omega_true[:, 1], #omega true
				# 					  xerr = [omega_true[:, 0] - binning[:,0], binning[:,1]-omega_true[:,0]],
				# 					  yerr = [omega_true[:,2], omega_true[:, 3]],
				# 					  color      = 'k',
    #                 fmt        = '.',
    #                 label      = r"$\omega^{\mathrm{true}}$"
		  #               )

		#China style -> Suppose \omega constant in bin and fit the time distribution
    names = ["sigma_0", "sigma_1", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "G", "DG", "DM", "tLL", "tUL"]
    pars_china = ipanema.Parameters.build(res.params, names)
    pars_china.lock()
    pars_china.add(dict(name="omega", value=0.3, min=0., max=1., free=True))
    w = []
    err = []
    eta = []
    nbins = 10
    eta_arr = np.array(data.df['eta'].array)
    sorted_eta = sorted(eta_arr)
    splitting  = np.array_split(sorted_eta, nbins)
    binning    = np.array([[s[0], s[-1]] for s in splitting])
    for eta_lo, eta_hi in binning:
      eta_cut = f"(eta >= {eta_lo} & (eta < {eta_hi}))"

			#eta
      binweights = data.df.query(eta_cut)['weight']
      sumW = np.sum(binweights)
      eta.append(np.sum(data.df.query(eta_cut)['eta'].array * binweights) / sumW)

      df = data.df.query(eta_cut)
      timed = ipanema.ristra.allocate(np.float64(df["time"].values))
      sigmatd = ipanema.ristra.allocate(np.float64(df["sigmat"].values))
      qd = ipanema.ristra.allocate(np.float64(df["q"].values))
      idd = ipanema.ristra.allocate(np.float64(df["id"].values))
      prob = 0*timed
      weightd = ipanema.ristra.allocate(np.float64(df["weight"].values))

      # data.allocate(weight='weight')
      # data.allocate(x=['q', 'id', 'time', 'sigmat'])
      # data.allocate(prob="0*weight")
      omega = ipanema.optimize(fcn_plot, pars_china, fcn_args=(timed, sigmatd, qd, idd, prob, weightd),
		                             method='minuit', verbose=True).params["omega"]
      w.append(omega.value)
      err.append(omega.stdev)

		#Calibration curve
    axplot.plot(x, y, label="calibration")

    axplot.errorbar(x = eta, #eta
									  y = w, #omega true
									  xerr = [eta - binning[:,0], binning[:,1]-eta],
									  yerr = err,
									  color      = 'k',
                    fmt        = '.',
                    label      = r"$\omega$"
		                )


    eta_range = [binning[0][0], binning[-1][1]]
    eta_interv = eta_range[1] - eta_range[0]
    eta_range[0] -= eta_interv * 0.02
    eta_range[1] += eta_interv * 0.02
    eta_range = [0., 0.55]
    omega_range = [0.1, 0.55]
    


    names = ['p0', 'p1']
    pars_plot = ipanema.Parameters.build(res.params, names)
    eta_bar = res.params["eta_bar"].value * np.ones_like(x)
    #TODO Mandar Hessian
    y_unc = uncertainty_wrapper(lambda p: calib_2(x, eta_bar, p), pars_plot)
    yl, yh = get_confidence_bands(y_unc)
    yl2, yh2 = get_confidence_bands(y_unc, sigma=2)
    axplot.fill_between(x, yh, yl, color = "red", alpha =0.5)
    axplot.fill_between(x, yh2, yl2, color = "indianred", alpha =0.5)

    axplot.set_xlabel(rf"$\eta^{{{tagger}}}$")
    axplot.set_ylabel(f"$\omega(\eta^{{{tagger}}})$")
    axplot.set_ylim(*omega_range)
    axplot.set_xlim(*eta_range)
    axplot.legend()
    fig.savefig(os.path.join(figs, "calibration.pdf"))
    # fig.savefig("test.pdf")

  return res.params

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
        description='Tagging Calibration for B+ and MC')
  parser.add_argument('--data', help='Input data.')

  parser.add_argument('--weight', help='Weight to apply to data')

  parser.add_argument('--time-res', 
											help='Path parameters calibration resolution')

  parser.add_argument('--tagger', choices = ['SS', 'IFT'],
											help='Branch to calibrate')

  parser.add_argument('--mode', help='Input data.')
  parser.add_argument('--version', help='Input data.')

  parser.add_argument('--year', help='Year of data taking')

  parser.add_argument('--model', default='linear')

  parser.add_argument('--output-json', 
											help='Location to save fit parameters')

  parser.add_argument('--output-plots', 
											help='Location to create plots')

  args = vars(parser.parse_args())

  # tp = "/scratch48/ramon.ruiz/Inclusive_Tagging/tagging-for-phi_s/tuples/Bs2DsPi/2017/Bs2DsPi_2017_selected_bdt_v1r0_gbweighted.root"
  # data = ipanema.Sample.from_root(tp)
  # tagger =  "SS"
  # mode = "Bs2DsPi"
  # weight = 'sigBsSW'
  # weight = 'sigBsSW*gbw'
  # figs = False
  # order = "single"
  # year = "2017"
  # data.time_res =  ipanema.Parameters.load(f"config/SS_cal_{year}.json") 
  # print(args)
  # exit()

  tp = args["data"]
  data = ipanema.Sample.from_root(tp)
  weight = args["weight"]
  tagger = args["tagger"]
  order = args["model"]
  year = args["year"]
  mode = args["mode"]
  data.time_res = ipanema.Parameters.load(args["time_res"])


  time_res = ipanema.Parameters.load(args["time_res"])
  # print(time_res)
  # exit()
  
	#TODO: time_res is not needed 
  pars = calibration_ss(data, mode, tagger, time_res, order, weight, figs=args["output_plots"])

  pars.dump(args['output_json'])

# vim: fdm=marker ts=2 sw=2 sts=2 sr noet
