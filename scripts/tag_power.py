# tag_power
#
#

__all__ = []
__author__ = ["Ramón Ángel Ruiz Fernández"]
__email__ = ["rruizfer@CERN.CH"]



import numpy as np
import ipanema
import complot
import uncertainties as unc
import argparse
import yaml
import os
import json

#Function by tagger and year only 
#By default the json in gitlab has all the years so one must identify the correct year
def get_pars(path, year):
  y = year if year != '201516' else '2016'
  data = json.load(open(str(path)))

  pars = ipanema.Parameters()
  pars_names = ["p0", "p1", "dp0", "dp1"]
  index = {
		        "2016" : 0,
		        "2017" : 1,
		        "2018" : 2	}
	
  for i, n in enumerate(pars_names):
    pars.add(dict(name = f"{n}",
                value = data['ResultSet'][index[y]]['Parameter'][i]["Value"],
                stdev = np.sqrt(data['ResultSet'][index[y]]['Parameter'][i]["Error"]**2  +
                            + data['ResultSet'][index[y]]['SystematicErrors'][0]['Values'][i]**2 +
                            + data['ResultSet'][index[y]]['SystematicErrors'][1]['Values'][i]**2 +
                            + data['ResultSet'][index[y]]['SystematicErrors'][2]['Values'][i]**2)
                ))

  for l, k in enumerate(pars_names):
    pars[k].correl = {}
    for m, k2 in enumerate(pars_names):
      pars[k].correl[k2] = data['ResultSet'][index[y]]['StatisticalCorrelationMatrix'][l][m]

  pars.add(dict(name="eta_bar",
					      value=data['Parameter'][0]["Value"],
                stdev = 0.))	


	  
  return pars

def get_pars_epm(path, year):
  line = []
  with open(path) as f:
    line = [line.split('=')[1].replace('\n', '') for line in f.readlines()]
  # data = np.loadtxt(open(str(path)))
  pars = ipanema.Parameters()
  pars.add(dict(name= 'p0',
					      value = float(line[1]) + float(line[2]),
							  stdev = float(line[4])
					      ))
  pars.add(dict(name= 'dp0',
					      value = 0.,
							  stdev = 0.
					      ))

  pars.add(dict(name= 'dp1',
					      value = 0.,
							  stdev = 0.
					      ))
  pars.add(dict(name= 'eta_bar',
					      value = float(line[1]),
							  stdev = 0.
					      ))

  pars.add(dict(name= 'p1', 
					      value = float(line[3]),
							  stdev = float(line[5])
					      ))
  pars['p0'].correl = {}
  pars['p0'].correl['p1'] = float(line[6])
  pars['p1'].correl = {}
  pars['p1'].correl['p0'] = float(line[6])
  return pars

def calibration(eta, dec, p):
  p = p.valuesdict()
  omega = (eta-p['eta_bar'])*(p['p1']+0.5*dec*p['dp1']) + (p['p0']+0.5*dec*p['dp0'])
  return omega
#
#
def tag_power(om, dec, w, norm):
  A = (1 + dec*(1.-2.*om))
  B = (1 - dec*(1.-2.*om))
  return np.sum(w*((A-B)**2)/(A+B)**2)/norm

def tag_power_comb(om_ss, dec_ss, #always ss first :)
									 om_os, dec_os,
									 w, norm):

  A = (1 + dec_os*(1.-2.*om_os))*(1 + dec_ss*(1.-2.*om_ss)) 
  B = (1 - dec_os*(1.-2.*om_os))*(1 - dec_ss*(1.-2.*om_ss)) 
  return np.sum(w*((A-B)**2)/(A+B)**2)/norm

def tag_power_err(eta, om,
									q, w, norm,
									pars):
  qD = q*q*(1. - 2. * om)
  deriv = np.array([qD, qD*(eta - pars["eta_bar"].value), 
													qD*0.5*np.sign(q), 
													qD*0.5*(eta-pars["eta_bar"].value)*np.sign(q)])
  pars_err = ipanema.Parameters.build(pars, ["p0", "p1", "dp0", "dp0"])
  cov = pars_err.cov()
  err = 0

  for i in range(len(pars_err.cov())):
    for j in range(len(pars_err.cov())):
	    err += np.sum(w*deriv[i])*np.sum(w*deriv[j]) * cov[i][j]
  err *= 16./norm**2
  return np.sqrt(err)

def comb_err(eta_ss, eta_os,
						 om_ss, om_os,
	           dec_ss, dec_os,
	           w, norm,
						 pars_ss, pars_os):
	#SS -> 0
	#OS -> 1
  D_ss = dec_ss * (1. - 2. *om_ss)
  D_os = dec_os * (1. - 2. *om_os)
  dw0 = (dec_ss * (D_ss + D_os)*(1.-D_os**2))/(1.+D_ss*D_os)**3
  dw1 = (dec_os * (D_ss + D_os)*(1.-D_ss**2))/(1.+D_ss*D_os)**3

  deriv_0 = np.array([dw0, dw0*(eta_ss-pars_ss['eta_bar'].value), 
										 dw0 * 0.5 * np.sign(dec_ss),
	                   dw0 * 0.5 * (eta_ss - pars_ss['eta_bar'].value)*np.sign(dec_ss)])

  deriv_1 = np.array([dw1, dw1*(eta_os-pars_os['eta_bar'].value), 
										dw1 * 0.5 * np.sign(dec_os),
	                   dw1 * 0.5 * (eta_os - pars_os['eta_bar'].value)*np.sign(dec_os)])

  pars_err = ipanema.Parameters.build(pars_ss, ["p0", "p1", "dp0", "dp0"])
  cov = pars_err.cov()
  err = 0

  for i in range(len(pars_err.cov())):
    for j in range(len(pars_err.cov())):
      err += np.sum(w*deriv_0[i])*np.sum(w*deriv_0[j]) * cov[i][j]

  pars_err = ipanema.Parameters.build(pars_os, ["p0", "p1", "dp0", "dp0"])
  cov = pars_err.cov()

  for i in range(len(pars_err.cov())):
    for j in range(len(pars_err.cov())):
      err += np.sum(w*deriv_1[i])*np.sum(w*deriv_1[j]) * cov[i][j]

  err *= 16./norm**2
  return np.sqrt(err)

  

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
        description='Tagging Calibration for B+ and MC')
  parser.add_argument('--calibration-SS', 
                      help='calibration-SS', required=False)
  parser.add_argument('--calibration-OS', 
                      help='calibration-OS', required=False)
  parser.add_argument('--calibration-IFT', 
                      help='calibration-IFT', required=False)
	
  parser.add_argument('--data', help="input data")


  parser.add_argument('--output', 
						help='Location to save fit parameters')

  args = vars(parser.parse_args())  

  
	#It is important to take into account Portability errors
  # SS_cal = [ipanema.Parameters.load(p) for p in args['calibraion_SS'].split(',')]
  # OS_cal = [ipanema.Parameters.load(p) for p in args['calibraion_OS'].split(',')]
  # IFT_cal = [ipanema.Parameters.load(p) for p in args['calibraion_IFT'].split(',')]

  years = ["201516", "2017", "2018"]
  # years = ["2018"]
  taggers = ["OS", "SSK", "IFT"]

  with open('config/tagger.yaml') as config:
	  config = yaml.load(config, Loader=yaml.FullLoader)
 
  eta_dict, q_dict = {}, {}
  for tagger in taggers:
	  eta_dict[tagger] = config[tagger]['branch']['eta']
	  q_dict[tagger] = config[tagger]['branch']['decision']

  branches = [config[tagger]['branch']['decision'] for tagger in taggers]
  branches += [config[tagger]['branch']['eta'] for tagger in taggers]
  branches += ["sigBsSW"]

	#Load samples
  data = {}
  for year in years:
	  # year = year if year!= "2016" else "201516"
	  data_tuple = f"tuples/Bs2JpsiPhi/{year}/v4r0.root"
	  data[year] = ipanema.Sample.from_root(data_tuple, branches=branches)

  
	# Load Calibration parameters
  dict_pars = {}
  for tagger in taggers:
    path =  f"calibrations_v4r0/ResultList-{tagger}-v4r0.json"
    for year in years:
      dict_pars[tagger+year] = get_pars(path, year)

  #To check -> #Warning	
  #Load Calibrations epm
  # dict_pars = {}
  # for tagger in taggers:
  #   for year in years:
  #     path =  f"test_epm_tagpower/2018/v4r0_{tagger}.py"
  #     dict_pars[tagger+year] = get_pars_epm(path, year)

  
  
  	
  # cov = dict_pars["OS2017"].cov()
  # print(len(dict_pars["OS2017"].cov()))
  # for i in range(len(cov)):
  #   for j in range(len(cov)):
  #     print(cov[i][j])
  # exit()
	# Main algorithm:
  exclude = {'OS' : 'SSK', 'SSK': 'OS'}
  for year in years:
    print("==============================================")
    print(f"Tagging power for year: {year}")
    print("==============================================")
    norm = data[year].df["sigBsSW"].sum()
    TP = {}
    omega, eta, errs, q = [], [], [], []
    for tagger in ["OS", "SSK"]:
      data_excl = data[year].df.query(f"{eta_dict[tagger]} > 0. & {eta_dict[tagger]} < 0.5 & {eta_dict[exclude[tagger]]} == 0.5")
      data_incl = data[year].df.query(f"{eta_dict[tagger]} > 0. & {eta_dict[tagger]} < 0.5 & {eta_dict[exclude[tagger]]} > 0.0 & {eta_dict[exclude[tagger]]} < 0.5")
      
      omegaExcl = calibration(data_excl[eta_dict[tagger]],
                              data_excl[q_dict[tagger]],
                              dict_pars[tagger+year])

      omegaIncl = calibration(data_incl[eta_dict[tagger]],
                              data_incl[q_dict[tagger]],
                              dict_pars[tagger+year])

      omega.append(omegaIncl.array)
      eta.append(data_incl[eta_dict[tagger]].array)
      q.append(data_incl[q_dict[tagger]].array)

      TP[tagger+year] = np.round(100*tag_power(
														              omegaExcl.array,
														              data_excl[q_dict[tagger]].array,
														              data_excl["sigBsSW"].array,
			                                    norm), 5)

      err = np.round(100*tag_power_err(data_excl[eta_dict[tagger]].array,
																		   omegaExcl.array,
																		   data_excl[q_dict[tagger]].array,
														          data_excl["sigBsSW"].array,
																		  norm,
																		  dict_pars[tagger+year]
									                     )
			                        , 5)
      errs.append(err)
			
      print(f"{tagger} {year} {TP[tagger+year]} +- {err}")
    
    TP["comb"+year] = np.round(100*tag_power_comb(omega[0], q[0],
                                                  omega[1], q[1],
                                                  data_incl["sigBsSW"].array, norm)
															, 5)

    err_comb = np.round(100*comb_err(eta[1], eta[0],
											               omega[1], omega[0],
											               q[1], q[0],
											              data_incl["sigBsSW"].array, norm,
											              dict_pars["SSK"+year], dict_pars["OS"+year])
												, 5) 

    print(f"OS+SSK {year} {TP['comb'+year]} +- {err_comb}")
    TP_comb = np.round(TP[f"OS{year}"] + TP[f"SSK{year}"] + TP[f"comb{year}"] ,3)
    total_err = np.round(np.sqrt(errs[0]**2 + errs[1]**2 + err_comb**2), 3)
    print("-------------------------------------------------------")
    print(f"Combination {year} {TP_comb} +- {total_err}")

    data_excl = data[year].df.query(f"{eta_dict['IFT']} > 0. & {eta_dict['IFT']} < 0.5")
    omegaIFT = calibration(data_excl[eta_dict["IFT"]],
                              data_excl[q_dict["IFT"]],
                              dict_pars["IFT"+year])

    TP["IFT"+year] = np.round(100*tag_power(
														              omegaIFT.array,
														              data_excl[q_dict["IFT"]].array,
														              data_excl["sigBsSW"].array,
			                                    norm), 5)

    err = np.round(100*tag_power_err(data_excl[eta_dict["IFT"]].array,
																		   omegaIFT.array,
																		   data_excl[q_dict["IFT"]].array,
														          data_excl["sigBsSW"].array,
																		  norm,
																		  dict_pars["IFT"+year]
									                     )
			                        , 5)

    print(f"{'IFT'} {year} {TP['IFT'+year]} +- {err}")
	

	








# vim: fdm=marker ts=2 sw=2 sts=2 sr noet
