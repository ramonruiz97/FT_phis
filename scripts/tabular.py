# tabular
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

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
        description='Tagging Calibration for B+ and MC')
  parser.add_argument('--input-201516', 
                      help='Input data.')
  parser.add_argument('--input-2017', 
                      help='Input data.')
  parser.add_argument('--input-2018', 
                      help='Input data.')

  parser.add_argument('--output', 
						help='Location to save fit parameters')


  args = vars(parser.parse_args())
  print(args)
  
  pars_201516 = [ipanema.Parameters.load(p) for p in args['input_201516'].split(',')]
  pars_2017 = [ipanema.Parameters.load(p) for p in args['input_2017'].split(',')]
  pars_2018 = [ipanema.Parameters.load(p) for p in args['input_2018'].split(',')]
	#201516
  p0_201516_value = [p['p0'].value for p in pars_201516]
  p0_201516_stdev = [p['p0'].stdev for p in pars_201516]

  p1_201516_value = [p['p1'].value for p in pars_201516]
  p1_201516_stdev = [p['p1'].stdev for p in pars_201516]

	#2017
  p0_2017_value = [p['p0'].value for p in pars_2017]
  p0_2017_stdev = [p['p0'].stdev for p in pars_2017]

  p1_2017_value = [p['p1'].value for p in pars_2017]
  p1_2017_stdev = [p['p1'].stdev for p in pars_2017]

	#2018
  p0_2018_value = [p['p0'].value for p in pars_2018]
  p0_2018_stdev = [p['p0'].stdev for p in pars_2018]

  p1_2018_value = [p['p1'].value for p in pars_2018]
  p1_2018_stdev = [p['p1'].stdev for p in pars_2018]
  
  table = []
  table.append(r"\begin{tabular}{lllllll} ")
  table.append(r"\toprule & mc $b_s\rightarrow j/\psi \phi$ & mc $b_s\rightarrow j/\psi \phi (pt cut)$  & mc $b_s\rightarrow d_s\pi$ (gbw) & mc $b_s\rightarrow d_s \pi$  & syst\\ ")
  table.append(r"\midrule ")
  table.append(r"201516\\ ")
  table.append(r"\midrule ")
  table.append(rf"$p0$ & {unc.ufloat(p0_201516_value[0], p0_201516_stdev[0]):.2uL}& {unc.ufloat(p0_201516_value[1], p0_201516_stdev[1]):.2uL} & {unc.ufloat(p0_201516_value[2], p0_201516_stdev[2]):.2uL} &  {unc.ufloat(p0_201516_value[3], p0_201516_stdev[3]):.2uL}\\ ")
  table.append(rf"$p1$ & {unc.ufloat(p1_201516_value[0], p1_201516_stdev[0]):.2uL}& {unc.ufloat(p1_201516_value[1], p1_201516_stdev[1]):.2uL} & {unc.ufloat(p1_201516_value[2], p1_201516_stdev[2]):.2uL} &  {unc.ufloat(p1_201516_value[3], p1_201516_stdev[3]):.2uL}\\ ")
  table.append(r"\midrule ")
  table.append(r"2017\\ ")
  table.append(r"\midrule ")
  table.append(rf"$p0$ & {unc.ufloat(p0_2017_value[0], p0_2017_stdev[0]):.2uL}& {unc.ufloat(p0_2017_value[1], p0_2017_stdev[1]):.2uL} & {unc.ufloat(p0_2017_value[2], p0_2017_stdev[2]):.2uL} &  {unc.ufloat(p0_2017_value[3], p0_2017_stdev[3]):.2uL}\\ ")
  table.append(rf"$p1$ & {unc.ufloat(p1_2017_value[0], p1_2017_stdev[0]):.2uL}& {unc.ufloat(p1_2017_value[1], p1_2017_stdev[1]):.2uL} & {unc.ufloat(p1_2017_value[2], p1_2017_stdev[2]):.2uL} &  {unc.ufloat(p1_2017_value[3], p1_2017_stdev[3]):.2uL}\\ ")
  table.append(r"\midrule ")
  table.append(r"2018\\ ")
  table.append(r"\midrule ")
  table.append(rf"$p0$ & {unc.ufloat(p0_2018_value[0], p0_2018_stdev[0]):.2uL}& {unc.ufloat(p0_2018_value[1], p0_2018_stdev[1]):.2uL} & {unc.ufloat(p0_2018_value[2], p0_2018_stdev[2]):.2uL} &  {unc.ufloat(p0_2018_value[3], p0_2018_stdev[3]):.2uL}\\ ")
  table.append(rf"$p1$ & {unc.ufloat(p1_2018_value[0], p1_2018_stdev[0]):.2uL}& {unc.ufloat(p1_2018_value[1], p1_2018_stdev[1]):.2uL} & {unc.ufloat(p1_2018_value[2], p1_2018_stdev[2]):.2uL} &  {unc.ufloat(p1_2018_value[3], p1_2018_stdev[3]):.2uL}\\ ")
  table.append(r"\midrule ")
  table.append(r"$<\eta>$ & \multicolumn{4}{c}{0.402} \\" )
  table.append(r"\bottomrule ")
  table.append(r"\end{tabular}")

  print(table)
  with open(args['output'], "w") as tex_file:
    tex_file.write("\n".join(table))
  tex_file.close()







# vim: fdm=marker ts=2 sw=2 sts=2 sr noet
