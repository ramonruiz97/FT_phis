# Test PID
#
#

__all__ = []
__author__ = ["Ramón Ángel Ruiz Fernández"]
__email__ = ["rruizfer@CERN.CH"]


import uproot3 as uproot
from ipanema import Sample, Parameters, Parameter
import complot
import matplotlib.pyplot as plt
import numpy as np
import argparse

from utils.helpers import trigger_scissors
from utils.plot import mode2tex, get_range, watermark, mode_tex, get_var_in_latex


tuple_Bs = "tuples/Bs2JpsiPhi/2017/v4r0.root"
# tuple_Bu = "tuples/Bu2JpsiKplus/2017/Bu2JpsiKplus_2017_v4r0_ipatia_combined_sweighted.root"
tuple_Bu = "tuples/Bu2JpsiKplus/2017/v4r9_gbweightedOS.root"

bs_branches = ["sigBsSW", "nTracks", "nPVs", "B_PT", "B_ETA"]
bu_branches = ["sigBuOSCombinationSW", "gbw", "nTracks", "nPVs", "B_PT", "B_ETA"]
weight_bu = ["gbw", "sigBuSW*gbw"]
# bu_branches += ["OSCombination_TAGDEC", "OSCombination_TAGETA"]

bs = Sample.from_root(f"{tuple_Bs}", branches=bs_branches).df
bu = Sample.from_root(f"{tuple_Bu}", branches=bu_branches).df

weight_bs = "sigBsSW"
weight_bu = ["gbw", "sigBuOSCombinationSW*gbw"]

branches_plot = ["nTracks", "B_PT", "B_ETA"]  

for j in branches_plot:
  x, y, pull = complot.compare_hist(bs.eval(f'{j}'), bu.eval(f"{j}"),
																	 bs.eval(weight_bs), bu.eval(weight_bu[0]),
																	 range=[bs[f"{j}"].min(), bu[f"{j}"].max()], bins=60, density=True)
  
  x2, y2, pull2 = complot.compare_hist(bs.eval(f'{j}'), bu.eval(f"{j}"),
																	 bs.eval(weight_bs), bu.eval(weight_bu[1]),
																	 range=[bs[f"{j}"].min(), bu[f"{j}"].max()], bins=60, density=True)

  fig, axplot, axpull = complot.axes_plotpull()


  axplot.fill_between(y.bins,y.counts,
											step="mid", facecolor='none',
										 edgecolor='C1', alpha=0.7, #hatch='xxx',
											label=f"gbw")

  axplot.fill_between(y2.bins,y2.counts,
											step="mid", facecolor='none',
										  edgecolor='C2', alpha=0.7, #hatch='xxx',
											label=f"sw*gbw")


  axplot.errorbar(x.bins, x.counts,
									yerr=x.yerr[0], xerr=x.xerr,
									color='k', fmt='.',
                  label=f"sWeighted Bs")

  axpull.fill_between(y.bins, x.counts/y.counts, 1, facecolor='C1')
  axpull.fill_between(y.bins, x.counts/y2.counts, 1, facecolor='C2')

  axpull.set_xlabel(f"{j}")
  axpull.set_ylabel(r"$N_{data}/N_{MC}$")
  axpull.set_ylim(-3, 4)
  axpull.set_yticks([-3 ,1, 4])
  axplot.legend()
  fig.savefig(f"tag_plots/{j}.pdf")


# vim: fdm=marker ts=2 sw=2 sts=2 sr noet
