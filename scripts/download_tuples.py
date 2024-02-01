DESCRIPTION = """
    This script downloads tuples from eos.
"""

__author__ = ['Ramon Ángel Ruiz Fernández']
__email__ = ['rruizfer@cern.ch']
__all__ = []


import os
import argparse


if __name__ == "__main__":
  p = argparse.ArgumentParser(description=DESCRIPTION)
  p.add_argument('--year', 
                   choices = ['201516', '2015', '2016', '2017', '2018'],
                   help='Year of data taking')

  p.add_argument('--mode', 
                   help='Full root file with huge amount of branches.')

  p.add_argument('--version', 
                 help='Full root file with huge amount of branches.')
  p.add_argument('--output', 
                 help='Full root file with huge amount of branches.')

  args = vars(p.parse_args())

  # Get the flags and that stuff
  v = args['version'].split("@")[0].split("bdt")[0]  # pipeline tuple version
  V = args['version'].replace('bdt', '')  # full version for phis-scq
  y = args['year']
  m = args['mode']
  output = args['output']
  EOSPATH = "/eos/lhcb/wg/B2CC/Bs2JpsiPhi-FullRun2" 
  local_path = f"{output}"

  if m=='Bs2DsPi':
      v = "v1r5"
      if y=="201516":
        temp1 = output.replace(".root", "_1.root")
        temp2 = output.replace(".root", "_2.root")
        eos_path_1 = f'{EOSPATH}/{v}/{m}/{m}_2015_selected_{v}.root'
        status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path_1} {temp1}")
        eos_path_2 = f'{EOSPATH}/{v}/{m}/{m}_2016_selected_{v}.root'
        status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path_2} {temp2}")
        os.system(f"hadd {output} {temp1} {temp2}")
        os.system(f"rm {temp1} & rm {temp2}")
      else:
        eos_path = f'{EOSPATH}/{v}/{m}/{m}_{y}_selected_{v}.root'
        status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path} {local_path}")

  if m=='MC_Bs2DsPi':
      v = "v1r0"
      if y=="201516": #TODO merge 2015 + 2016
        y="2016"
        eos_path = f'{EOSPATH}/{v}/{m}/{m}_{y}_str28r2_sim09h_ldst_pidcorrected_selected_preselected_{v}.root'
        status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path} {local_path}")
      if y=="2017":
        eos_path = f'{EOSPATH}/{v}/{m}/{m}_{y}_str29r2_sim09g_dst_pidcorrected_selected_preselected_{v}.root'
        print(eos_path)
        status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path} {local_path}")
      if y=="2018":
        eos_path = f'{EOSPATH}/{v}/{m}/{m}_{y}_str34_sim09g_dst_pidcorrected_selected_preselected_{v}.root'
        status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path} {local_path}")
      else:
        print(f"Unrecognized years for {m} and {v}")
        
  elif m=='Bu2JpsiKplus': #Esto es una gitanada pero de momento Bu solo v1r1
      # eos_path = f'{EOSPATH}/v1r1/{m}/{y}/{m}_{y}_selected_bdt_sw_v1r1.root'
      # status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path} {local_path}")
      if y=="201516":
        temp1 = output.replace(".root", "_1.root")
        temp2 = output.replace(".root", "_2.root")
        eos_path_1 = f'{EOSPATH}/tuples/v4r0/{m}/2015/v4r0_{m}_2015_selected.root'
        status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path_1} {temp1}")
        eos_path_2 = f'{EOSPATH}/tuples/v4r0/{m}/2016/v4r0_{m}_2016_selected.root'
        status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path_2} {temp2}")
        os.system(f"hadd {output} {temp1} {temp2}")
        os.system(f"rm {temp1} & rm {temp2}")
      else:
        eos_path = f'{EOSPATH}/tuples/v4r0/{m}/{y}/v4r0_{m}_{y}_selected.root'
        status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path} {local_path}")

  elif m=='MC_Bu2JpsiKplus': #Esto es una gitanada pero de momento Bu solo v1r1
      # eos_path = f'{EOSPATH}/v1r1/{m}/{y}/{m}_{y}_selected_bdt_sw_v1r1.root'
      # status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path} {local_path}")
      if y=="201516":
        y="2016"
      eos_path = f'{EOSPATH}/tuples/v4r0/{m}/{y}/v4r0_{m}_{y}_selected.root'
      status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path} {local_path}")
  
  elif m=='Bs2JpsiPhi':
      if y=="201516":
        y="2016"
      path = f"/scratch49/marcos.romero/sidecar/{y}/{m}/v4r0_sWeight.root"
      status = os.system(f"cp {path} {local_path}")
      # eos_path = f'{EOSPATH}/{v}/{m}/{y}/{m}_{y}_selected_bdt_sw_{v}.root'
      # status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path} {local_path}")

  elif m=='MC_Bs2JpsiPhi':
      if y=="201516":
        temp1 = output.replace(".root", "_1.root")
        temp2 = output.replace(".root", "_2.root")
        eos_path_1 = f'{EOSPATH}/tuples/v4r0/{m}/2015/v4r0_{m}_2015_selected.root'
        status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path_1} {temp1}")
        eos_path_2 = f'{EOSPATH}/tuples/v4r0/{m}/2016/v4r0_{m}_2016_selected.root'
        status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path_2} {temp2}")
        os.system(f"hadd {output} {temp1} {temp2}")
        os.system(f"rm {temp1} & rm {temp2}")
      else:
        path = f"/scratch49/marcos.romero/sidecar/{y}/{m}/v4r0_sWeight.root"
        status = os.system(f"cp {path} {local_path}")
        
  else:
      print("Mode not recognized")

      # eos_path = f'{EOSPATH}/{v}/{m}/{y}/{m}_{y}_selected_bdt_{v}.root'
      # status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path} {local_path}")

  # if status:
  #     print(f"Trying wo version in name")
  #     eos_path = f'{EOSPATH}/{v}/{m}/{y}/{m}_{y}_selected_bdt_sw.root'
  #     status = os.system(f"xrdcp -f root://eoslhcb.cern.ch/{eos_path} {local_path}")









