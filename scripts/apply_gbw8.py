import numpy as np
import pandas
import argparse, yaml
from hep_ml import reweight
from hep_ml.metrics_utils import ks_2samp_weighted
import uproot3 as uproot
import matplotlib as plt

with open(r"/scratch48/ramon.ruiz/Inclusive_Tagging/tagging-for-phi_s/config/branches.yaml") as file:
  BRANCHES = yaml.load(file, Loader=yaml.FullLoader)


def gbWeight(
    original_file, original_tree, original_cut, 
    target_file, target_tree, target_cut, 
    output_file, output_tree,  mode, year, version):
    """
    Original branches -> Branches to keep in Tree 
    Target branches -> Branches to keep in  Tree
    ob -> Branches to rw8
    tb -> Branches to rw8
    of -> Original file w/ original branches
    tf -> Target file w/ target branches
    """

    original_weight = "sigBsSW"
    target_weight = "sigBsSW"
    if "Bu" in mode:
      original_weight = "sigBuSW"
    if "MC_Bs2DsPi" in mode:
      original_weight = "time/time"
    
    #Parameters gbw8
    nEstimators =  20
    Rate =  0.1
    nDepth = 3

    original_branches = list(BRANCHES[f'{mode}'].keys())
    ob =  ["B_PT", "nPVs", "nTracks", "B_ETA"] #Warning Gitanada
    if "Bs2DsPi" in mode:
      ob =  ["B_PT", "nPVs", "nTracks", "BETA"] #Warning Gitanada
    tb = ["B_PT", "nPVs", "nTracks", "B_ETA"]
    target_branches = list(BRANCHES['Bs2JpsiPhi'].keys())

    if "Prompt" in original_file: 
        original_branches = list(BRANCHES['Bs2DsPi_Prompt'].keys())
        target_branches = list(BRANCHES['Bs2DsPi'].keys())
        target_branches.append('gbw')
        ob =  ["B_IPCHI2_OWNPV", "B_PT", "nLong"] #Warning Gitanada
        tb = ob


    of = uproot.open(original_file)[original_tree].pandas.df(branches=original_branches, flatten=False)
    tf = uproot.open(target_file)[target_tree].pandas.df(branches=target_branches, flatten=False)

    if 'Bs2DsPi' in original_file:
      #This also works for DsPi_Prompt bcs the target file
      #DsPi data is already cut in B_PT
      of = of.query('(B_PT>2e3) & B_PVFitDs_M_1 > 5300 & B_PVFitDs_M_1 < 5600')
      print(f'Cut applied on {original_file}')

    if ("SW" or "sw") not in original_weight:
      original_weight = "time/time"
    if ("Bu" in mode):
      original_weight = "sigBuOSCombinationSW"

    ow = of.eval(original_weight)
    tw = tf.eval(target_weight)
    
    print("Settings:\n")
    print(f'Branches used to weight on {original_file}', ob)
    print(f'weights original {original_file}', original_weight)
    print(f'Branches used to weight on {target_file}', tb)
    print(f'weights target {target_file}', target_weight)

    reweighter = reweight.GBReweighter(n_estimators=nEstimators, 
                                       learning_rate=Rate, max_depth=nDepth, 
                                       min_samples_leaf=1000, 
                                       gb_args={'subsample': 1.0})

    reweighter.fit(of[ob], tf[tb], original_weight=ow, target_weight=tw)
    gb_weights = reweighter.predict_weights(of[ob])

    # Add weights to tuple
    of['gbw'] = gb_weights
    with uproot.recreate(f"{output_file}") as f:
      _branches = {}
      for k, v in of.items():
          if 'int' in v.dtype.name:
              _v = np.int32
          elif 'bool' in v.dtype.name:
              _v = np.int32
          else:
              _v = np.float64
          _branches[k] = _v
      mylist = list(dict.fromkeys(_branches.values()))
      f[original_tree] = uproot.newtree(_branches)
      f[original_tree].extend(of.to_dict(orient='list'))
    return of['gbw'].values


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--original-file', help='The file name of origin file')
    parser.add_argument('--original-tree', help='The tree name of origin tree')
    parser.add_argument('--target-file',   help='Name of target file')
    parser.add_argument('--target-tree',   help='The tree name of origin tree')
    # parser.add_argument('--original-weight', default="",   help='Weight original')
    # parser.add_argument('--target-weight',    default= "", help='Weight Target')
    parser.add_argument('--year',          help='Year')
    parser.add_argument('--mode',          help='mode of target decay')
    parser.add_argument('--version',          help='mode of target decay')
    parser.add_argument('--output-file',   help='Name of the output file')
    parser.add_argument('--output-tree', default="DecayTree",   help='Name of the output file')
    parser.add_argument('--original-cut', default="",   help='Cut on origin')
    parser.add_argument('--target-cut',    default= "", help='Cut on target')
    args = parser.parse_args()
    print(args)
    gbWeight(**vars(args))
