#Template for Bs taking SS, OS and IFT Calibrations


__author__ = ["Ramón Ángel Ruiz Fernández"]
__email__ = ["rruizfer@cern.ch"]





RootFile = "${input_tuple}"
TupleName = "DecayTree"
CalibrationMode = "Bs"
Nmax = -1


UseWeight = 1
BranchWeight = "${sweight}"


#OS Taggers
import ${OSComb_Calibration}
OS_Combination_Use = 1
OS_Combination_TypeDec = "Int_t"
OS_Combination_BranchDec = "OS_Combination_DEC"
OS_Combination_TypeProb = "Double_t"
OS_Combination_BranchProb = "OS_Combination_ETA"

#IFT taggers
import ${IFT_Calibration}
Inclusive_Use   = 1
Inclusive_TypeDec = "Int_t"
Inclusive_BranchDec = "B_IFT_InclusiveTagger_TAGDEC"
Inclusive_TypeProb = "Double_t"
Inclusive_BranchProb = "B_IFT_InclusiveTagger_TAGETA"

import ${SS_Calibration}
SS_Kaon_Use   = 1
SS_Kaon_TypeDec = "Int_t"
SS_Kaon_BranchDec = "B_SSKaonLatest_TAGDEC"
SS_Kaon_TypeProb = "Double_t"
SS_Kaon_BranchProb = "B_SSKaonLatest_TAGETA"


#SS + OS Combination
PerformOfflineCombination_OSplusSS = 1
OS_Combination_InComb = 1
SS_Kaon_InComb = 1
