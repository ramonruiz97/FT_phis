# Combination Template
#
#

__all__ = []
__author__ = ["Ramon Ángel Ruiz Fernández"]
__email__ = ["rruizferl@cern.ch"]


###########################
# OS + SS tagger branch 
# Combination 
###########################

OS_Kaon_Use            = 1
OS_Muon_Use            = 1
OS_Electron_Use        = 1
OS_Charm_Use           = 1
VtxCharge_Use          = 1


SS_Kaon_Use           = 1

SS_Kaon_TypeDec       = "Int_t"
SS_Kaon_BranchDec     = "B_SSKaonLatest_TAGDEC"
SS_Kaon_BranchDec     = "B_SSKaonLatest_TAGDEC"
SS_Kaon_TypeProb      = "Double_t"
SS_Kaon_BranchProb    = "B_SSKaonLatest_TAGETA"

OS_Kaon_BranchDec      = "B_OSKaonLatest_TAGDEC"
OS_Kaon_BranchProb     = "B_OSKaonLatest_TAGETA"
OS_Kaon_TypeDec        = "Int_t"
OS_Kaon_TypeProb       = "Double_t"

OS_Muon_BranchDec      = "B_OSMuonLatest_TAGDEC"
OS_Muon_BranchProb     = "B_OSMuonLatest_TAGETA"
OS_Muon_TypeDec        = "Int_t"
OS_Muon_TypeProb       = "Double_t"

OS_Electron_BranchDec  = "B_OSElectronLatest_TAGDEC"
OS_Electron_BranchProb = "B_OSElectronLatest_TAGETA"
OS_Electron_TypeDec    = "Int_t"
OS_Electron_TypeProb   = "Double_t"

OS_Charm_BranchDec     = "B_OSCharm_TAGDEC"
OS_Charm_BranchProb    = "B_OSCharm_TAGETA"
OS_Charm_TypeDec       = "Int_t"
OS_Charm_TypeProb      = "Double_t"

VtxCharge_BranchDec    = "B_OSVtxCh_TAGDEC"
VtxCharge_BranchProb   = "B_OSVtxCh_TAGETA"
VtxCharge_TypeDec      = "Int_t"
VtxCharge_TypeProb     = "Double_t"


# config current tuple {{{

RootFile = "${input_tuple}"
TupleName = "DecayTree"
Selection = ""
Nmax = -1

CalibrationMode = "Bs"
CalibrationLink = "MISTAG"
CalibrationModel = "POLY"
DoCalibrations = 0
CalibrationDegree = 1
UseNewtonRaphson = 0




# OS_Kaon_InOSComb = 1
# OS_Muon_InOSComb = 1
# OS_Electron_InOSComb = 1
# OS_Charm_InOSComb = 1
# VtxCharge_InOSComb = 1
# PerformOfflineCombination_OS = 1
# WriteCalibratedMistagBranches = 1
# OS_Combination_Write = 1

PerformOfflineCombination_OSplusSS = 1
OS_Kaon_InComb = 1
OS_Muon_InComb = 1
OS_Electron_InComb = 1
OS_Charm_InComb = 1
VtxCharge_InComb = 1
SS_Kaon_InComb = 1

WriteCalibratedMistagBranches = 1
Combination_Write = 1
CalibratedOutputFile = "${output_tuple}"

# }}}


# vim: fdm=marker ts=2 sw=2 sts=2 sr et
