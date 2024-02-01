# template for SS Dspi tagging calibration
#
#

__all__ = []
__author__ = ["Ramón Ángel Ruiz Fernández"]
__email__ = ["rruizfer@cern.ch"]


#Sample Configuration
datafile  = "${input_tuple}" 
TupleName = "DecayTree"
#In principle no needed selection for epm -> Chinesse configuration below
# Selection = "tagss_eta_new>0&&tagss_eta_new<=0.5&&B_LOKI_DTF_CTAU_CONV>=0.3&&B_LOKI_DTF_CTAU_CONV<=10&&tagss_dec_new!=0"
Selection = ""
Nmax        = -1
CalibrationMode = "Bs"
DoCalibrations  = 1
BranchID                  = "${idvar}"
UseWeight                 =  1
WeightFormula             = "${sweight}"


#Time Configuration 
UseTau                    = 1
UseTauErr                 = 1
TypeTau                   = "Double_t"
BranchTau                 = "${timevar}"
BranchTauErr              = "${timeerrvar}"
TauUnits  = "ps"
DeltaM = ${DeltaM}
DeltaMErr = ${DeltaMerr}
Lifetime = ${Lifetime}
LifetimeErr = ${LifetimeErr}
DeltaGamma = ${DeltaGamma}
DeltaGammaErr = ${DeltaGammaErr}
ResolutionGaussian1_Frac = 1
ResolutionGaussian1_A = ${p0}
ResolutionGaussian1_B = ${p1}
ResolutionGaussian1_C = 0
ResolutionGaussian2_Frac = 0

# DeltaM = 17.765
# DeltaMErr = 0.006
# Lifetime = 1.516
# LifetimeErr = 0.006
# DeltaGamma = 0.082
# DeltaGammaErr = 0.005
# ResolutionGaussian1_A =  0.0028
# ResolutionGaussian1_B = 1.1311

#Tagging Configuration
SS_Kaon_Use           = 1
SS_Kaon_TypeDec       = "Int_t"
SS_Kaon_BranchDec     = "B_SSKaonLatest_TAGDEC"
SS_Kaon_BranchDec     = "B_SSKaonLatest_TAGDEC"
SS_Kaon_TypeProb      = "Double_t"
SS_Kaon_BranchProb    = "B_SSKaonLatest_TAGETA"

#Statistics Input
# UseNewtonRaphson = 0


