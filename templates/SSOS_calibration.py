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


#Tagging Configuration
SS_Kaon_Use           = 1
SS_Kaon_TypeDec       = "Int_t"
SS_Kaon_BranchDec     = "B_SSKaonLatest_TAGDEC"
SS_Kaon_BranchDec     = "B_SSKaonLatest_TAGDEC"
SS_Kaon_TypeProb      = "Double_t"
SS_Kaon_BranchProb    = "B_SSKaonLatest_TAGETA"


OS_Kaon_Use = 1
OS_Muon_Use = 1
OS_Electron_Use = 1
OS_Charm_Use = 1
VtxCharge_Use = 1

OS_Kaon_BranchDec = "B_OSKaonLatest_TAGDEC"
OS_Kaon_BranchProb = "B_OSKaonLatest_TAGETA"
OS_Kaon_TypeDec = "Int_t"
OS_Kaon_TypeProb = "Double_t"

OS_Muon_BranchDec = "B_OSMuonLatest_TAGDEC"
OS_Muon_BranchProb = "B_OSMuonLatest_TAGETA"
OS_Muon_TypeDec = "Int_t"
OS_Muon_TypeProb = "Double_t"

OS_Electron_BranchDec = "B_OSElectronLatest_TAGDEC"
OS_Electron_BranchProb = "B_OSElectronLatest_TAGETA"
OS_Electron_TypeDec = "Int_t"
OS_Electron_TypeProb = "Double_t"

OS_Charm_BranchDec = "B_OSCharm_TAGDEC"
OS_Charm_BranchProb = "B_OSCharm_TAGETA"
OS_Charm_TypeDec = "Int_t"
OS_Charm_TypeProb = "Double_t"

VtxCharge_BranchDec = "B_OSVtxCh_TAGDEC"
VtxCharge_BranchProb = "B_OSVtxCh_TAGETA"
VtxCharge_TypeDec = "Int_t"
VtxCharge_TypeProb = "Double_t"
#Statistics Input
UseNewtonRaphson = 0 #0 is Minuit

