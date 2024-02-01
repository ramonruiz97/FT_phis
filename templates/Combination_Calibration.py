# template for Bu2JpsiKplus tagging calibration
#
#

__all__ = []
__author__ = ["Ramón Ángel Ruiz Fernández"]
__email__ = ["rruizfer@cern.ch"]


# OS tagger {{{


# Input configuration {{{
# This is the file/directory that you want to run
# Multiple files can be specified by setting NumFiles = N
# and then setting RootFile_1, RootFile_2, ..., RootFile_N


# OS_Combination_Use = 1
# OS_Combination_TypeDec = "Int_t"
# OS_Combination_BranchDec = "OS_Combination_DEC"
# OS_Combination_TypeProb = "Double_t"
# OS_Combination_BranchProb = "OS_Combination_ETA"


Combination_Use = 1
Combination_TypeDec = "Int_t"
Combination_BranchDec = "Combination_DEC"
Combination_TypeProb = "Double_t"
Combination_BranchProb = "Combination_ETA"

NumFiles = 1
RootFile_1 = "${input_tuple}"

TupleName = "DecayTree"
Selection = ""
Nmax = -1  # Events to run, -1 means all

# }}}

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

# calibration settings {{{

CalibrationMode = "Bs"
DoCalibrations = 1
CalibrationLink = "MISTAG"
CalibrationDegree = 1
CalibrationModel = "POLY"
UseNewtonRaphson = 0

# }}}

# branch names and types {{{

BranchID = "${idvar}"
UseWeight = 1
WeightFormula = "${sweight}"

# }}}


# vim: fdm=marker ts=2 sw=2 sts=2 sr et
