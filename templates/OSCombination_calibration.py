
__all__ = []
__author__ = ["Ramón Ángel Ruiz Fernández"]
__email__ = ["rruizfer@cern.ch"]


#Sample Configuration
datafile  = "${input_tuple}" 
TupleName = "DecayTree"
Selection = ""
Nmax        = -1
CalibrationMode = "Bs"
DoCalibrations  = 1
BranchID                  = "${idvar}"
UseWeight                 =  1
WeightFormula             = "${sweight}"


OS_Combination_Use = 1
OS_Combination_TypeDec = "Int_t"
OS_Combination_BranchDec = "OS_Combination_DEC"
OS_Combination_TypeProb = "Double_t"
OS_Combination_BranchProb = "OS_Combination_ETA"
#Statistics Input
UseNewtonRaphson = 0 #0 is Minuit

