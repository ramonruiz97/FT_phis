### Tagging for Phi_s

You need your eos credentials to download from eos needed tuples:

```
kinit {user}@CERN.CH
```


You also need to have EspressoPerformanceMonitor in your work directory.

```
source installer.sh
```


For running OS_Cal and SS_Cal from DsPi you have to call the following rule: 
```
snakemake output/Combination_Plan/201{5,6,7,8}/Bs2JpsiPhi/v1r0 -j -F
```

For running OS_Cal from Bu and SS_Cal from DsPi you have to run:
````
snakemake output/Combination_baseline/201{5,6,7,8}/Bs2JpsiPhi/v1r0 -j -F
```` 



