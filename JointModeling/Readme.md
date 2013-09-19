## Joint Modeling Example Code ##

The code in `Example_JointModeling_Code.R` fits a joint model for longitudinal and survival data using the Monte Carlo ECM (MC-ECM) algorithm. 

Details of the type of model and algorithm being used can be found in:

+ Hsieh, F., Tseng, Y.K. and Wang, J.-L. (2006) Joint modeling of survival and longitudinal data: likelihood approach revisited. *Biometrics*, 62-4, pp. 1037--1043.

Available at e.g., <www.jstor.org/stable/4124524‎>, <www.stat.ncu.edu.tw/teacher/Tsengyk/biom_570%5B001-007%5D.pdf>.

The file `Example_JointModeling_Code.R` is entirely self-contained and can be run via `source('Example_JointModeling_Code.R')` (from the appropriate working directory). The file loads in the data from `Aids.RData`, and will output a file `Results.RData` containing the results of the fit.



