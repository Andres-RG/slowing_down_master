library(EWSmethods)

set.seed(125) #seed to ensure reproducible results
skylark_data <- data.frame(time = seq(1:50), abundance = rnorm(50,mean = 100,sd=20), trait = rnorm(50,mean=40,sd=5)) #dummy skylark dataset
#----
## ROLLING WINDOW
ews_metrics <- c("SD","ar1","skew") #the early warning signal metrics we wish to compute
roll_ews <- uniEWS(data = skylark_data[,1:2], metrics =  ews_metrics, method = "rolling", winsize = 50) #lets use a rolling window approach
roll_ews$EWS$cor #return the Kendall Tau correlations for each EWS metric
plot(roll_ews,  y_lab = "Skylark abundance")
#----
## EXPANDING WINDOW
exp_ews <- uniEWS(data = skylark_data[,1:2], metrics =  ews_metrics, method = "expanding", burn_in = 10, threshold = 2,  tail.direction = "one.tailed") #lets use a rolling window approach

head(exp_ews$EWS) #return the head of the EWS dataframe
#>   time metric.score metric.code rolling.mean rolling.sd threshold.crossed
#> 1   10    0.0000000         ar1    0.0000000         NA                 0
#> 2   11   -0.7071068         ar1   -0.3535534  0.5000000                 0
#> 3   12   -0.9223669         ar1   -0.5431579  0.4825449                 0
#> 4   13   -0.1500572         ar1   -0.4448827  0.4403012                 0
#> 5   14    0.8428688         ar1   -0.1873324  0.6906950                 0
#> 6   15   -1.8668306         ar1   -0.4672488  0.9229121                 0
#>   count.used        str
#> 1  108.91249         NA
#> 2   93.56813 -0.7071068
#> 3  109.56960 -0.7858522
#> 4  103.92341  0.6695997
#> 5  114.29655  1.4915428
#> 6   80.79726 -1.5164845
plot(exp_ews, y_lab = "Skylark abundance")
#----
trait_metrics <- c("SD", "ar1", "trait")
exp_ews_trait <- uniEWS(data = skylark_data[,1:2], metrics =  trait_metrics, trait = skylark_data$trait, method = "expanding", burn_in = 10, threshold = 2, tail.direction = "one.tailed")
plot(exp_ews_trait, y_lab = "Skylark abundance", trait_lab = "Body mass (g)", trait_scale = 5)
#----
##APROXIMACION MULTIVARIADA
set.seed(123)

octopus_spp_data <- matrix(nrow = 50, ncol = 5)
octopus_spp_data <- as.data.frame(cbind("time"=seq(1:50),sapply(1:dim(octopus_spp_data)[2], function(x){octopus_spp_data[,x] <- rnorm(50,mean=500,sd=200)}))) #create our hypothetical, uncollapsing ecosystem

oct_exp_ews <- multiEWS(data = octopus_spp_data, method = "expanding", threshold = 2, tail.direction = "one.tailed") #lets use an expanding window approach
plot(oct_exp_ews)

#----
##EWSNet
options(timeout = max(600, getOption("timeout"))) #due to possible internet issues, increase the timeout options from 60 seconds to 600
ewsnet_reset(remove_weights = FALSE, auto = TRUE)
#> Model weights downloaded
#> ewsnet_init(envname = "EWSNET_env", auto = TRUE) #prepares your workspace using 'reticulate' and asks to install Anaconda (if no appropriate Python found) and/or a Python environment before activating that environment with the necessary Python packages
#> Attention: may take up to 10 minutes to complete
#> EWSNET_env successfully found and activated. Necessary python packages installed

install.packages("reticulate")
library(reticulate)

reticulate::py_config() #confirm that "EWSNET_env" has been loaded
#> python:         /Users/ul20791/Library/r-miniconda-arm64/envs/EWSNET_env/bin/python
#> libpython:      /Users/ul20791/Library/r-miniconda-arm64/envs/EWSNET_env/lib/libpython3.10.dylib
#> pythonhome:     /Users/ul20791/Library/r-miniconda-arm64/envs/EWSNET_env:/Users/ul20791/Library/r-miniconda-arm64/envs/EWSNET_env
#> version:        3.10.14 | packaged by conda-forge | (main, Mar 20 2024, 12:51:49) [Clang 16.0.6 ]
#> numpy:          /Users/ul20791/Library/r-miniconda-arm64/envs/EWSNET_env/lib/python3.10/site-packages/numpy
#> numpy_version:  1.26.4
#> 
#> NOTE: Python version was forced by use_python() function

py_packages <- reticulate::py_list_packages() #list all packages currently loaded in to "EWSNET_env"
head(py_packages)
#>           package  version              requirement     channel
#> 1         absl-py    2.1.0            absl-py=2.1.0        pypi
#> 2       alabaster   0.7.16         alabaster=0.7.16        pypi
#> 3      astunparse    1.6.3         astunparse=1.6.3        pypi
#> 4           babel   2.15.0             babel=2.15.0        pypi
#> 5           bzip2    1.0.8              bzip2=1.0.8 conda-forge
#> 6 ca-certificates 2024.2.2 ca-certificates=2024.2.2 conda-forge

skylark_ewsnet <- ewsnet_predict(skylark_data$abundance, scaling = TRUE, ensemble = 25, envname = "EWSNET_env") #perform EWSNet assessment using white noise and all 25 models. The envname should match ewsnet_init()

skylark_ewsnet
#>                pred no_trans_prob smooth_trans_prob critical_trans_prob
#> 1 Smooth Transition    0.01444911         0.9638591          0.02169183