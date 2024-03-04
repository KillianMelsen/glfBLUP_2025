# !!! IMPORTANT: SET MKL_NUMTHREADS=1 and MKL_DYNAMIC=TRUE in ~/.profile !!!
# Total runtime of 56 minutes for the 250 datasets.
Sys.getenv("MKL_NUM_THREADS")
Sys.getenv("MKL_DYNAMIC")

# ### Part 1 (56 minutes):
# 
# # wd <- "C:/Users/Killian/Desktop/gfblup-methodological-paper"
# wd <- "~/gfblup_methodology"
# start <- Sys.time()
# 
# setwd(wd)
# rm(list = setdiff(ls(), "start"))
# tictoc::tic.clearlog()
# source("analyses_scripts/hyper/univariate_hyper.R")
# 
# setwd(wd)
# rm(list = setdiff(ls(), "start"))
# tictoc::tic.clearlog()
# source("analyses_scripts/hyper/lsBLUP_hyper.R")
# 
# setwd(wd)
# rm(list = setdiff(ls(), "start"))
# tictoc::tic.clearlog()
# source("analyses_scripts/hyper/siBLUP_hyper.R")
# 
# setwd(wd)
# rm(list = setdiff(ls(), "start"))
# tictoc::tic.clearlog()
# source("analyses_scripts/hyper/gfBLUP_hyper.R")
# 
# end <- Sys.time()
# cat("\n\nDone!")
# end - start



### Part 2 ( minutes):

# wd <- "C:/Users/Killian/Desktop/gfblup-methodological-paper"
wd <- "~/gfblup_methodology"
start <- Sys.time()

setwd(wd)
rm(list = setdiff(ls(), "start"))
tictoc::tic.clearlog()
source("analyses_scripts/hyper/lsBLUP_hyper_RF075.R")

setwd(wd)
rm(list = setdiff(ls(), "start"))
tictoc::tic.clearlog()
source("analyses_scripts/hyper/siBLUP_hyper_RF075.R")

end <- Sys.time()
cat("\n\nDone!")
end - start
