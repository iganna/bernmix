
op <- options(digits.secs = 6)


library(parallel)
library('poibin')
path_to_distributions = "~/storage/projects/bernmix/timing/data_poibin/"
files = list.files(path_to_distributions)

gettime <-function(f)
{
  
  #print(f)
  d = data.matrix(read.table(paste(c(path_to_distributions, f),  sep = "", collapse = "")))
  probs = as.vector(d[1,])
  n = length(probs)
  
  weights = rep(1, n)
  kk = 0:n
  
  start_time <- Sys.time()

  pdf_poibin = dpoibin(kk=kk, pp=probs, wts=weights)
  
  end_time <- Sys.time()
  runtime = end_time - start_time
  t = as.numeric(runtime, units = "secs")
  write.table(t,
              file = paste(c('~/storage/projects/bernmix/timing/res_poibin_R/time_',
                             f, '.txt'), sep = "", collapse = ""),
              col.names = F, row.names = F)
  return(runtime[[1]])
}


# Calculate the number of cores
no_cores <- 10

# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl,"path_to_distributions")
clusterEvalQ(cl, library('poibin'))

# x = sapply(files, gettime)

x = parSapply(cl, files, gettime)

# write.table(x=x, file = '~/storage/projects/bernmix/timing/res_for_timing/runtime_RES1.txt')
