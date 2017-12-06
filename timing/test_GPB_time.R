
op <- options(digits.secs = 6)


library(parallel)
library('GPB')
# library('poibin')
path_to_distributions = "~/storage/projects/bernmix/timing/data_for_timing/"
files = list.files(path_to_distributions)

gettime <-function(f)
{
  
  #print(f)
  d = data.matrix(read.table(paste(c(path_to_distributions, f),  sep = "", collapse = "")))
  probs = as.vector(d[1,])
  bval = round(as.vector(d[2,]))
  n = length(probs)
  
  aval = rep(0, n)
  weights = rep(1, n)
  
  kk = 0:sum(bval)
  # if (length(bval) != length(probs))
  #   print('aaaa')

  start_time <- Sys.time()
  pdf = dgpb(kk=kk, pp=probs, aval=aval, bval=bval, wts=weights)
  
  # kk = 0:length(probs)
  # pdf_poibin = dpoibin(kk=kk, pp=probs, wts=weights)
  
  end_time <- Sys.time()
  runtime = end_time - start_time
  
  write.table(runtime,
              file = paste(c('~/storage/projects/bernmix/timing/res_for_timing/time_',
                            f), sep = "", collapse = ""),
              col.names = F, row.names = F)
  return(runtime[[1]])
}


# Calculate the number of cores
no_cores <- 10

# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl,"path_to_distributions")
clusterEvalQ(cl, library('GPB'))

# x = sapply(files, gettime)

x = parSapply(cl, files, gettime)

# write.table(x=x, file = '~/storage/projects/bernmix/timing/res_for_timing/runtime_RES1.txt')
