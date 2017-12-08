
n_range = round(10 ^ seq(2, 5, 1/5))
r = c()
id = 0
for (i in 1:5)
{
  res = c()
  for(n in n_range)
  {
    id = id+1
    r = rbind(r, c(n, id))
    print(c(n, id))
    probs = runif(n)
    
    write.table(x = rbind(probs),
                file = paste(c('~/storage/projects/bernmix/timing/data_poibin/disrt_',
                               sprintf("%05i", id) ), sep = "", collapse = ""),
                col.names = F, row.names = F)
  }
}


write.table(x=r, file = '~/storage/projects/bernmix/timing/idx_poibin.txt')



