library('GPB')

n_range = round(10 ^ seq(1, 2.8, 0.2))
m_range = round(10 ^ seq(3, 6, 0.2))
op <- options(digits.secs = 6)
runtime_mx = matrix(, nrow = 0, ncol = length(n_range) * length(m_range))
for (i in 1:20)
{
  res = c()
  for(n in n_range)
    for (m in m_range)
    {
      print(c(n, m, i))
      probs = runif(n)
      
      a = rep(0, n)
      b = runif(n)
      d = m/(sum(b) - sum(a))
      aval = round(a * d)
      bval = round(b * d)
      bval[bval == max(bval)] = bval[bval == max(bval)] + (m-sum(bval))
      print(sum(bval))
      weights = rep(1, n)
      
      kk = 0:sum(bval)
      start_time <- Sys.time()
      pdf = dgpb(kk=kk, pp=probs, aval=aval, bval=bval, wts=weights)
      end_time <- Sys.time()
      runtime = end_time - start_time
      res = c(res, runtime[[1]])
    }
  runtime_mx = rbind(runtime_mx, res)
}

write.table(x=runtime_mx, file = 'runtime_mx.txt')



