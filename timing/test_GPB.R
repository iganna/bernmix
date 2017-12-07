library('GPB')

n_range = c(50, 100, 500, 1000, 5000, 10000)
m_range = round(10 ^ seq(3, 7, 1/5))
op <- options(digits.secs = 6)
runtime_mx = matrix(, nrow = 0, ncol = length(n_range) * length(m_range))
id = 0
for (i in 1:20)
{
  res = c()
  for(n in n_range)
    for (m in m_range)
    {
      if(n >= m)
        next
      
      id = id+1
      print(c(n, m, i))
      probs = runif(n)
      
      a = rep(0, n)
      b = runif(n)
      d = m/(sum(b) - sum(a))
      aval = round(a * d)
      bval = round(b * d)
      if(min(bval) == 0)
      {
        print('aaaa')
        bval[bval == min(bval)] = 1
      }
      while(sum(bval) != m)
      {
        idx = 1:n
        idx = idx[bval == max(bval)]
        bval[idx[1]] = bval[idx[1]] - 1 * sign(sum(bval) - m)
      }
      
      
      print(sum(bval))
      weights = rep(1, n)
      
      write.table(x = rbind(probs, bval), 
                  file = paste(c('~/storage/projects/bernmix/timing/data_for_timing/disrt_', 
                                 sprintf("%05i", id) ), sep = "", collapse = ""),
                  col.names = F, row.names = F)
      

    }

}

#write.table(x=runtime_mx, file = 'runtime_mx.txt')



