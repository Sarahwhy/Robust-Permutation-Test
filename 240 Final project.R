# Script for STAT240 Project

# Two sample permutation test in mean when underlying distribution are different
# Estimate rejection probability by simulation
n = 5000
m1 = 500
m2 = 500
p_val = c()

for (i in 1:n){
  x = rnorm(m1,mean = 0, sd = 1)
  y = rnorm(m2,mean = 0, sd = 5)
  perm = 999
  ts = median(x) - median(y)
  reps <- numeric(perm)
  for (j in 1:perm){
    all = c(x,y)
    k = 1:length(all)
    p = sample(k, size=length(x), replace=FALSE)
    permsample1 = all[p]
    permsample2 = all[-p]
    reps[j] = median(permsample1) - median(permsample2)
  }
  pvalue = mean(abs(reps)>=ts)
  p_val = c(p_val,pvalue)
}
# Probability of wrong prediction
hist(p_val)
length(p_val[p_val<0.05])/length(p_val)


# Subsampling method



# Chung et al's method
n = 10000
m1 = 51
m2 = 101
p_val = c()

generate_prob <-function(m){
  result = c()
  for (i in 1:m){
    result = c(result,pbinom((m-1)/2, size=m, prob=(i-1)/m)-pbinom((m-1)/2, size=m, prob=i/m))
  }
  return(result)
}


test_stat <-function(x,y){
  p = length(x)/(length(x)+length(y))
  x_sorted = sort(x)
  y_sorted = sort(y)
  median_x = median(x)
  median_y = median(y)
  temp_x = (x_sorted-median_x)**2
  temp_y = (y_sorted-median_y)**2
  prob_x = generate_prob(length(x))
  prob_y = generate_prob(length(y))
  var = sqrt(sum(prob_x*temp_x)*length(x)/p+sum(prob_y*temp_y)*length(y)/(1-p))
  test_stat_result = (length(x)+length(y))**0.5*(median_x-median_y)/var
  return(test_stat_result)
}

for (i in 1:n){
  x = rnorm(m1,mean = 0, sd = 1)
  y = rnorm(m2,mean = 0, sd = sqrt(5))
  perm = 999
  ts = test_stat(x,y)
  reps <- numeric(perm)
  for (j in 1:perm){
    all = c(x,y)
    k = 1:length(all)
    p = sample(k, size=length(x), replace=FALSE)
    permsample1 = all[p]
    permsample2 = all[-p]
    reps[j] = test_stat(permsample1,permsample2)
  }
  pvalue = mean(reps>=ts)
  p_val = c(p_val,pvalue)
}

# Probability of wrong prediction
hist(p_val, main="Histogram of P-value for Studentized Median Test P~N(0,1) Q~N(0,5)",breaks=20)
length(p_val[p_val<0.05])/length(p_val)


