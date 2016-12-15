# Script for STAT240 Project

# Two sample permutation test in mean when underlying distribution are different
# Estimate rejection probability by simulation
traditional_perm_test_mean<-function(m1,m2,n){
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
  return(p_val)
}

# Chung et al's method

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

Chung_robust_test_median<-function(m1,m2,n){
  p_val = c()
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
  return (p_val)
}



# PB test by Krishnamoorthy et al

Krishnamoorthy <- function(m1,m2,n){
  p_val = c()
  for (i in 1:n){
    M = 1000
    x = rnorm(m1,mean = 0, sd = 1)
    y = rnorm(m2,mean = 0, sd = sqrt(5))
    S_x = var(x)
    S_y = var(y)
    T_N0 = m1*mean(x)**2/S_x + m2*mean(y)**2/S_y -(m1*mean(x)/S_x+m2*mean(y)/S_y)**2/(m1/S_x+m2/S_y)
    count = 0
    for (j in 1:M){
      Z1 = rnorm(1,mean = 0, sd = 1)
      Z2 = rnorm(1,mean = 0, sd = 1)
      chi_1 =  rchisq(1,m1-1)
      chi_2 =  rchisq(1,m2-1)
      T_NB = Z1**2*(m1-1)/chi_1 + Z2**2*(m2-1)/chi_2 - ((sqrt(m1)*Z1*(m1-1)/(sqrt(S_x)*chi_1))+(sqrt(m2)*Z2*(m2-1)/(sqrt(S_y)*chi_2)))**2/((m1*(m1-1))/(S_x*chi_1)+(m2*(m2-1))/(S_y*chi_2))
      if (T_NB>T_N0){count = count+1}
    }
    p_val = c(p_val, count/M)
  }
  return (p_val)
} 


# Lix et al (2005) robust method by trimmed statistics
trimmed_mean <- function(x,gamma){
  g = floor(length(x)*gamma)
  x = x[order(x)]
  result = sum(x[g+1:length(x)-g])/(length(x)-2*g)
  return (result)
}

winsoried_var  <- function(x,gamma){
  g = floor(length(x)*gamma)
  x = x[order(x)]
  x[1:g] = rep(x[g+1],g)
  x[length(x)-g+1:length(x)] = rep(x[length(x)-g],g)
  return (var(x))
}


Lix <- function(m1,m2,n,gamma){
  p_val = c()
  for (i in 1:n){
    x = rnorm(m1,mean = 0, sd = 1)
    y = rt(m2,5)
    x_trim = trimmed_mean(x,gamma)
    y_trim = trimmed_mean(y,gamma)
    x_std = sqrt(winsoried_var(x,gamma))
    y_std = sqrt(winsoried_var(y,gamma))
    h1 = m1-2*floor(gamma*m1)
    h2 = m2-2*floor(gamma*m2)
    test_stat = (x_trim - y_trim)**2*((m1-1)*x_std/(h1*(h1-1))+(m2-1)*y_std/(h2*(h2-1)))**(-1)
    A1 = x_std/m1
    A2 = y_std/m2
    A = 1+1/2*(1/(m1-1)*(A1/(A1+A2))**2+1/(m2-1)*(A2/(A1+A2))**2)
    B = 3/4*(A**(-2)*A1**2/(m1-1)+A**(-2)*A2**2/(m2-1))
    p_val = c(p_val,pchisq(test_stat,1,A+B))
  } 
  return (p_val)
} 

p_val_1 = Lix(50,50,10000,0.2)
p_val_2 = Lix(50,100,10000,0.2)
p_val_3 = Lix(300,300,10000,0.2)
p_val_4 = Lix(500,500,10000,0.2)
p_val_5 = Lix(1000,1000,10000,0.2)
p_val_6 = Lix(10000,50000,10000,0.2)

# Histogram of P-values 
par(mfrow=c(2,3))
hist(p_val_1, main="m=50,n=50",breaks=20)
hist(p_val_2, main="m=50,n=100",breaks=20)
hist(p_val_3, main="m=300,n=300",breaks=20)
hist(p_val_4, main="m=500,n=500",breaks=20)
hist(p_val_5, main="m=1000,n=1000",breaks=20)
hist(p_val_6, main="m=1000,n=5000",breaks=20)
