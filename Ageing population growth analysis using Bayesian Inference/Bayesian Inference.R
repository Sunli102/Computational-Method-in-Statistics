########################
### read in data set ###
########################

dat=read.csv('/Users/lisun/GSU/Computational/final project/datasets/Population_age65+/API_SP.POP.65UP.TO.ZS_DS2_en_csv_v2_10517738.csv', 
             header = F,row.names = 1 ) 


#### clean data
# str(dat)          # 267 obs. of  63 variables
data1 = dat[-c(1:2), -c(1:3,63)] 
colnames(data1) = as.character(unlist(data1[1,]))   # name column names
data1 = data1[-1, ]             # remove column name row 
data1 = data1[, -59]            # remove last empty column year 2018

# remove the variables with missing values 
id_no_miss = which(apply(is.na(data1),1,sum)==0)
data= data1[id_no_miss,]
# str(data)       # 237 obs. of  58 variables



### likelihood function 
lik_log = function(th, pre, pos)
{ a1=th[1] ; b1=th[2] ; a2=th[3] ;  b2=th[4] 
  sum((a1-1)*log(pre)-(1/b1)*pre - a1*log(b1) - log(gamma(a1))) +
  sum((a2-1)*log(pos)-(1/b2)*pre - a2*log(b2) - log(gamma(a2)))}

### prior function 
prior_log = function(th, th0)
{
  a1=th[1] ; b1=th[2] ; a2=th[3] ;  b2=th[4]
  a01 = th0[1]; b01 = th0[2]; a02 = th0[3]; b02 = th0[4]
  if (a1 <= 0 | a2 <= 0 | b1 <= 0 | b2 <= 0 )
  {
    return(0)
  }else
  {  result =
    -log(a01*sqrt(2*pi)) -((a1-a01)^2/(2*a01^2))
    -log(b01*sqrt(2*pi)) -((b1-b01)^2/(2*b01^2))
    -log(a02*sqrt(2*pi)) -((a2-a02)^2/(2*a02^2))
    -log(b02*sqrt(2*pi)) -((b2-b02)^2/(2*b02^2))
  return(result)
  }
}

### posterior function 
post_log = function(th, th0,pre, pos){prior_log(th, th0) + lik_log(th, pre, pos)}


# Here is what does the MCMC (Metropolis method): 
my.MCMC = function(nit, th0, pre, pos)
{
  th = th0
  results_log = matrix(0, nrow=nit , ncol=4) 
  accepted_log = matrix(th0, nrow=1 , ncol=4) 
  rejected_log = matrix(th0, nrow=1 , ncol=4) 
  results_log[1,] = th0 
  
  initial.time = proc.time()
  for(it in 2:nit)
  {
    cand = th + rnorm(4,  mean = 0, sd=0.1) 
    res = post_log(cand,th0, pre,pos) - post_log(th,th0, pre,pos)
    
    if (runif(1) < exp(res))
    {
      th=cand
      accepted_log = rbind(accepted_log, cand)
    }else
    {
      rejected_log = rbind(rejected_log, cand)
    }
    results_log[it,] = th
  }
  
  time =  proc.time() - initial.time; time 
  my_list = list('results_log'= results_log, 'accepted_log' = accepted_log,
       'rejected_log'=rejected_log, 'time' = time)
  return(my_list)
}


# Starting values 
a01 = 10; b01 = 3; a02 = 11; b02 = 5
th0=c(a01, b01 ,a02, b02)
nit=50000

###################################
### compare the first 2 years #####
###   year 1960 - year 1961   #####
###################################
pre.1 = data[,1]
pos.1 = data[,2]

### draw a picture 
xlim =c(min(pre.1),max(pre.1)) 
par(mfrow=c(2 ,1)) 
hist (pre.1,breaks = 100 , col='red',xlim=xlim, 
      main = 'Histogram of year 1960', xlab = 'ratio')
hist (pos.1,breaks = 100 , col='red',xlim=xlim,
      main = 'Histogram of year 1961', xlab = 'ratio')


# Call MCMC function
result.1 = my.MCMC(nit = nit, th0 = th0, pre = pre.1, pos = pos.1)

results_log.1 = result.1$results_log
time.1 = result.1$time; time.1        # 22.004


# Take a peek at what we got edit (results) 
par(mfrow=c(4 ,1)) 
plot (results_log.1 [ ,1], col = 'red') ;
plot(results_log.1 [ ,2], col ='red') ;
plot (results_log.1 [ ,3], col = 'red') ;
plot(results_log.1 [ ,4], col = 'red') 


# omit burn-in data
res=results_log.1[(0.25*nit):nit,] 
a1s.1 = res[,1];b1s.1 = res[,2] ;a2s.1 = res[,3];b2s.1 = res[,4] 

# plot and histogram 
par(mfrow=c(4 ,2)) 
plot (a1s.1, col = rgb(1,0,0)) ;hist(a1s.1, col = rgb(1,0,0))
plot (b1s.1, col = rgb(0.8,0,0)) ;hist(b1s.1, col = rgb(0.8,0,0))
plot (a2s.1, col = rgb(0.6,0,0)) ;hist(a2s.1, col = rgb(0.6,0,0))
plot (b2s.1, col = rgb(0.4,0,0)) ;hist(b2s.1, col = rgb(0.4,0,0)) 



# calculate the simulated estimators
a1.1 = mean(a1s.1);a1.1    # 4.472469  4.430022
b1.1 = mean(b1s.1);b1.1    # 1.06927   1.079413
a2.1 = mean(a2s.1);a2.1    # 4.516944  4.48286
b2.1 = mean(b2s.1);b2.1    # 1.057707  1.066407
a01.1 = mean(pre.1)^2/var(pre.1);a01.1  # 3.533503
b01.1 = var(pre.1)/mean(pre.1);b01.1    # 1.340233
a02.1 = mean(pos.1)^2/var(pos.1);a02.1  # 3.487251
b02.1 = var(pos.1)/mean(pos.1);b02.1    # 1.36551
# > 3.533503*1.340233 -3.487251*1.36551
# [1] -0.02615879

# compare the mean difference between the first 2 years
par(mfrow=c(2 ,1)) 
mu1 = a1s.1*b1s.1; mu2 = a2s.1*b2s.1
plot (mu1-mu2,col = 'red') 
hist (mu1-mu2, main = 'Histogram of the mean difference between 1960-1961', col = 'red') 
mean_dif = mean(mu1-mu2)
abline(v = mean_dif, col = 'blue',lwd = 2)
mean((mu1-mu2) < 0)   #0.4892604 0.4869203  0.5153196





###################################
### compare the last 2 years ######
###   year 2016 - year 2017  ######
###################################

pre.2 = data[,57]
pos.2 = data[,58]

### draw a picture 
xlim =c(min(data),max(data)) 
par(mfrow=c(2 ,1)) 
hist (pre.2,breaks = 100 , col='green',xlim=xlim,
      main = 'Histogram of year 2016', xlab = 'ratio')
hist (pos.2,breaks = 100 , col='green',xlim=xlim,
      main = 'Histogram of year 2017', xlab = 'ratio')


# Call MCMC function
result.2 = my.MCMC(nit = nit, th0 = th0, pre = pre.2, pos = pos.2)

results_log.2 = result.2$results_log
time.2 = result.2$time; time.2     # 31.389



# Take a peek at what we got edit (results) 
par(mfrow=c(4 ,1)) 
plot(results_log.2[,1], col = 'green');plot(results_log.2[,2],col = 'green');
plot(results_log.2[,3], col = 'green');plot(results_log.2[,4],col = 'green') 


# omit burn-in data
res=results_log.2 [(0.25*nit):nit,] 
a1s.2 = res[,1];b1s.2 = res[,2] ;a2s.2 = res[,3];b2s.2 = res[,4] 

# pliot final data set and its histogram
par(mfrow=c(4 ,2)) 
plot(a1s.2, col = rgb(0,1,0));hist(a1s.2, col = rgb(0,1,0))
plot(b1s.2, col = rgb(0,0.8,0));hist(b1s.2, col = rgb(0,0.8,0))
plot(a2s.2, col = rgb(0,0.6,0));hist(a2s.2, col = rgb(0,0.6,0))
plot(b2s.2, col = rgb(0,0.4,0)) ;hist(b2s.2, col = rgb(0,0.4,0)) 




# calculate the simulated estimators
a1.2 = mean(a1s.2);a1.2    # 2.279482    2.281352
b1.2 = mean(b1s.2);b1.2    # 3.747456    3.745995
a2.2 = mean(a2s.2);a2.2    # 2.4489      2.442178
b2.2 = mean(b2s.2);b2.2    # 3.499152    3.509284
a01.2 = mean(pre.2)^2/var(pre.2);a01.2  # 2.12209
b01.2 = var(pre.2)/mean(pre.2);b01.2    # 3.989415
a02.2 = mean(pos.2)^2/var(pos.2);a02.2  # 2.129495
b02.2 = var(pos.2)/mean(pos.2);b02.2    # 4.063988
# # > 2.12209*3.989415-2.129495*4.063988
# # [1] -0.1883444


# compare the mean difference between the last 2 years
par(mfrow=c(2 ,1)) 
mu1 = a1s.2*b1s.2; mu2 = a2s.2*b2s.2
plot (mu1-mu2,col = 'green') 
hist (mu1-mu2, main = 'Histogram of the mean difference between 2016-2017', col = 'green') 
mean_dif = mean(mu1-mu2)  
abline(v = mean_dif, col = 'blue',lwd = 2)
mean((mu1-mu2) < 0)   # 0.5225391  0.5220927 0.4999333



# plot and Histogram of the mean difference between 1960-2017
par(mfrow=c(2 ,1)) 
mu1 = a1s.1*b1s.1; mu4 = a2s.2*b2s.2
plot (mu1-mu4) 
hist (mu1-mu4, main = 'Histogram of the mean difference between 1960-2017') 
mean_dif = mean(mu1-mu4)
abline(v = mean_dif, col = 'blue',lwd = 2)
mean((mu1-mu4) < 0)


#########################################
###### data analysis for year 1960 ######
#########################################
### algorithm for sampling 
accepted = result.1$accepted_log
rejected = result.1$rejected_log

nrow(accepted)
#[1] 7598
nrow(rejected)
#[1] 42403

### calculate the simulated estimators
a1.1 = mean(a1s.1);a1.1
# [1] 4.39136
b1.1 = mean(b1s.1);b1.1   
# [1] 1.08786


### plot data sampling
plot(rejected[,1], rejected[,2], type = 'p', col = 'red', 
     main = 'MCMC sampling for a and b (all samples)', xlab = 'a', ylab = 'b' )
lines(accepted[,1],accepted[,2], type = 'l', col = 'blue')

plot(rejected[1:100,1], rejected[1:100,2], type = 'p', col = 'red',
     main = 'MCMC sampling for a and b (first 100 samples)', xlab = 'a', ylab = 'b')
lines(accepted[1:100,1],accepted[1:100,2], type = 'l', col = 'blue')


### the traces of a and b and the histogram of the traces
par(mfrow=c(2 ,2))
plot (a1s.1, col = rgb(1,0,0), main = 'trace for a',xlab = 'iteration') 
hist(a1s.1, col = rgb(1,0,0), main = 'histogram of a', xlab = 'a')
plot (b1s.1, col = rgb(0.8,0,0),main = 'trace for b',xlab = 'iteration') 
hist(b1s.1, col = rgb(0.8,0,0),main = 'histogram of b', xlab = 'b')

### Histogram and gamma funciton curve
xlim = c(min(pre.1), max(pre.1))
hist (pre.1,col= rgb(1,0, 0),xlim=xlim, breaks = 100, prob = T, 
      main = 'Histogram and simulated gamma function curve for year 1960', 
      xlab = 'ratio')

x = rgamma(nit, shape = a1.1, scale = b1.1)
d = density(x)
lines(d$x,d$y,col = 'blue', type = 'l', lwd = 3)








xlim = c(min(pre.1), max(pre.1))
hist (pre.1,col= rgb(0,0,0.5, 0.2),xlim=xlim, breaks = 100, prob = T)

x = rgamma(nit, shape = a1.1, scale = b1.1)
d = density(x)
lines(d$x,d$y,col = 'blue', type = 'l')

x1 = rgamma(nit, shape = a01.1, scale = b01.1)
d1 = density(x1)
lines(d1$x,d1$y,col = 'green', type = 'l')




xlim = c(min(pos.1), max(pro.1))
hist (pos.1,col= rgb(0.5,0,0, 0.2),xlim=xlim, breaks = 100, prob = T)

x = rgamma(nit, shape = a2.1, scale = b2.1)
d = density(x)
lines(d$x,d$y,col = 'blue', type = 'l')

x1 = rgamma(nit, shape = a02.1, scale = b02.1)
d1 = density(x1)
lines(d1$x,d1$y,col = 'green', type = 'l')




xlim = c(min(pre.2), max(pre.2))
hist (pre.2,col= rgb(0,0.5,0, 0.2),xlim=xlim, breaks = 100, prob = T)

x = rgamma(nit, shape = a1.2, scale = b1.2)
d = density(x)
lines(d$x,d$y,col = 'blue', type = 'l')

x1 = rgamma(nit, shape = a01.2, scale = b01.2)
d1 = density(x1)
lines(d1$x,d1$y,col = 'green', type = 'l')






xlim = c(min(pos.2), max(pos.2))
hist (pos.2,col= rgb(0.5,0.5,0.5, 0.2),xlim=xlim, breaks = 100, prob = T)

x = rgamma(nit, shape = a2.2, scale = b2.2)
d = density(x)
lines(d$x,d$y,col = 'blue', type = 'l')

x1 = rgamma(nit, shape = a02.2, scale = b02.2)
d1 = density(x1)
lines(d1$x,d1$y,col = 'green', type = 'l')



 

