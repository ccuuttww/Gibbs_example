library(ProbBayes)

#beta binomial Gibbs
gibbs_betabin <- function(n, a, b, p = 0.5, iter){
  x <- matrix(0, iter, 2)
  for(k in 1:iter){
    y <- rbinom(1, size = n, prob = p)
    p <- rbeta(1, y + a, n - y + b )
    x[k, ] <- c(y,p)
  }
  x
}


a=b=5
n=20
k=13
vec<-gibbs_betabin(n, a, b, p = 0.5, iter = 100000)

#marginal distribution of y
hist(vec[,1][20000:100000],20,freq = FALSE)

#marginal distribution of p
hist(vec[,2][20000:100000],20,freq = FALSE)

#traceplot
plot(ts(vec[20000:100000,1]))



#marginal probability of Y = 13
length(which(vec[20000:100000,1]==k))
length(vec[20000:100000,1])
print(length(which(vec[20000:100000,1]==k))/length(vec[20000:100000,1]))



#marginal probability of Y =13  with exactly formula
#actually this is direct formula more information 
#please visit https://en.wikipedia.org/wiki/Beta-binomial_distribution
#the part of (while the marginal distribution m(k|μ, M) is given by)
part1 = gamma(n+1)/(gamma(k+1)*gamma(n-k+1))
part2 = gamma(a+b)/(gamma(a)*gamma(b))
part3 = (gamma(a+k)*gamma(n+b-k))/gamma(n+a+b)
print(part1*part2*part3)



###########################################################################


#gamma gamma poisson Gibbs
gibbs_Gamma_Gamma_poisson <- function(sum_data, n, b, iter){
x <- matrix(0, iter, 2)


for(k in 1:iter){    
      lambda <-rgamma(1,sum_data,n+b)  
      b <- rgamma(1,2,lambda+1)
      x[k, ] <- c(lambda,b)
}

x
}


n<-20
b<-2
iter <- 10000
set.seed(265)
sum_data<-sum(rpois(20,15))
vec <- gibbs_Gamma_Gamma_poisson(sum_data,n,b,iter)

#lambda marginal   
hist(vec[,1][4000:10000])
mean(vec[,1][4000:10000])
#b marginal
hist(vec[,2][4000:10000])
mean(vec[,2][4000:10000])


###########################################################################


#Bivariate normal distribuion Gibbs
gibbs<-function (n, rho)
   {
      mat <- matrix(ncol = 2, nrow = n)
      x <- 0
      y <- 0
      mat[1, ] <- c(x, y)
      for (i in 2:n) {
         x <- rnorm(1, rho * y, sqrt(1 - rho^2)) #
         y <- rnorm(1, rho * x, sqrt(1 - rho^2)) #
         mat[i, ] <- c(x, y)
         }
      mat
      }
bvn<-gibbs(10000,0.97)
mean(bvn[,1][4000:10000])
var(bvn[,1][4000:10000])
hist(bvn[,1][4000:10000],40)
mean(bvn[,2][4000:10000])
var(bvn[,2][4000:10000])
hist(bvn[,2][4000:10000],40)




# Gibbs Sampling Normal distribution with two unknown 
#
# Normal prior for mu
# mu ~ N(0,100)
#
# Gamma prior for precision 1/sigama^2
# 1/sigma^2 ~ Gamma(5,8)
#
# proecedure
# mu        |1/sigama^2,x1,......,xn
# 1/sigama^2|mu        ,x1,......,xn
# until converge

Generate_data_Gibbs_Normal <- function(n)
{
  #generate 53 data sample
  mu  <- 5
  var <- 2
  
  x <- rnorm(n,mu,sqrt(var))
  x
}


Normal_Gibbs_Nmu_Gphi <- function(p_a,p_b,p_mu,p_sigama_square,n,iter,data)
{
  sim_mu    <- rep(0,iter)
  sim_var   <- rep(1,iter)
  data_mean <- mean(data)
  for(i in 2:iter){
    
    # mu        |1/sigama^2,x1,......,xn
    part_a     <- (p_mu/p_sigama_square)  +  (n*data_mean/sim_var[i-1]) #p_sigama_square prior variance
    part_b     <- (1/p_sigama_square)     +  (n/sim_var[i-1])
    norm_mu    <- (part_a)/(part_b)
    
    part_c     <-  n/sim_var[i-1]
    part_d     <-  1/p_sigama_square
    norm_std   <-  1/sqrt(part_c + part_d)
    
    sim_mu[i]  <-  rnorm(1,norm_mu,norm_std)
    
    
    # 1/sigama^2|mu        ,x1,......,xn
    a          <- p_a + (n/2)
    b          <- 0.5*sum((data-sim_mu[i])^2)+p_b
    phi        <- rgamma(1,a,b)
    sim_var[i] <- 1/phi #become estimate of variance since phi = 1/sigama^2
    
  }
  output <- matrix(c(sim_mu,sim_var),iter,2)
  output
}


data  <-  Generate_data_Gibbs_Normal(53)
vec   <-  Normal_Gibbs_Nmu_Gphi(5,8,0,100,length(data),30000,data)

mean(vec[,1])
mean(vec[,2])


