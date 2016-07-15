
#Integration via monte carlo simulation:

lambda = 0.5

a=0.25
b=0.5
n_samples=10000000 # 10 million samples

proposal_values <- runif(min=a,max=b,n=n_samples) #obtaining proposal sampled values
hist(proposal_values)

z=proposal_values

samples <- (1/lambda)*exp((-z)/lambda)  #obtaining samples from the negative exponetial distribution

hist(samples)

integral_resolution <-  ( (b-a) /n_samples ) * sum(samples) #Obtaining the integral solution

print(round(integral_resolution,3))

