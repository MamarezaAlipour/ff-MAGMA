library( ff-MAGMA)
options( echo=FALSE)
test.for.zero.flag<- 1

#
##### generate test data
#

genCovMat = function(x, theta, lambda) {
  distanceMatrix<- rdist(x,x)
  Sigma<- exp( -distanceMatrix/theta ) + diag(x=lambda, nrow=nrow(distanceMatrix))
  return(Sigma)
}

#generate observation locations
n=500
x = matrix(runif(2*n), nrow=n)

#generate observations at the locations
trueTheta = .2
trueLambda = .1
Sigma = genCovMat(x, trueTheta, trueLambda)

U = magmaChol(Sigma, nGPUs=1)
y = t(U)%*%as.vector(rnorm(n))

#
######set MLE computation parameters
#

testThetas = seq(from=trueTheta/2, to=2*trueTheta, length=20)
par.grid=list(theta=testThetas)
guessLambda = trueLambda

#
##### test using distance matrix
#

print("testing using distance matrix")

set.seed(1)
out1 = mKrig.MLE(x, y, lambda=guessLambda, par.grid=par.grid)
lambda.MLE = out1$lambda.MLE
theta.MLE = out1$cov.args.MLE$theta

#perform mKrig at MLE parameters
out1 = mKrig(x, y, lambda=lambda.MLE, theta=theta.MLE)
print("finished default case")

set.seed(1)
out2 = mKrigMAGMA.MLE(x, y, lambda=guessLambda, par.grid=par.grid, nGPUs=1)
lambda.MLE = out2$lambda.MLE
theta.MLE = out2$cov.args.MLE$theta

#perform mKrig at MLE parameters
out2 = mKrigMAGMA(x, y, lambda=lambda.MLE, theta=theta.MLE, nGPUs=1)

print("finished accelerated case with compact distance matrix")

set.seed(1)
out3 = mKrigMAGMA.MLE(x, y, lambda=guessLambda, par.grid=par.grid, 
                      cov.args= list(Distance="rdist"), nGPUs=1)
lambda.MLE = out3$lambda.MLE
theta.MLE = out3$cov.args.MLE$theta

#perform mKrig at MLE parameters
out3 = mKrigMAGMA(x, y, lambda=lambda.MLE, theta=theta.MLE, nGPUs=1)

#
##### test comatibility with other fields functions
#

temp1<- predict( out1)
temp2<- predict( out2)
temp3 = predict( out3)
test.for.zero( temp1, temp2, tag="predict compatibility: default versus accelerated")
test.for.zero( temp2, temp3, tag="predict compatibility: dist versus rdist")

#
##### test SE
#

temp1 = predictSE(out1)
temp2 = predictSE(out2)
temp3 = predictSE(out3)

test.for.zero( temp1, temp2, tag="predictSE compatibility: default versus accelerated")
test.for.zero( temp2, temp3, tag="predictSE compatibility: dist versus rdist")





cat("all done with mKrigMAGMA.MLE tests", fill=TRUE)
options( echo=TRUE)