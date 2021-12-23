library( ff-MAGMA)
options( echo=FALSE)
test.for.zero.flag<- 1

set.seed(1)

# test data
data( ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]

#
# Omit the NAs
good<- !is.na( y)
x<- x[good,]
y<- y[good]
compactDistMat = dist(x)
distMat = dist(x, x)

##### test using distance matrix
print("testing using distance matrix")

mKrig(x,y, cov.function = "stationary.fast.cov", lambda=2.0, theta=1.5) -> out1

mKrigMAGMA(x,y, cov.args= list(Covariance="Exponential", Distance="dist"), 
           lambda=2.0, theta=1.5, nGPUs=1) -> out2

#NOTE: compact distance matrix should not be used by user for fields compatibility reasons
mKrigMAGMA(x,y, cov.args= list(Covariance="Exponential", Distance="dist"), 
           lambda=2.0, theta=1.5, nGPUs=1, distMat=compactDistMat) -> out3

mKrigMAGMA(x,y, cov.args= list(Covariance="Exponential", Distance="dist"), 
           lambda=2.0, theta=1.5, nGPUs=1, distMat=distMat) -> out4

temp1<- predict( out1)
temp2<- predict( out2)
temp3 = predict( out3)
temp4 = predict( out4)
test.for.zero( temp1, temp2, tag="predict: no MAGMA versus MAGMA")
test.for.zero( temp2, temp3, tag="predict: no distance matrix versus compact distance matrix")
test.for.zero( temp2, temp4, tag="predict: no distance matrix versus distance matrix")

##### test SE
print("testing using predictSE")

temp1 = predictSE(out1)
temp2 = predictSE(out2)
#temp3 = predictSE(out3)
temp4 = predictSE(out4)

test.for.zero( temp1, temp2, tag="predictSE: no MAGMA versus MAGMA")
#test.for.zero( temp2, temp3, tag="predictSE: no distance matrix versus compact distance matrix")
test.for.zero( temp2, temp4, tag="predictSE: no distance matrix versus distance matrix")





cat("all done with mKrigMAGMA tests", fill=TRUE)
options( echo=TRUE)

