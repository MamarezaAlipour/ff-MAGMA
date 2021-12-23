library( ff-MAGMA)
options( echo=FALSE)
test.for.zero.flag<- 1
data(ozone2)
y<- ozone2$y[16,]
x<- ozone2$lon.lat

#
# Omit the NAs

good<- !is.na( y)
x<- x[good,]
y<- y[good]

#####test that stationary.fast.cov and stationary.cov return the same result:

#with x1 == x2:

x1<- x[1:20,]
compactDistMat = dist(x1)
distMat = dist(x1, x1)
look<- stationary.cov(x1, theta=4)
look2 <- stationary.fast.cov(x1, theta=4, distMat = compactDistMat)
look3 <- stationary.fast.cov(x1, theta=4, distMat = distMat)
test.for.zero( look, look2, tag="stationary.cov versus stationary.fast.cov compact distMat")
test.for.zero( look, look3, tag="stationary.cov versus stationary.fast.cov matrix distMat")

#with x1 != x2:

x1<- x[1:20,]
x2=x[1:10,]
distMat = dist(x1, x2)
look<- stationary.cov(x1, x2, theta=4)
look2 <- stationary.fast.cov(x1, x2, theta=4, distMat = distMat)
test.for.zero( look, look2, tag="stationary.cov versus stationary.fast.cov asymmetric distMat")

##### test for correct value when using C argument:

Ctest<- rnorm(10)

#with x1 == x2:

x1 = x[1:10,]
compactDistMat = dist(x1)
distMat = dist(x1, x1)

temp<- stationary.cov( x1, C= Ctest, 
                       Covariance= "Wendland", 
                       k=2, dimension=2, theta=1.5 )

temp2 = stationary.fast.cov( x1, C= Ctest, 
                        Covariance= "Wendland", 
                        k=2, dimension=2, theta=1.5 )

temp3 = stationary.fast.cov( x1, C= Ctest, 
                             Covariance= "Wendland", 
                             k=2, dimension=2, theta=1.5, distMat=compactDistMat )

temp4 = stationary.fast.cov( x1, C= Ctest, 
                             Covariance= "Wendland", 
                             k=2, dimension=2, theta=1.5, distMat=distMat )

test.for.zero(temp, temp2, tag="stationary.cov vs stationary.fast.cov with C set")
test.for.zero(temp, temp3, tag="stationary.cov vs stationary.fast.cov with C set, compact distMat")
test.for.zero(temp, temp4, tag="stationary.cov vs stationary.fast.cov with C set, matrix distMat")

#with x1 != x2:

x1 = x
x2 = x[1:10,]

distMat = dist(x1, x1)

temp<- stationary.cov( x1, x2, C= Ctest, 
                       Covariance= "Wendland", 
                       k=2, dimension=2, theta=1.5 )

temp2 = stationary.fast.cov( x1, x2, C= Ctest, 
                             Covariance= "Wendland", 
                             k=2, dimension=2, theta=1.5 )

temp3 = stationary.fast.cov( x1, x2, C= Ctest, 
                             Covariance= "Wendland", 
                             k=2, dimension=2, theta=1.5, distMat=distMat )

test.for.zero(temp, temp2, tag="stationary.cov vs stationary.fast.cov with C set and asymmetric distMat not given")
test.for.zero(temp, temp3, tag="stationary.cov vs stationary.fast.cov with C set and asymmetric distMat given")

#
##### test covariance functions for onlyUpper=FALSE
#

x1=x[1:10,]
out1 = stationary.cov(x1, x1)
out2 = stationary.fast.cov(x1, x1)
out3 = Exp.fast.cov(x1, x1)

test.for.zero( out1, out2, tag="onlyUpper=FALSE: stationary.cov versus stationary.fast.cov")
test.for.zero( out2, out3, tag="onlyUpper=FALSE: stationary.fast.cov versus Exp.fast.cov")

#
##### test covariance functions for onlyUpper=TRUE
#
out2 = stationary.fast.cov(x1, onlyUpper=TRUE)
out3 = Exp.fast.cov(x1, onlyUpper=TRUE)

test.for.zero( out2[upper.tri(out2)], out3[upper.tri(out3)], tag="onlyUpper=TRUE: stationary.fast.cov with Exponential versus Exp.fast.cov")