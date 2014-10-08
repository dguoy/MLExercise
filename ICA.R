library(fastICA)
library(audio)

sampleSize <- 50000
inputSize <- 2
alpha <- 0.0001

s1 <- load.wave("/Users/jp61130/Downloads/010000100mix1.wav")
s2 <- load.wave("/Users/jp61130/Downloads/010000100mix2.wav")

X <- matrix(c(s1, s2), nrow=inputSize, byrow=T)
meanX <- rowMeans(X)
X <- X - meanX
sigma <- X %*% t(X) / sampleSize
sigma.svd <- svd(sigma)
ZCAWhite <- sigma.svd$u %*% diag(1 / sqrt(sigma.svd$d + 0.001)) %*% t(sigma.svd$u)
X <- ZCAWhite %*% X

theta <- runif(inputSize * inputSize, min = -1, max = 1)

checkNumericalGradient(theta, costFunction, gradFunction)

optTheta <- optim(theta,
		function(theta) -costFunction(theta, X),
		function(theta) -gradFunction(theta, X),
		method = "L-BFGS-B", control = list(trace = 3, maxit = 2000))$par

W <- matrix(optTheta, nrow=inputSize)
S <- W %*% X

S[1, ] <- S[1, ] / max(abs(S[1, ]))
S[2, ] <- S[2, ] / max(abs(S[2, ]))

s1[1:sampleSize] <- S[1, ]
s2[1:sampleSize] <- S[2, ]

save.wave(s1, "/Users/jp61130/Downloads/s1.wav")
save.wave(s2, "/Users/jp61130/Downloads/s2.wav")

#***************************************************************************************************
costFunction <- function(theta, X) {
	sampleSize <- ncol(X)
	W <- matrix(theta, nrow=inputSize)
	sum(log(sigmoidGrad(W %*% X))) +
		sampleSize * log(det(W))
}
gradFunction <- function(theta, X) {
	sampleSize <- ncol(X)
	W <- matrix(theta, nrow=inputSize)
	gradW <- matrix(0, nrow=inputSize, ncol=inputSize)
	solveW <- solve(t(W))
	for(i in 1:sampleSize) {
		gradW <- gradW + (1 - 2 * sigmoid(W %*% X[, i])) %*% X[, i] + solveW
	}
	as.vector(gradW)
}

sigmoid <- function(z) 1 / (1 + exp(-z))
sigmoidGrad <- function(z) sigmoid(z) * (1 - sigmoid(z))
