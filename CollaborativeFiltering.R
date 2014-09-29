library(R.matlab)
library(ggplot2)

movies <- readMat("data/ex8_movies.mat")
movieParams <- readMat("data/ex8_movieParams.mat")

theta <- 0.005 * runif(movieParams$num.movies * movieParams$num.features + movieParams$num.users * movieParams$num.features)

optTheta <- optim(theta,
		function(theta) costFunction(theta, movies$Y, movies$R, movieParams$num.users, movieParams$num.movies, movieParams$num.features),
		function(theta) gradFunction(theta, movies$Y, movies$R, movieParams$num.users, movieParams$num.movies, movieParams$num.features),
		method = "L-BFGS-B", control = list(trace = 3, maxit = 2000))$par

thetaMovie <- matrix(optTheta[1 : (movieParams$num.movies*movieParams$num.features)], movieParams$num.movies, movieParams$num.features)
movieCluster <- kmeans(thetaMovie, 100)$cluster

#************************************************ Function ***********************************************************
costFunction <- function(theta, Y, R, numUsers, numMovies, numFeatures, lambda=1.5) {
	thetaMovie <- matrix(theta[1 : (numMovies*numFeatures)], numMovies, numFeatures)
	thetaUser <- matrix(theta[(numMovies*numFeatures+1) : length(theta)], numUsers, numFeatures)
	(1 / 2) * sum((thetaMovie %*% t(thetaUser) - Y)^2 * R) +
			(lambda / 2) * (sum(thetaMovie^2) + sum(thetaUser^2))
}
gradFunction <- function(theta, Y, R, numUsers, numMovies, numFeatures, lambda=1.5) {
	thetaMovie <- matrix(theta[1 : (numMovies*numFeatures)], numMovies, numFeatures)
	thetaUser <- matrix(theta[(numMovies*numFeatures+1) : length(theta)], numUsers, numFeatures)

	gradMovie <- ((thetaMovie %*% t(thetaUser) - Y) * R) %*% thetaUser + lambda * thetaMovie
	gradUser <- t((thetaMovie %*% t(thetaUser) - Y) * R) %*% thetaMovie + lambda * thetaUser

	c(as.vector(gradMovie), as.vector(gradUser))
}
