library(R.matlab)
library(ggplot2)

movies <- readMat("data/ex8_movies.mat")
movieParams <- readMat("data/ex8_movieParams.mat")

J <- function(Y, R, X, Theta, lambda) {
	return (0.5 * sum((X %*% t(Theta) - Y)^2 * R) + (lambda / 2) * (sum(X^2) + sum(Theta^2)))
}

alpha <- 0.0005
lambda <- 1.5
num_users <- movieParams$num.users
num_movies <- movieParams$num.movies
num_features <- movieParams$num.features
# X <- movieParams$X
# Theta <- movieParams$Theta
X <- matrix(rnorm(num_movies * num_features), num_movies, num_features)
Theta <- matrix(rnorm(num_users * num_features), num_users, num_features)
Y <- movies$Y
R <- movies$R

for(i in 1:1000) {
	X <- X - alpha * (((X %*% t(Theta) - Y) * R) %*% Theta + lambda * X)
	Theta <- Theta - alpha * (t((X %*% t(Theta) - Y) * R) %*% X + lambda * Theta)
	print(J(Y, R, X, Theta, lambda))
}

movie_cluster <- kmeans(X, 100)$cluster

