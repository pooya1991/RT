#stsc_manopt
# Here we are defining a cost function and our goal is to minimize the cost function on a specefic manifold
# "Stiefel" manifold.The Stiefel manifold is the set of all orthonormal k-frames real coordinate space.
# Manifold St(n,p) of n x p matrices whose columns are orthonormal.

# library for optimization on the steifel manifold
library(rstiefel)

library(ManifoldOptim)

get_rotation_matrix <- function(X,C){
    cost <- function(R){
        Z <- X %*% R
        M <- apply(Z,2,max)
        return(sum((Z/M)^2))
    # This function siumlates a random orthonormal matrix from the uniform distribution on 
    # the Stiefel manifold. an C*C orthonormal matrix
    manifold <- rustiefel(C, C)
    opt <- optStiefel(cost, manifold, Vinit=rustiefel(N, P),
                      method="curvilinear",
                      searchParams=list(rho1=0.1, rho2=0.9, tau=1),tol=1e-4) ### need work.!!!!
    list(opt , cost(opt))
    }
}