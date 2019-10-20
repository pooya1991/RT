# stsc_np.py

Given_rotation <- function( i, j , theta , dimension){
    g <- diag(dimension)
    c <- cos(theta)
    s <- sin(theta)
    g[i,i] <- c
    g[j,j] <- c
    if(i >j){
        g[i,j] <- s
        g[j,i] <- -s
    }else if (i<j){
        g[i,j] <- -s
        g[j,i] <- s
    }else {
        stop("i and j must be different")
    } 
    g
}


Given_rotation_gradient <- function(i , j , theta , dimension){
    g <- matrix(0 , dimension, dimension)
    c <- cos(theta)
    s <- sin(theta)
    g[i,i] <- -s
    g[j,j] <- -s
    if(i >j){
        g[i,j] <- c
        g[j,i] <- -c
    }else if (i<j){
        g[i,j] <- -c
        g[j,i] <- c
    }else {
        stop("i and j must be different")
    } 
    g
}


generate_U_list <- function(ij_list , theta_list , dimension){
    k <- 0
    g_u_list <- list()
    for(i in 1:length(theta_list)){
        k <- k+1
        g_u_list[[k]]<- Given_rotation(ij_list[i,1], ij_list[i,2] , theta_list[i], dimension)
    }
    g_u_list
}


generate_V_list <- function(ij_list , theta_list , dimension){
    k <- 0
    g_v_list <- list()
    for (i in 1:length(theta_list)){
        k <- k+1
        g_v_list[[k]] <- Given_rotation_gradient(ij_list[i,1], ij_list[i,2] , theta_list[i], dimension)
    }
    g_v_list
}

get_U_ab <- function(a,b,U_list , K){
    I <- diag(length(U_list[[1]][,1]))
    if(a==b){
        if(a < K & a!=0){
            U_list[a]
           }else {
                return(I)}}  
     else if(a>b){
           return(I)}
     else {
            Reduce("%*%" , U_list[a:b] , I)
        }
}


get_A_matrix <- function(X , U_list , V_list , k , K){
    Ul <- get_U_ab(1, k, U_list, K)
    V  <- V_list[[k]]
    Ur <-  get_U_ab(k + 1, K, U_list, K)
    
    X %*% Ul %*% V %*% Ur
}  

get_rotation_matrix <- function(X,C){
    k <- 0
    ij_list = list()
    for (i in 1:length(C)){
        for ( j in 1:length(C)){
            if (i<j){
                k <- k+1
                ij_list[[k]] <- c(i , j)
            }
        }
    }
    
    ij_list <- do.call(rbind, ij_list)
    K <- nrow(ij_list)
    
    cost_and_grad <- function(theta_list){
        U_list <- generate_U_list(ij_list, theta_list, C)
        V_list <- generate_V_list(ij_list, theta_list ,C)
        R <- Reduce( "%*%" ,U_list , diag(C) )
        Z <- X %*% R
        mi <- apply( Z, 2, which.max)
        for (i in 1:length(mi)){
            M[i] <- aa[i ,mi[i]]}
        cost <- sum((Z/M)^2)
        grad <- vector(mode='numeric' , length=K)
        for (k in 1:K){
            A <- get_A_matrix(X, U_list, V_list, k, K)
            for (i in length(mi)){
                MA[i] <- A[i,mi[i]]}
            tmp <-  (Z / (M^2)) * A
            tmp <-  tmp - ((Z^2) / (M^3)) * MA
            tmp = 2 * sum(tmp)
            grad[k] = tmp}

        list(cost, grad)
    }
    theta_list_init <- vector(mode = "numeric" , length = C*(C-1)/2)
    opt <- optim(theta_list_init, cost_and_grad, method = "CG")
    

    list (opt$value , Reduce("%*%",generate_U_list(ij_list, opt$par , C) , diag(C)))
}
