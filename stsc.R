calc_rotation <- function(i, j, theta, .dim) {
    if (i == j) stop("i and j must be different", call = FALSE) 
    m <- diag(.dim)
    m[i, i] <- cos(theta)
    m[j, j] <- cos(theta)
    m[max(i, j), min(i, j)] <- sin(theta)
    m[min(i, j), max(i, j)] <- -sin(theta)
    m
}

calc_rotation_gradient <- function(i , j , theta , .dim) {
    if (i == j) stop("i and j must be different", call = FALSE) 
    m <- diag(0 , .dim)
    m[i, i] <- -sin(theta)
    m[j, j] <- -sin(theta)
    m[max(i, j), min(i, j)] <- cos(theta)
    m[min(i, j), max(i, j)] <- -cos(theta)
    m
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

get_rotation_matrix <- function(X, p){
    ij_list <- purrr::cross2(1:p, 1:p, .filter = `>=`)
    i <- purrr::map(ij_list, 1)
    j <- purrr::map(ij_list, 2)
    K <- length(ij_list)
    cost_and_grad <- function(theta_list) {
        U_list <- purrr::pmap(list(i = i, j = j, theta = theta_list), calc_rotation, .dim = p)
        V_list <- purrr::pmap(list(i = i, j = j, theta = theta_list), calc_rotation_gradient, .dim = p)
        R <- Reduce( "%*%" ,U_list , diag(p) )
        Z <- X %*% R
        # get the index of maximum element in each row
        mi <- apply( Z, 1, which.max)
        # i need the values of elements of maximum
        for (i in 1:length(mi)){
            M[i] <- Z[i ,mi[i]]}
        cost <- sum((Z/M)^2)
        grad <- vector(mode='numeric' , length=K)
        for (k in 1:K){
            A <- get_A_matrix(x, U_list, V_list, k, K)
            for (i in length(mi)){
                MA[i] <- A[i,mi[i]]}
            tmp <-  (Z / (M^2)) * A
            tmp <-  tmp - ((Z^2) / (M^3)) * MA
            tmp = 2 * sum(tmp)
            grad[k] = tmp}
        
        list(cost, grad)
    }
    theta_list_init <- vector(mode = "numeric" , length = p*(p-1)/2)
    opt <- optim(theta_list_init, cost_and_grad, method = "CG")
    
    
    list (opt$value , Reduce("%*%",generate_U_list(ij_list, opt$par , p) , diag(p)))
}
