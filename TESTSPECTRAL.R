nn <- 7
size <- dim(dist_mat)[1]
# intialize the k-nearest neighbor matrix
knn_mat <- matrix(0, nrow = size , ncol = size)
neighbor_index <- matrix(0, nrow = size , ncol = nn)
# find the 7 nearest neighbors for each point, considering local scale
for (i in 1:size){
    neighbor_index[i,] <- order(dist_mat[i,])[2:(nn + 1)]
    knn_mat[i,][neighbor_index] <- 1 
}

colnames(neighbor_index) <- seq(1 ,nn , 1)
row.names(neighbor_index) <- colnames(dist_mat)
# k'th position is chosen( here I consider 1th nearest neighbor)
# I find the local_scale.The local_scale is the distance between each Id and its 1th nearest neighbor
local_scale1 <- c()
for (i in 1:size){
    local_scale1 [i] <- dist_mat[i , neighbor_index[ i,nn]]
}

local_scale1 <- matrix(local_scale1 , nrow = size , ncol = 1)

# calculating sigma i * sigma j
sigma <- local_scale1 %*% t(local_scale1)

# Affinity matrix by considering local_scale ( 1th nearest neighbor)
affinitymatirx <- exp(-(dist_mat*dist_mat)/sigma)
diag(affinitymatirx) <- 0
# A diagonal matrix which diagonal elements are degree of each variable in affinity matrix
deg <- colSums(affinitymatirx)
D <- diag(1 / sqrt(deg))
# laplacian matrix
L <- diag(deg) - affinitymatirx
# Normalized laplacian matrix
N_L <- D %*% L %*% D

# calculating eigenvalues and eigenvectors of laplacian matrix
eigenvalues <- eigen(N_L)$values
largest_gap <- abs(diff(eigenvalues))


plot(1:10 , eigenvalues[1:10] , main = "largest eigenvalues of matrix" , type="p" , col="blue" ,
     pch = 21 , bg="blue" )


gaps <- abs(diff(eigenvalues))
ndx <- order(gaps, decreasing = T)[1:nn]
gaps[ndx]
ndx

plot( 1:length(eigenvalues) , sort(eigenvalues , decreasing = F) , main = "largest eigenvalues of matrix" , type="p",
      col="blue" ,pch = 16 , bg="red" )
