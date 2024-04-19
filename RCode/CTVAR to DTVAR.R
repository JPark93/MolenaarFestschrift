# Creating CT-VAR to DT-Var
library(qgraph)
QtoPsi=function(Q = NULL, base.A = NULL, DeltaT = 1.0){
  A_hashtag = kronecker(base.A, diag(rep(1, 4))) + kronecker(diag(rep(1, 4)), base.A)
  Qtorow = matrix(Q, ncol = 1)
  Psi = round(matrix(solve(A_hashtag) %*% (expm(A_hashtag * DeltaT) - 
                                       diag(rep(1, 4^2))) %*% Qtorow, ncol = 4),2)
  return(Psi)
}
  base.A=matrix(0, 4, 4)
  diag(base.A) = -0.5
  base.A[1,2]=-0.3
  base.A[3,1]=-0.3
  base.A[2,4]=0.3

# Graph of OU-Model - A Matrix
par(mfrow=c(1, 3))
  qgraph(base.A, diag=TRUE, layout = 'circle', fade = FALSE,
         theme = 'gray',
         edge.labels = base.A, edge.label.position=0.75)
# Graph of VAR(1) - Delta-t = 1.0
  delt1 = round(expm(base.A * 1.0), 2)
  qgraph(delt1, diag=TRUE, layout = 'circle', fade = FALSE, 
         theme = 'gray',
         edge.labels = delt1, edge.label.position=0.75)
# Graph of VAR(1) - Delta-t = 10.0
  delt10 = round(expm(base.A * 10.0), 2)
  qgraph(delt10, diag=TRUE, layout = 'circle', fade = FALSE, 
         theme = 'gray',
         edge.labels = delt10, edge.label.position=0.75)

  
  
# Graph of OU-Model - Q Matrix
  base.Q=diag(1, 4)
  qgraph(base.Q, diag=TRUE, layout = 'circle', fade = FALSE,
         theme = 'gray',
         edge.labels = base.Q, edge.label.position=0.75)
# Graph of VAR(1) - Delta-t = 1.0 - Psi
  delt1.Psi = QtoPsi(Q = base.Q, base.A = base.A, DeltaT = 1.0)
  qgraph(delt1.Psi, diag=TRUE, layout = 'circle', fade = FALSE,
         theme = 'gray',
         edge.labels = delt1.Psi, edge.label.position=0.75)
# Graph of VAR(1) - Delta-t = 10.0 - Psi
  delt10.Psi = QtoPsi(Q = base.Q, base.A = base.A, DeltaT = 10.0)
  qgraph(delt10.Psi, diag=TRUE, layout = 'circle', fade = FALSE,
         theme = 'gray',
         edge.labels = delt10.Psi, edge.label.position=0.75)
  