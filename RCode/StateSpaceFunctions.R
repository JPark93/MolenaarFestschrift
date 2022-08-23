# This is a function for simulating data following any arbitrary linear state-space models
# Input arguments:
# a0 = Mean of initial latent variable values at time t_0
# P0 = Covariance matrix of initial latent variable values at time t_0
# Q = Process noise covariance matrix
# R = Measurement error covariance matrix
# Phi = Transition matrix relating the values of the latent variables at time t-1 to those at time t.
# Z1 = matrix of regression coefficients linking the time-varying covariates to the latent variables
# Z2 = matrix of regression coefficients linking the time-varying covariates to the observed variables
# Lambda = Factor loading matrix linking the latent variables to the observed variables
# alpha = Vector of intercepts for dynamic model
# tau = Vector of intercepts for measurement model
# nt = Number of time points per unit
# np = Number of units (e.g., subjects)
# ne = Number of latent variables
# ny = Number of manifest variables
# nx = Number of time-varying covariates/predictors
# npad = Number of time points you want to simulate to wash out hypothetical initial conditions before keeping for further analyses
# x = Matrix of covariates containing np*(nt + npad) x nx covariate values in long format; with the first subject's
# values going into the first (nt + npad) rows, followed by the second subject's data, and so on.
# isWrite = A binary flag stating whether you want to save the simulated data to an external file. By default,
# the observed data will be saved to your working directory, with the name of 'StateSpaceObs.csv', and
# the latent variabe values will be saved to the same directory with the name of 'StateSpaceLatent.csv'.

simStateSpaceData = function(a0, P0, Q,R,Phi,Lambda,Z1=NULL,Z2=NULL,alpha,tau,nt,np,ne,ny,nx,npad,x=NULL){
ist=npad+1   
ntt=nt+npad
Qs = chol(Q)  #Cholesky (matrix square root) of Q
if (det(R) > 0){Rs = chol(R)} else {Rs = matrix(0,ny,ny)} #Cholesky (matrix square root) of R

a=matrix(0,ntt,ne) #Place holder for latent variable values for one unit (e.g., subject)
y=matrix(0,ntt,ny) #Place holder for manifest variable scores for one unit (e.g., subject)
yall=matrix(0,nt*np,ny+nx) # Place holder for latent variables for all units (e.g., all subjects)
all = matrix(0,nt*np,ne) #Place holder for manifest variable scores for all units (e.g., all subjects)

for (j in 1:np){
a[1,1:ne] = matrix(a0 + rnorm(ne)%*%chol(P0),ncol=ne)
etmp=t(rnorm(ny)%*%Rs)
y[1,1:ny] = Lambda%*%matrix(a[1,1:ne],ncol=1)+etmp+tau
  if (!is.null(Z2)) {y[1,1:ny] = y[1,1:ny]  + Z2%*%matrix(x[(1+(j-1)*nt),],ncol=ny)}
  for (i in 2:ntt) {
    ztmp=t(rnorm(ne)%*%Qs)
    etmp=t(rnorm(ny)%*%Rs)
    atmp=Phi%*%matrix(a[i-1,1:ne],ncol=1)+ztmp+ alpha 
    if (!is.null(Z1)) {atmp = atmp + Z1%*%matrix(x[(i+(j-1)*nt),],ncol=1)}
    a[i,1:ne]=matrix(atmp,ncol=ne)
    ytmp=Lambda%*%atmp+etmp+tau
    if (!is.null(Z2)) {ytmp = ytmp + Z2%*%matrix(x[(i+(j-1)*nt),],ncol=1)}
    y[i,1:ny]=matrix(ytmp,ncol=ny)
  }
  #Stack all subjects' data together
if (nx>0){
  yall[(1+(j-1)*nt):(j*nt),1:(ny+ny)] = matrix(cbind(y[(ist:ntt),1:ny],x[ist:ntt]),nt,ny+nx)
}else{
  yall[(1+(j-1)*nt):(j*nt),1:ny] = matrix(y[(ist:ntt),1:ny],nt,ny)
}
  all[(1+(j-1)*nt):(j*nt),1:ne] = a[(ist:ntt),1:ne]
}

return(list(obsData = yall, stateData = all))

}
