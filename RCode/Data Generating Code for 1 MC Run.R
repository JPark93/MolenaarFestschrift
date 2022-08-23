library(Matrix); 
source('~/work/StateSpaceFunctions.R')
thedir = '~/scratch/'
# Making relevant directories for everything to save in one place
dir.create(file.path(paste0(thedir,'/JUL22/', sep = '')), showWarnings = F)
dir.create(file.path(paste0(thedir,'/JUL22/OUModels/', sep = '')), showWarnings = F)
dir.create(file.path(paste0(thedir,'/JUL22/AMats/', sep = '')), showWarnings = F)
# Simulating N = 50 (NG1 = 25; NG2 = 25) ----
ne = 4 # Variables
allY = NULL
set.seed(010593+x)
# Large EF condition
sd = 0.0; d = 'Large'; y = NULL
maxT <- 10000# Total number of time points
Deltat = .1 #Resolution of Data Generation
ne = 4 #Number of latent variables
ny = 4 #Number of observed variables
TimeSeq  <- seq(0, (maxT-1)*Deltat,by=Deltat) #Time indices, with, Delta t = 0.1
mu = c(0,0,0,0) #Home base values
a0 = c(2,-1, .5,-.5)#c(-3,1.5, 3,5) #Initial condition means
P0 = diag(1,ne) #Initial condition co-variance matrix
value = .9 # 0.9 for the large condition
A1 =  -matrix(c(0.5, rnorm(1,value,sd), 0.0, 0.0,
       0.0, 0.5, 0.0, rnorm(1,-value,sd),
       rnorm(1,value,sd), 0.0, 0.5, 0.0,
       0.0, 0.0, 0.0, 0.5), 4, 4,byrow=TRUE)
A2 = -matrix(c(0.5, rnorm(1,-value,sd), 0.0, 0.0,
      0.0, 0.5, 0.0, 0.0,
      rnorm(1,value,sd), 0.0, 0.5, 0.0,
      0.0, 0.0, rnorm(1,-value,sd), 0.5), 4, 4,byrow=TRUE)
amt = c((x*25-24):(x*25),(x*25-24):(x*25)+2500)
# Simulating 50 subjects; total MC Study simulated 5,000 subjects
for(i in amt){
      if(i < (5000/2)+1) A = A1 else A = A2
      write.csv(A, file=paste0(thedir,"/JUL22/AMats/AMAT SUBJ", i,"OU_cohen ",d,".csv"), row.names = FALSE)
      #b = B x mu 
      b = -A %*%matrix(mu,ncol=1)
      # G - Matrix square-root of process noise covariance matrix
      G = diag(1,ne)
      # Maybe add covariances/variances
      Q = G%*%t(G) # Process noise covariance matrix
      #Components related to the discrete-time solution of the OU model
      A_del = expm(A*Deltat)
      b_del= solve(A)%*%(expm(A*Deltat)-diag(1,ne))%*%b 
      #Process noise covariance matrix, Psi
      A_hashtag = kronecker(A,diag(rep(1,ne))) + kronecker(diag(rep(1,ne)),A)
      Qtorow = matrix(Q,ncol=1)
      Psi = matrix(solve(A_hashtag)%*%(expm(A_hashtag*Deltat)-diag(rep(1,ne^2)))%*%Qtorow,ncol=ne)
      #Measurement-related parameters
      Lambda = diag(1,ne) #Factor loading matrix
      R = diag(1e-5,ne) #Measurement error covariance matrix
      thedat = simStateSpaceData(a0 = a0, P0=P0, Q=Psi,R=R,Phi =A_del,Lambda = Lambda,
                                 alpha=b_del,tau = rep(0,ny),nt = maxT,np = 1,ne = ne,ny = ny,nx=0,npad=50)
      y = data.frame(ID = rep(i,each=maxT),
                     Time = TimeSeq,# rep((0:(maxT-1))*Deltat,1),
                     y1 = thedat$obsData[,1],
                     y2 = thedat$obsData[,2],
                     y3 = thedat$obsData[,3],
                     y4 = thedat$obsData[,4])
    write.table(y,file=paste0(thedir,"/JUL22/OUModels/SUBJ", i,"OU_cohen ",d,".txt"),sep=",")
}