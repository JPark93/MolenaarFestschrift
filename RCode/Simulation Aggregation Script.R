library(ggplot2); library(mclust)
library(dplyr); library(Matrix)
library(cowplot); library(matrixcalc)
std.var = function(Phihat = NULL, Psihat = NULL){
  true.phi = matrix(NA, dim(Phihat)[1], dim(Phihat)[2])
  I = diag(1, dim(Phihat)[1] * dim(Phihat)[2])
  vecP0 = solve(I - kronecker(Phihat, Phihat)) %*% vec(Psihat)
  P0 = matrix(vecP0, dim(Phihat)[1], dim(Phihat)[2])
  for(j in 1:dim(Phihat)[1]){
    for(k in 1:dim(Phihat)[2]){
      true.phi[j,k] = Phihat[j,k] * (sqrt(P0[k,k])/sqrt(P0[j,j]))
    }
  }
  Psi = cov2cor(Psihat)
  return(list(Phi = true.phi, Psi = Psi))
}
thedirs = c("/Users/JonathanPark/Desktop/JUL22/")
dirs = c("Last14","Last100")#,"All")
EFs = c('SMALL','MEDIUM','LARGE','EXPL')
ne = 4; njobs = 100;
N = 50; Sim.Res = NULL
gimme.crit = 0.75; gvar.crit = 0.00
bias.var = function(matA = NULL, matB = NULL, 
                    rel = FALSE, nonzero = TRUE,
                    RMSE = TRUE, sep.diag = TRUE){
  odiag = function(x){x[col(x) != row(x)]}
  # Mat B is always the reference matrix
  if(is.null(matA) | is.null(matB)){
    object = data.frame(Bias = NA, rBias = NA, RMSE = NA, SD = NA,
                      dBias = NA, odBias = NA,
                      drBias = NA, odrBias = NA)
    row.names(object) = 'Values'
  }else{
  if(is.list(matA) & is.list(matB)){
    bias=list();rbias=list()
    dbias=list();drbias=list()
    odbias=list();odrbias=list()
    for(i in 1:length(matA)){
      if(nonzero == TRUE){
        logic = matB[[i]]!=0
        odlogic = !c(as.matrix(matB[[i]])) %in% diag(as.matrix(matB[[i]])) & c(as.matrix(matB[[i]])) != 0
      }else logic = matrix(TRUE, nrow(matB[[i]]), ncol(matB[[i]]))
      bias[[i]] = (matA[[i]] - matB[[i]])[logic]
      dbias[[i]] = diag(as.matrix(matA[[i]])) - diag(as.matrix(matB[[i]]))
      odbias[[i]] = c(as.matrix(matA[[i]]))[odlogic] - c(as.matrix(matB[[i]]))[odlogic]
    }
    if(rel == TRUE){
      for(i in 1:length(matA)){
        if(nonzero == TRUE){
          logic = matB[[i]]!=0
        }else logic = matrix(TRUE, nrow(matB[[i]]), ncol(matB[[i]]))
        rbias[[i]] = ((matA[[i]] - matB[[i]])/matB[[i]])[logic]
        drbias[[i]] = (diag(as.matrix(matA[[i]])) - diag(as.matrix(matB[[i]])))/diag(as.matrix(matB[[i]]))
        odrbias[[i]] = (c(as.matrix(matA[[i]]))[odlogic] - c(as.matrix(matB[[i]]))[odlogic])/
          c(as.matrix(matB[[i]]))[odlogic]
      }
      rbias = median(unlist(rbias), na.rm = TRUE)
    }
      if(RMSE == TRUE){
        rmse = sqrt(median(unlist(bias)^2, na.rm = TRUE))  
      }
      sdB = sd(unlist(bias), na.rm = TRUE)
      bias = median(unlist(bias), na.rm = TRUE)
    }else{
      if(nonzero == TRUE){
        logic = matB!=0
      }else logic = matrix(TRUE, nrow(matB), ncol(matB))
      bias = median(matA - matB, na.rm = TRUE)
      sdB = sd(matA - matB, na.rm = TRUE)
      rmse = sqrt(median(bias, na.rm = TRUE)^2)
      rbias = median(((matA - matB)/matB)[logic], na.rm = TRUE)
    }
  object = data.frame(Bias = bias, rBias = rbias, RMSE = rmse, SD = sdB,
                      dBias = median(unlist(dbias),na.rm=TRUE), odBias = median(unlist(odbias),na.rm=TRUE),
                      drBias = median(unlist(drbias),na.rm=TRUE), odrBias = median(unlist(odrbias),na.rm=TRUE))
  row.names(object) = 'Values'
  }
  return(object)
}
for(w in thedirs){
for(z in dirs){
  for(y in EFs){
    for (theDeltat in c(0.1,0.5,1.0,10)){
      estBs = list(); estBg = list()
      estPs = list(); estPg = list()
      ari.list = list(); I = diag(1,ne)
      aris = NULL; A=list(); Q=list(); transA = list()
      for(l in 1:njobs){
          amt = c((l*25-24):(l*25),(l*25-24):(l*25)+2500)
          if(!file.exists(paste0(w,'/',z,'/',y,'/',theDeltat,'/job ',l,'/'))){
            message(paste0('Job ', l, '. Config: ', z, ' DelT = ', theDeltat))
            next
          }else{
          ari.list[[l]] = cbind(rep(1:2, each = N/2),
                                read.csv(paste0(w,'/',z,'/',y,'/',theDeltat,'/job ', l,'/',gimme.crit,'N',gvar.crit, 
                                                '/graphical VAR/Subgroup Solutions/Group Membership.csv'),header=FALSE)[,2],
                                read.csv(paste0(w,'/',z,'/',y,'/',theDeltat,'/job ', l,'/',gimme.crit,'N',gvar.crit,
                                                '/summaryFit.csv'))[c(1,12,23,34,45,47:50,2:11,
                                                                      13:22,24:33,35:44,46),14])
          aris = rbind(aris,cbind(adjustedRandIndex(ari.list[[l]][,1],ari.list[[l]][,2]), 
                                  adjustedRandIndex(ari.list[[l]][,1],ari.list[[l]][,3])))
          filedir = paste0(w,'/',z,'/',y,'/',theDeltat,'/job ',l,'/',gimme.crit,'N',gvar.crit,'/')
          #SMC: As good hygiene, remove all files in that directory
          for (m in 1:N){
            A[[amt[m]]] = read.csv(paste0(w,'/AMats/AMAT SUBJ', amt[[m]], 'OU_cohen ',y,'.csv'))
            A_hashtag = kronecker(as.matrix(A[[amt[m]]]),diag(rep(1,4))) + kronecker(diag(rep(1,ne)),as.matrix(A[[amt[m]]]))
            # Discrete time equivalent process noise covariance matrices
            Q[[amt[m]]] = matrix(solve(A_hashtag)%*%(expm(A_hashtag*theDeltat)-diag(rep(1,ne^2)))%*%matrix(diag(1,ne),ncol=1),ncol=ne)
            if(!file.exists(paste0(filedir,'individual/subj',m,'Betas.csv'))){
              estBs[[amt[m]]] = NA
              estPs[[amt[m]]] = NA
            }else{
              temp=as.matrix(read.csv(paste0(filedir,'individual/subj',m,'Betas.csv'))[2:(ne*2+1)])
              estBs[[amt[m]]] = solve(I - temp[,(ne+1):(ne*2)]) %*% temp[,1:ne]
              # Reading in sVAR Psi matrix
              temp1 = as.matrix(read.csv(paste0(filedir,'individual/subj',m,'PsiUnstd.csv'))[(ne+1):(ne*2),(ne+1):(ne*2)+1])
              # Turning sVAR Psi to VAR psi*
              estPs[[amt[m]]] = solve(I - temp[,(ne+1):(ne*2)]) %*% temp1 %*% t(solve(I - temp[,(ne+1):(ne*2)])) 
            }
            # Reading in phi and psi on VAR metric from gVAR
            if(!file.exists(paste0(filedir,'graphical VAR/Individual Solutions/Subj',m,'beta.csv'))){
              estBg[[amt[m]]] = NA
              estPg[[amt[m]]] = NA
            }else{
              estBg[[amt[m]]] = read.csv(paste0(filedir,'graphical VAR/Individual Solutions/Subj',m,'beta.csv'),header=FALSE)[,2:(ne+1)]
              estPg[[amt[m]]] = solve(read.csv(paste0(filedir,'graphical VAR/Individual Solutions/Subj',m,'kappa.csv'),header=FALSE))  
            }
          }
      }
        }
      for(i in 1:length(A)){
        if(is.null(A[[i]])){
          transA[[i]] = matrix(NA, 4, 4)
        }else transA[[i]] = as.matrix(expm(as.matrix(A[[i]])*theDeltat))
      }
      for(i in 1:length(estBg)){
        if(is.null(estBg[[i]])){
          estBg[[i]] = matrix(NA, 4, 4)
        }else estBg[[i]] = estBg[[i]]
      }
      for(i in 1:length(estBs)){
        if(is.null(estBs[[i]])){
          estBs[[i]] = matrix(NA, 4, 4)
        }else estBs[[i]] = estBs[[i]]
      }
      for(i in 1:length(estPg)){
        if(is.null(estPg[[i]])){
          estPg[[i]] = matrix(NA, 4, 4)
        }else estPg[[i]] = estPg[[i]]
      }
      for(i in 1:length(estPs)){
        if(is.null(estPs[[i]])){
          estPs[[i]] = matrix(NA, 4, 4)
        }else estPs[[i]] = estPs[[i]]
      }
      for(i in 1:length(Q)){
        if(is.null(Q[[i]])){
          Q[[i]] = matrix(NA, 4, 4)
        }else Q[[i]] = Q[[i]]
      }
      TPhi = list()
      TPsi = list()
      for(i in 1:length(transA)){
        if(is.null(Q[[i]]) | is.null(transA[[i]]) | 
           any(is.na(Q[[i]])) | any(is.na(transA[[i]]))){
          TPhi[[i]] = matrix(NA, 4, 4)
          TPsi[[i]] = matrix(NA, 4, 4)
          next
          }else{
            TPhi[[i]] = std.var(transA[[i]], Q[[i]])$Phi
            TPsi[[i]] = std.var(transA[[i]], Q[[i]])$Psi
          }
      }
      plot.writer(estimate = estBg, true = TPhi,
            REL = TRUE, TDIR = paste0('~/Desktop/Figures/scgVAR/Phi/',z, y,theDeltat,'/'),
            TITLE = '', TITLE.REL = '',
            TITLE.HEAT = '', TITLE.HEAT.2 = '')
      plot.writer(estimate = estBs, true = TPhi,
            REL = TRUE, TDIR = paste0('~/Desktop/Figures/GIMME/Phi/',z, y,theDeltat,'/'),
            TITLE = '', TITLE.REL = '',
            TITLE.HEAT = '', TITLE.HEAT.2 = '')
      plot.writer(estimate = estPg, true = TPsi,
            REL = TRUE, TDIR = paste0('~/Desktop/Figures/scgVAR/Psi/',z, y,theDeltat,'/'),
            TITLE = '', TITLE.REL = '',
            TITLE.HEAT = '', TITLE.HEAT.2 = '')
      plot.writer(estimate = estPs, true = TPsi,
            REL = TRUE, TDIR = paste0('~/Desktop/Figures/GIMME/Psi/',z, y,theDeltat,'/'),
            TITLE = '', TITLE.REL = '',
            TITLE.HEAT = '', TITLE.HEAT.2 = '')
      biasBforG = bias.var(estBg, TPhi, rel = TRUE, nonzero = TRUE, RMSE = TRUE)  
      biasBforS = bias.var(estBs, TPhi, rel = TRUE, nonzero = TRUE, RMSE = TRUE)
      biasPforG = bias.var(estPg, TPsi, rel = TRUE, nonzero = TRUE, RMSE = TRUE)  
      biasPforS = bias.var(estPs, TPsi, rel = TRUE, nonzero = TRUE, RMSE = TRUE)
      Sim.Sum = data.frame(BiasB = c(biasBforG$Bias, biasBforS$Bias),
                           dBiasB = c(biasBforG$dBias, biasBforS$dBias),
                           odBiasB = c(biasBforG$odBias, biasBforS$odBias),
                           BiasP = c(biasPforG$Bias, biasPforS$Bias),
                           dBiasP = c(biasPforG$dBias, biasPforS$dBias),
                           odBiasP = c(biasPforG$odBias, biasPforS$odBias),
                           rBiasB = c(biasBforG$rBias,biasBforS$rBias),
                           drBiasB = c(biasBforG$drBias,biasBforS$drBias),
                           odrBiasB = c(biasBforG$odrBias,biasBforS$odrBias),
                           rBiasP  = c(biasPforG$rBias, biasPforS$rBias),
                           drBiasP  = c(biasPforG$drBias, biasPforS$drBias),
                           odrBiasP  = c(biasPforG$odrBias, biasPforS$odrBias),
                           ARI = colMeans(aris, na.rm = TRUE),
                           sdB = c(biasBforG$SD,biasBforS$SD), sdP = c(biasPforG$SD,biasPforS$SD),
                           rmseB = c(biasBforG$RMSE,biasBforS$RMSE), rmseP = c(biasPforG$RMSE,biasPforS$RMSE),
                   Model = factor(rep(1:2,each = 1), labels = c('scgVAR','S-GIMME')),
                   DelT = factor(rep(1, times = 2), labels = theDeltat),
                   EFS = factor(rep(1, times = 2), labels = y),
                   Cond = factor(rep(1, times = 2), labels = z),
                   Setting = factor(rep(1, times = 2), labels = w))
      Sim.Res = rbind(Sim.Res, Sim.Sum)
      message(paste0(theDeltat,y,z))
    }
  }
}
}
Sim.Res %>% mutate_if(is.numeric, round, digits = 2)
beepr::beep(5)
Cond.labs = c('T = 14', 'T = 100')
names(Cond.labs) = c('Last14', 'Last100')
EFS.labs = c('Small','Large','Small','Large')
names(EFS.labs) = c('SMALL', 'MEDIUM','LARGE','EXPL')
p1 = ggplot(data = Sim.Res[Sim.Res$EFS == 'SMALL' | 
                             Sim.Res$EFS == 'MEDIUM',], 
            aes(x=DelT,y=rmseP,color=Model,group=Model)) + 
    geom_point() + geom_line() + theme_bw() +
    geom_hline(yintercept=0.00, linetype="dashed", color = "red", size=1) + 
    labs(x=expression(""*Delta*t*""), 
         y=expression(" "*RMSE*" in "*Psi*"")) +
         # y = 'ARI') +
  lims(y = c(-0.0,0.5)) +
  facet_grid(cols = vars(Cond), rows = vars(EFS), 
             labeller = labeller(Cond = Cond.labs, EFS = EFS.labs), 
             scales = 'free_y') #+ coord_cartesian(ylim=c(-0.3, 0.05))

p2 = ggplot(data = Sim.Res[Sim.Res$EFS == 'LARGE' | 
                             Sim.Res$EFS == 'EXPL',], 
            aes(x=DelT,y=rmseP,color=Model,group=Model)) + 
    geom_point() + geom_line() + theme_bw() +
    geom_hline(yintercept=0.00, linetype="dashed", color = "red", size=1) + 
    labs(x=expression(""*Delta*t*""), 
         y=expression(" "*RMSE*" in "*Psi*"")) +
         # y = 'ARI') +
  lims(y = c(-0.0,0.5)) +
  facet_grid(cols = vars(Cond), rows = vars(EFS), 
             labeller = labeller(Cond = Cond.labs, EFS = EFS.labs), 
             scales = 'free_y') #+ coord_cartesian(ylim=c(-0.3, 0.05))

plot_grid(p1, p2, labels = c("           Stable", 
                             "Near Unstable"), nrow = 2,
          label_size = 12, hjust = -3.25)

# Cond.labs = c('T = 14', 'T = 100')
# names(Cond.labs) = c('Last14', 'Last100')
# EFS.labs = c('Small','Medium','Large','Explosive')
# names(EFS.labs) = c('SMALL', 'MEDIUM','LARGE','EXPL')
# ggplot(data = Sim.Res, aes(x=DelT,y=rBiasB,color=Model,group=Model)) + 
#     geom_point() + geom_line() + theme_bw() +
#     geom_hline(yintercept=0.00, 
#                linetype="dashed", color = "red", size=1) + 
#   labs(x=expression(""*Delta*t*""), 
#        y=expression(" "*rBias*" in "*Phi*"")) +
#        #y='Adjusted Rand Index (ARI)') +
#   # lims(y = c(-1.0,0.0)) +
#   facet_grid(cols = vars(Cond), rows = vars(EFS), 
#              labeller = labeller(Cond = Cond.labs, EFS = EFS.labs),
#              scales = 'free_y') + coord_cartesian(ylim=c(-1.2, 0.10))
