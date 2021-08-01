library(survsim)
library(lme4)
library(ggplot2)
library(survival)
library(penalized)
library(plyr)

######################################
### Mise en forme pour DATA SHARED ###
######################################

Calise = function(Feat, strata_num)
{
    # block diagonal bg(X1,X2,....,XK)
    Feat = as.matrix(Feat)
    Strat = unique(strata_num)
    nstr = length(Strat)
    n    = nrow(Feat)
    p    = ncol(Feat)
    Fin  = matrix(0, n, p*nstr)
    i0   = 1
    for (k in 1:nstr)
    {
        #browser()
        ID_S = which(strata_num==Strat[k])
        Fin[ID_S, ((k-1)*p+1):(k*p)] = as.matrix(Feat[ID_S,])
    }
    return(Fin)
}

#########################################
### POISSON IPD NETWORK Meta-Analysis ###
#########################################

OPTIM_LAMBDA_HETmult_penalized_L2L1 <- function(DESGN_MAT ,OBS ,OFFSET, VAR_HET , STRATA_HET ,ID_UnPen, ID_PEN_L1, LAMBDA_L1 ,AW_L1=NULL ,NORM="Ridge" , lambda_init = 1e-10 , COEF_INIT_PEN = NULL, COEF_INIT_UPEN = NULL , MaxItPen = MaxItPen, epsilonPen = epsilonPen ){
    
    Iter_Max = 20 # Maximum iterations of the algorithm
    MaxItPen = MaxItPen
    epsilonPen = epsilonPen
    
    MAT_PEN_L1 = DESGN_MAT[,ID_PEN_L1]
    MAT_PEN_L1 = apply(MAT_PEN_L1,2,FUN=function(x){ return((x - mean(x))/sd(x)) })
    
    if (is.null(AW_L1)){ AW_L1 = rep(1,dim(MAT_PEN_L1)[2])}
    
    # FORMULA FOR UNPENALIZED VARIABLES
    DESGN_MATunpen = DESGN_MAT
    DESGN_MATunpen$OFF = OFFSET
    FORM = paste("~",paste(names(DESGN_MATunpen)[ID_UnPen], collapse = " + "),"+ offset(OFF)",sep = " ")
    
    N_HET = length(VAR_HET) # NUMBER OF VARIABLES WITH BETWEEN TRIAL VARIABILITY
    NOBS = dim(DESGN_MAT)[1] # NUMBER OF OBERVATION
    ID_HET = list()
    
    DF = rep(0,N_HET) # INITIALISATION OF DEGREE OF FREEDOM
    DFvect= DF
    S2 = rep(0,N_HET)
    S2vect= S2
    if( length(lambda_init)==N_HET ){ LN = lambda_init }else{ LN = rep(lambda_init[1],N_HET) }
    Lvect = LN
    VN = 1/(LN)
    Vvect = 1/(LN)
    
    # BUILDING BLOCK DIAGONAL MATRIX FOR VARIABLES WITH BETWEEN TRIAL VARIABILITY #
    HET_MAT_PEN = c()
    for (k in 1:N_HET){
        HIER_MAT = Calise(DESGN_MAT[[ VAR_HET[k] ]], DESGN_MAT[[ STRATA_HET[k] ]])
        if(k==1){ID_HET[[k]] = 1:dim(HIER_MAT)[2]}else{ID_HET[[k]] = dim(HET_MAT_PEN)[2] + 1:dim(HIER_MAT)[2]}
        HET_MAT_PEN = cbind( HET_MAT_PEN , HIER_MAT)
    }
    HET_MAT_PEN = as.matrix(HET_MAT_PEN) # DESGIN MATRIX WITH L2 PENALIZED COLUMN
    
    DESGN_MAT_INTERCEPT = as.matrix(cbind(rep(1,dim(DESGN_MATunpen)[1]),DESGN_MATunpen[,ID_UnPen],HET_MAT_PEN,MAT_PEN_L1)) #DESIGN MATRIX FOR COMPUTATION OF THE HESSIAN
    Nvar_UnPEN = length(ID_UnPen)+1
    
    PEN_WEIGHTS = rep(0,dim(HET_MAT_PEN)[2])
    for (k in 1:N_HET){
        PEN_WEIGHTS[ ID_HET[[k]] ] = LN[k]
    }
    ID_COL_NUL_HET = which(apply(HET_MAT_PEN==0,2,prod)==1)
    
    ### UPDATE HET_MAT_PEN with L1 PENALIZED COVARIATES ###
    ID_COL_L2 = 1:dim(HET_MAT_PEN)[2]
    ID_COL_L1 = dim(HET_MAT_PEN)[2] + 1:dim(MAT_PEN_L1)[2]
    HET_MAT_PEN = cbind(HET_MAT_PEN,MAT_PEN_L1)
    PEN_WEIGHTS_L2 = c( PEN_WEIGHTS, rep(0,dim(MAT_PEN_L1)[2]) )
    PEN_WEIGHTS_L1 = c( 0*PEN_WEIGHTS, AW_L1*LAMBDA_L1 )
    
    HET_MAT_PEN = as.matrix(HET_MAT_PEN) # FINAL DESIGN MATRIX FOR PENALIZED (L1 AND/OR L2) CORVARIATES

    ### INITIALISATION OF PARAMETERS ###
    startbeta = COEF_INIT_PEN # Penalized parameters
    startgamma = COEF_INIT_UPEN # Non Penalized parameters
    
    ### LAUNCH SOLVER ###
    if ( (!is.null(COEF_INIT_PEN)) &  (!is.null(COEF_INIT_UPEN)) ){
    DSP = penalized(response = OBS, unpenalized=as.formula(FORM), penalized = HET_MAT_PEN , data=DESGN_MATunpen, lambda1=PEN_WEIGHTS_L1, lambda2=PEN_WEIGHTS_L2, model="poisson",standardize=FALSE,trace = FALSE, startbeta =startbeta , startgamma = startgamma,maxiter=MaxItPen, epsilon =epsilonPen)
    }else{
       DSP = penalized(response = OBS, unpenalized=as.formula(FORM), penalized = HET_MAT_PEN , data=DESGN_MATunpen, lambda1=PEN_WEIGHTS_L1, lambda2=PEN_WEIGHTS_L2, model="poisson",standardize=FALSE,trace = FALSE,maxiter=MaxItPen, epsilon =epsilonPen)
    }

    for (i in 1:Iter_Max){
        
        ### HESSIAN MATRIX COMPUTATION ###
        ACTIVE_SET_ID = c(1:Nvar_UnPEN, Nvar_UnPEN + ID_COL_L2, Nvar_UnPEN + length(ID_COL_L2) + which(coefficients(DSP, "penalized")[ID_COL_L1]!=0))
        W= sparseMatrix(i=1:NOBS, j=1:NOBS, x = exp(linear.predictors(DSP)))
        HESS = - (t(DESGN_MAT_INTERCEPT)%*%W%*%DESGN_MAT_INTERCEPT)
        H = try(solve(HESS[ACTIVE_SET_ID,ACTIVE_SET_ID] - diag(c(rep(0,Nvar_UnPEN),PEN_WEIGHTS_L2))[ACTIVE_SET_ID,ACTIVE_SET_ID]  )%*%HESS[ACTIVE_SET_ID,ACTIVE_SET_ID],silent=TRUE)
        if(class(H)=="try-error"){break}
        
        ### UPDATE Degree of freedom and L2 penalization parameters value ###
        for (k in 1:N_HET){
            DF[k] = sum(diag(H[Nvar_UnPEN+ID_HET[[k]],Nvar_UnPEN+ID_HET[[k]]]))
            VN[k] = sum(coefficients(DSP, "penalized")[ID_HET[[k]]]^2)/DF[k]
            LN[k] =  1/(VN[k])
            S2[k] = sum(coefficients(DSP, "penalized")[ID_HET[[k]]]^2)
        }
        
        Lvect = rbind(Lvect,LN)
        Vvect = rbind(Vvect,VN)
        DFvect = rbind(DFvect,DF)
        S2vect = rbind(S2vect,S2)
        
        PEN_WEIGHTS_L2 = rep(0,dim(HET_MAT_PEN)[2])
        
        for (k in 1:N_HET){
            PEN_WEIGHTS_L2[ ID_HET[[k]] ] = LN[k]
        }
        
        ### LAUNCH SOLVER WITH NEW  PENALIZATION PARAMETERS ###
        DSP = penalized(response = OBS, unpenalized=as.formula(FORM), penalized = HET_MAT_PEN , data=DESGN_MATunpen, lambda1=PEN_WEIGHTS_L1, lambda2=PEN_WEIGHTS_L2, model="poisson",standardize=FALSE,startbeta=coefficients(DSP, "penalized"),startgamma=coefficients(DSP, "unpenalized"),trace = FALSE,maxiter=MaxItPen, epsilon =epsilonPen)
        
        #print(i)
        
        ### CHECK CONVERGENCE ###
        if( prod ( ((abs(Lvect[i,] - Lvect[i+1,])/Lvect[i,]) <1e-3) | ( (abs(Vvect[i,] - Vvect[i+1,]) )<1e-6 )  ) ){break}
        
    }
    
    ACTIVE_SET_ID = c(1:Nvar_UnPEN, Nvar_UnPEN + ID_COL_L2, Nvar_UnPEN + length(ID_COL_L2) + which(coefficients(DSP, "penalized")[ID_COL_L1]!=0))
    W= sparseMatrix(i=1:NOBS, j=1:NOBS, x = exp(linear.predictors(DSP)))
    HESS = - (t(DESGN_MAT_INTERCEPT)%*%W%*%DESGN_MAT_INTERCEPT)
    H = try(solve(HESS[ACTIVE_SET_ID,ACTIVE_SET_ID] - diag(c(rep(0,Nvar_UnPEN),PEN_WEIGHTS_L2))[ACTIVE_SET_ID,ACTIVE_SET_ID]  )%*%HESS[ACTIVE_SET_ID,ACTIVE_SET_ID],silent=TRUE)
    DF = sum(diag(H))
    
    return(list(HR_EST = coefficients(DSP,"all"),
    COEF_PEN_EST = coefficients(DSP,"penalized"),
    COEF_UPEN_EST = coefficients(DSP,"unpenalized"),
    BIC = -2*loglik(DSP) + DF*log(NOBS),
    DF = DF,
    LIK = -2*loglik(DSP),
    ID_L1 = Nvar_UnPEN + ID_COL_L1 ,
    Vvect=Vvect,Lvect=Lvect,DFvect=DFvect))
    
}


MODEL_SELECTION_AdaptLASSO <- function(DESGN_MAT ,OBS ,OFFSET, VAR_HET , STRATA_HET ,ID_UnPen, ID_PEN_L1 ,NORM="Ridge" , lambda_init = 1e-5 , COEF_INIT_PEN = NULL, COEF_INIT_UPEN = NULL){
    
    MaxItPen = 1e5
    epsilonPen = 1e-5
    
    UPEN_EST <- OPTIM_LAMBDA_HETmult_penalized_L2L1(DESGN_MAT, OBS, OFFSET, VAR_HET, STRATA_HET, ID_UnPen, ID_PEN_L1, 1e-5, AW_L1=NULL, NORM="Ridge" , lambda_init = lambda_init , COEF_INIT_PEN = NULL, COEF_INIT_UPEN = NULL, MaxItPen = MaxItPen, epsilonPen = epsilonPen)
    
    
    AW_L1 = 1/abs(UPEN_EST$HR_EST[UPEN_EST$ID_L1])

    L1_SEQ = exp(c(500,200,50,20,10,8,6,5,4.5,seq(4,-5,-0.05)))
    
    BIC_SEQ = c()
    LIK_SEQ = c()
    DF_SEQ = c()
    COEF_SEQ = c()
    COEF_SEQ_PEN = c()
    VAR_SEQ = c()
    
    for (L1 in L1_SEQ){
        
        PEN_EST <- OPTIM_LAMBDA_HETmult_penalized_L2L1(DESGN_MAT, OBS, OFFSET, VAR_HET, STRATA_HET, ID_UnPen, ID_PEN_L1, L1 , AW_L1, NORM="Ridge" , lambda_init = lambda_init , COEF_INIT_PEN = COEF_INIT_PEN, COEF_INIT_UPEN = COEF_INIT_UPEN, MaxItPen = MaxItPen, epsilonPen = epsilonPen)
        
        
        COEF_INIT_UPEN = PEN_EST$COEF_UPEN_EST
        COEF_INIT_PEN = PEN_EST$COEF_PEN_EST
        lambda_init = PEN_EST$Lvect[dim(PEN_EST$Lvect)[1],]
        
        AW_R =AW_L1*0
        AW_R[which(PEN_EST$HR_EST[PEN_EST$ID_L1]==0)] = 1e10
        
        COEF_SEQ_PEN = rbind(COEF_SEQ_PEN,PEN_EST$HR_EST)
        
        if(L1==L1_SEQ[1]){AW_Rold = AW_R}
        
        if ((prod(AW_R==AW_Rold)!=1)|(L1==L1_SEQ[1])){
            
            if (sum(AW_R)==0){
                RELAX_EST <- UPEN_EST
            }else{
                RELAX_EST <- OPTIM_LAMBDA_HETmult_penalized_L2L1(DESGN_MAT, OBS, OFFSET, VAR_HET, STRATA_HET, ID_UnPen, ID_PEN_L1, 1 , AW_R, NORM="Ridge" , lambda_init = lambda_init, COEF_INIT_PEN = COEF_INIT_PEN, COEF_INIT_UPEN = COEF_INIT_UPEN , MaxItPen = MaxItPen, epsilonPen = epsilonPen)
            }
            
        
            BIC_SEQ = c(BIC_SEQ,RELAX_EST$BIC)
            LIK_SEQ = c(LIK_SEQ,RELAX_EST$LIK)
            DF_SEQ = c(DF_SEQ,RELAX_EST$DF)
            COEF_SEQ = rbind(COEF_SEQ,RELAX_EST$HR_EST)
            VAR_SEQ = rbind(VAR_SEQ,RELAX_EST$Vvect[dim(RELAX_EST$Vvect)[1],])
            AW_Rold = AW_R
            
        }
        
    }
    
    return(list(BIC=BIC_SEQ, LIK = LIK_SEQ, DF =DF_SEQ , COEF=COEF_SEQ , VAR = VAR_SEQ, COEF_PEN=COEF_SEQ_PEN,L1_SEQ=L1_SEQ) )
    
}

