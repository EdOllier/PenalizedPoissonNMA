set.seed(123456)
source("CODE_FUNCTION.r")

#########################
### Load example DATA ###
#########################

# DATA_SIM_POISSON_CON is a dataframe (collapse format) with the following columns:
# Y : Number of events
# LTAR : Logarithm of the time at risk
# P : time period
# ID_T : Trial ID
# TTT : Treatment
# ARM : 0 "control arm", 1 "experimental arm"
# SEX : Binary covariates
# AGE : Binary covariates
# HR_XvsY : Treatement contrast variable comparing treatment Y to treatment X
# SEX_HR_XvsY : Interaction between variable SEX and Treatement contrast variable comparing treatment Y to treatment X
# AGE_HR_XvsY : Interaction between variable AGE and Treatement contrast variable comparing treatment Y to treatment X

load("DATA_EXAMPLE.RData")
#DATA_SIM_POISSON_CON = DATA_SIM_POISSON_CON[,-c(6,24:34)]

###############################################################
### Construction of the design matrix for the lasso solver  ###
###############################################################

nt=5 # 5 trials per comparison
N_Time_P = 4 # Follow-up divided in 4 time periods

N_COMP = length(levels(DATA_SIM_POISSON_CON$COMP)) # total number of comparisons
N_TTT = length(levels(DATA_SIM_POISSON_CON$TTT)) # number of treatments
N_TRIAL = rep(nt,times=N_COMP)

# Construction of an indicator variable for each time periods (reference = first time period)
MAT_PER_CON = matrix(0,dim(DATA_SIM_POISSON_CON)[1],N_Time_P-1)
    
for (k in 1:(N_Time_P-1)){
    MAT_PER_CON[,k] = (DATA_SIM_POISSON_CON$P==(k+1))*1
}
    
MAT_PER_CON = as.data.frame(MAT_PER_CON)


ID_COL_HR_CON = 11:14 #column ID of treatement contrast variable
ID_COL_HR_CON_SEX = 15:18 #column ID of interaction between variable SEX and treatement contrast variable
ID_COL_HR_CON_AGE = 19:22 #column ID of interaction between variable AGE and treatement contrast variable

#COMP_LEVELS_CON = levels(DATA_SIM_POISSON_CON$COMP)

# INTERACTION CONTRAST-CONTRAST
INTER_MAT_HR_CON = model.matrix(as.formula(paste("~-1 + (",paste(names(DATA_SIM_POISSON_CON)[ID_COL_HR_CON], collapse = " + "),")^2",sep="")), data = DATA_SIM_POISSON_CON)[,-(1:length(ID_COL_HR_CON))]
INTER_MAT_HR_CON = INTER_MAT_HR_CON[,-which(apply(INTER_MAT_HR_CON,2,sum)==0)]

# INTERACTION COVARIATE-CONTRAST
MAT_INTER_TTT_SEX_CON = DATA_SIM_POISSON_CON[,ID_COL_HR_CON_SEX]
MAT_INTER_TTT_AGE_CON = DATA_SIM_POISSON_CON[,ID_COL_HR_CON_AGE]
    
# VAR AVEC HETERO
NAMES_VAR_CON = c("INTER",names(DATA_SIM_POISSON_CON)[ID_COL_HR_CON])
NAMES_STRATA_CON = c(rep("ID_T", 1+length(names(DATA_SIM_POISSON_CON)[ID_COL_HR_CON]) ))
    
# Design matrix of TIME EPRIOD by Treatement CONTRAST INTERACTION
INTER_MAT_HR_P_CON = model.matrix(as.formula(paste("~-1 + (",paste(names(DATA_SIM_POISSON_CON)[ID_COL_HR_CON], collapse = " + "),")*(",paste(names(MAT_PER_CON), collapse = " + "),")",sep="")), data = cbind(DATA_SIM_POISSON_CON,MAT_PER_CON))[,-(1:(length(ID_COL_HR_CON)+dim(MAT_PER_CON)[2]))]


DESIGN_MAT_UnPL1L2_CON = cbind(
    INTER=rep(1,dim(DATA_SIM_POISSON_CON)[1]),
    DATA_SIM_POISSON_CON[,c(4,ID_COL_HR_CON)],
    MAT_PER_CON,INTER_MAT_HR_CON,
    DATA_SIM_POISSON_CON[,7:8],
    MAT_INTER_TTT_SEX_CON,
    MAT_INTER_TTT_AGE_CON,
    INTER_MAT_HR_P_CON)

############################
### LAUNCH of the SOLVER ###
############################
    
SEL_AdLasso <- MODEL_SELECTION_AdaptLASSO(DESIGN_MAT_UnPL1L2_CON ,
    OBS=DATA_SIM_POISSON_CON$Y ,
    OFFSET = DATA_SIM_POISSON_CON$LTAR,
    VAR_HET = NAMES_VAR_CON,
    STRATA_HET = NAMES_STRATA_CON,
    ID_UnPen = c(3:(2+(N_TTT-1)+(N_Time_P-1))),
    ID_PEN_L1 = (3+(N_TTT-1)+(N_Time_P-1)):dim(DESIGN_MAT_UnPL1L2_CON)[2],
    NORM="Ridge" ,
    lambda_init = 1e-5)

###################################
### Selected model based on BIC ###
###################################

ID_BIC_MIN = which(SEL_AdLasso$BIC==min(SEL_AdLasso$BIC))[1]
COEF_SELECT = SEL_AdLasso$COEF[ID_BIC_MIN,] #Coefficients of the selected model
SD_EST = sqrt(SEL_AdLasso$VAR[ID_BIC_MIN,]) #estimation of between trial variability (SD)


