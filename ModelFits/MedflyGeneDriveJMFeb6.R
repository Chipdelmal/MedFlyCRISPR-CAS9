#########################################################################
## Population suppressing gene drive system in D. melanogaster         ##
## Model fitting code: John Marshall (john.marshall@berkeley.edu)      ##
#########################################################################

# Clear all stored parameters:
rm(list=ls())

# Load required packages:
library(car)
library(truncnorm)
library(ggplot2)
library(coda)
library(grid)

#########################################################################
## Load the data:                                                      ##
#########################################################################

## Read in the data from drive experiments:
setwd("/Users/chipdelmal/Documents/Github/MedFlyCRISPR-CAS9/ModelFits/")

Data <- read.csv("HomingTraLocusGenerationalData.csv")
Data
head(Data)

gens <- Data$Generation

A_R1_M_Pos <- Data$A_R1_M_Pos; A_R1_M_Neg <- Data$A_R1_M_Neg; A_R1_M_All <- Data$A_R1_M_All;
A_R1_F_Pos <- Data$A_R1_F_Pos; A_R1_F_Neg <- Data$A_R1_F_Neg; A_R1_F_All <- Data$A_R1_F_All;
A_R2_M_Pos <- Data$A_R2_M_Pos; A_R2_M_Neg <- Data$A_R2_M_Neg; A_R2_M_All <- Data$A_R2_M_All;
A_R2_F_Pos <- Data$A_R2_F_Pos; A_R2_F_Neg <- Data$A_R2_F_Neg; A_R2_F_All <- Data$A_R2_F_All;
A_R3_M_Pos <- Data$A_R3_M_Pos; A_R3_M_Neg <- Data$A_R3_M_Neg; A_R3_M_All <- Data$A_R3_M_All;
A_R3_F_Pos <- Data$A_R3_F_Pos; A_R3_F_Neg <- Data$A_R3_F_Neg; A_R3_F_All <- Data$A_R3_F_All;
A_R4_M_Pos <- Data$A_R4_M_Pos; A_R4_M_Neg <- Data$A_R4_M_Neg; A_R4_M_All <- Data$A_R4_M_All;
A_R4_F_Pos <- Data$A_R4_F_Pos; A_R4_F_Neg <- Data$A_R4_F_Neg; A_R4_F_All <- Data$A_R4_F_All;
A_R5_M_Pos <- Data$A_R5_M_Pos; A_R5_M_Neg <- Data$A_R5_M_Neg; A_R5_M_All <- Data$A_R5_M_All;
A_R5_F_Pos <- Data$A_R5_F_Pos; A_R5_F_Neg <- Data$A_R5_F_Neg; A_R5_F_All <- Data$A_R5_F_All;

B_R1_M_Pos <- Data$B_R1_M_Pos; B_R1_M_Neg <- Data$B_R1_M_Neg; B_R1_M_All <- Data$B_R1_M_All;
B_R1_F_Pos <- Data$B_R1_F_Pos; B_R1_F_Neg <- Data$B_R1_F_Neg; B_R1_F_All <- Data$B_R1_F_All;
B_R2_M_Pos <- Data$B_R2_M_Pos; B_R2_M_Neg <- Data$B_R2_M_Neg; B_R2_M_All <- Data$B_R2_M_All;
B_R2_F_Pos <- Data$B_R2_F_Pos; B_R2_F_Neg <- Data$B_R2_F_Neg; B_R2_F_All <- Data$B_R2_F_All;
B_R3_M_Pos <- Data$B_R3_M_Pos; B_R3_M_Neg <- Data$B_R3_M_Neg; B_R3_M_All <- Data$B_R3_M_All;
B_R3_F_Pos <- Data$B_R3_F_Pos; B_R3_F_Neg <- Data$B_R3_F_Neg; B_R3_F_All <- Data$B_R3_F_All;
B_R4_M_Pos <- Data$B_R4_M_Pos; B_R4_M_Neg <- Data$B_R4_M_Neg; B_R4_M_All <- Data$B_R4_M_All;
B_R4_F_Pos <- Data$B_R4_F_Pos; B_R4_F_Neg <- Data$B_R4_F_Neg; B_R4_F_All <- Data$B_R4_F_All;
B_R5_M_Pos <- Data$B_R5_M_Pos; B_R5_M_Neg <- Data$B_R5_M_Neg; B_R5_M_All <- Data$B_R5_M_All;
B_R5_F_Pos <- Data$B_R5_F_Pos; B_R5_F_Neg <- Data$B_R5_F_Neg; B_R5_F_All <- Data$B_R5_F_All;

C_R1_M_Pos <- Data$C_R1_M_Pos; C_R1_M_Neg <- Data$C_R1_M_Neg; C_R1_M_All <- Data$C_R1_M_All;
C_R1_F_Pos <- Data$C_R1_F_Pos; C_R1_F_Neg <- Data$C_R1_F_Neg; C_R1_F_All <- Data$C_R1_F_All;
C_R2_M_Pos <- Data$C_R2_M_Pos; C_R2_M_Neg <- Data$C_R2_M_Neg; C_R2_M_All <- Data$C_R2_M_All;
C_R2_F_Pos <- Data$C_R2_F_Pos; C_R2_F_Neg <- Data$C_R2_F_Neg; C_R2_F_All <- Data$C_R2_F_All;
C_R3_M_Pos <- Data$C_R3_M_Pos; C_R3_M_Neg <- Data$C_R3_M_Neg; C_R3_M_All <- Data$C_R3_M_All;
C_R3_F_Pos <- Data$C_R3_F_Pos; C_R3_F_Neg <- Data$C_R3_F_Neg; C_R3_F_All <- Data$C_R3_F_All;
C_R4_M_Pos <- Data$C_R4_M_Pos; C_R4_M_Neg <- Data$C_R4_M_Neg; C_R4_M_All <- Data$C_R4_M_All;
C_R4_F_Pos <- Data$C_R4_F_Pos; C_R4_F_Neg <- Data$C_R4_F_Neg; C_R4_F_All <- Data$C_R4_F_All;
C_R5_M_Pos <- Data$C_R5_M_Pos; C_R5_M_Neg <- Data$C_R5_M_Neg; C_R5_M_All <- Data$C_R5_M_All;
C_R5_F_Pos <- Data$C_R5_F_Pos; C_R5_F_Neg <- Data$C_R5_F_Neg; C_R5_F_All <- Data$C_R5_F_All;

D_R1_M_Pos <- Data$D_R1_M_Pos; D_R1_M_Neg <- Data$D_R1_M_Neg; D_R1_M_All <- Data$D_R1_M_All;
D_R1_F_Pos <- Data$D_R1_F_Pos; D_R1_F_Neg <- Data$D_R1_F_Neg; D_R1_F_All <- Data$D_R1_F_All;
D_R2_M_Pos <- Data$D_R2_M_Pos; D_R2_M_Neg <- Data$D_R2_M_Neg; D_R2_M_All <- Data$D_R2_M_All;
D_R2_F_Pos <- Data$D_R2_F_Pos; D_R2_F_Neg <- Data$D_R2_F_Neg; D_R2_F_All <- Data$D_R2_F_All;
D_R3_M_Pos <- Data$D_R3_M_Pos; D_R3_M_Neg <- Data$D_R3_M_Neg; D_R3_M_All <- Data$D_R3_M_All;
D_R3_F_Pos <- Data$D_R3_F_Pos; D_R3_F_Neg <- Data$D_R3_F_Neg; D_R3_F_All <- Data$D_R3_F_All;
D_R4_M_Pos <- Data$D_R4_M_Pos; D_R4_M_Neg <- Data$D_R4_M_Neg; D_R4_M_All <- Data$D_R4_M_All;
D_R4_F_Pos <- Data$D_R4_F_Pos; D_R4_F_Neg <- Data$D_R4_F_Neg; D_R4_F_All <- Data$D_R4_F_All;
D_R5_M_Pos <- Data$D_R5_M_Pos; D_R5_M_Neg <- Data$D_R5_M_Neg; D_R5_M_All <- Data$D_R5_M_All;
D_R5_F_Pos <- Data$D_R5_F_Pos; D_R5_F_Neg <- Data$D_R5_F_Neg; D_R5_F_All <- Data$D_R5_F_All;

#########################################################################
## Set up the model:                                                   ##
#########################################################################

## Likelihood calculation for model spread of homing allele through a population:

# For a WD heterozygote, the gametes produced would be:
#  * D at a proportion of (1+e)/2, where "e" is the homing efficiency
#  * In-frame resistant allele, R, present in the gametes at proportion rhoR/2
#  * Out-of-frame resistant allele, B, present in the gametes at proportion rhoB/2
#  * Wild-type gametes, W, in the gametes with proportion (1-e-rhoR-rhoB)/2

# The ratio of in-frame resistant alleles (R) to out-of-frame alleles (B) should be
# ~1:2, i.e. rhoR ~= 0.5*rhoB.

# Expected fitness costs are:
# XY individuals:
#  * WW: fitness = 1 (by definition, all other fitnesses measured relative to this)
#  * WD, WR, WB, DD, DR, DB, RR, RB, BB: fitness = 1 (don't care about tra)

# XX individuals:
#  * WW: fitness = 1 (by definition, all other fitnesses measured relative to this)
#  * RR: fitness = 1 (R has no fitness cost because it's an in-frame indel)
#  * WR, WB, DR, RB: fitness = 1 (only need one copy of tra)
#  * WD, DD, DB, BB: fitness = 0 (sterile male)

logLike_HomingMod <- function(prHome, prCleave, prR, sR, Male_Pos_Data, Male_Neg_Data, Male_Total_Data,
                             Female_Pos_Data, Female_Neg_Data, Female_Total_Data) {

  Total_Data <- Male_Total_Data + Female_Total_Data

  rhoR <- prCleave*(1-prHome)*prR
  rhoB <- prCleave*(1-prHome)*(1-prR)
  eD   <- prCleave*prHome

  WW_XX <- (Female_Neg_Data[1] / Total_Data[1])
  WD_XX <- (Female_Pos_Data[1] / Total_Data[1])
  WR_XX <- 0
  WB_XX <- 0
  DD_XX <- 0
  DR_XX <- 0
  DB_XX <- 0
  RR_XX <- 0
  RB_XX <- 0
  BB_XX <- 0

  WW_XY <- (Male_Neg_Data[1] / Total_Data[1])
  WD_XY <- (Male_Pos_Data[1] / Total_Data[1])
  WR_XY <- 0
  WB_XY <- 0
  DD_XY <- 0
  DR_XY <- 0
  DB_XY <- 0
  RR_XY <- 0
  RB_XY <- 0
  BB_XY <- 0

  Male_Pos_Pred <- (DD_XY + DR_XY + DB_XY + WD_XY) + (DD_XX + DB_XX + WD_XX)
  Male_Neg_Pred <- (WW_XY + WR_XY + WB_XY + RR_XY + RB_XY + BB_XY) + (BB_XX)

  Female_Pos_Pred <- DR_XX
  Female_Neg_Pred <- WR_XX + WB_XX + RR_XX + RB_XX + WW_XX

	for (i in 2:length(Total_Data)) {

    # Calculate allele frequencies among male parents:
	  D_Freq_Males <- (0.5*(1+eD)*WD_XY[i-1] + 1*DD_XY[i-1] + 0.5*DR_XY[i-1] + 0.5*DB_XY[i-1])
	  R_Freq_Males <- (0.5*rhoR*WD_XY[i-1] + 0.5*DR_XY[i-1] + 0.5*WR_XY[i-1] + RR_XY[i-1] + 0.5*RB_XY[i-1])
	  B_Freq_Males <- (0.5*rhoB*WD_XY[i-1] + 0.5*DB_XY[i-1] + 0.5*WB_XY[i-1] + 0.5*RB_XY[i-1] + BB_XY[i-1])
	  W_Freq_Males <- (0.5*(1-eD-rhoR-rhoB)*WD_XY[i-1] + 0.5*WR_XY[i-1] + 0.5*WB_XY[i-1] + WW_XY[i-1])

	  # Calculate allele frequencies among female parents:
	  D_Freq_Females <- (0.5*DR_XX[i-1])
	  R_Freq_Females <- (0.5*DR_XX[i-1] + 0.5*WR_XX[i-1] + RR_XX[i-1] + 0.5*RB_XX[i-1])
	  B_Freq_Females <- (0.5*WB_XX[i-1] + 0.5*RB_XX[i-1])
	  W_Freq_Females <- (0.5*WR_XX[i-1] + 0.5*WB_XX[i-1] + WW_XX[i-1])

	  # Calculate allele frequencies among offspring (unnormalized):
	  WW <- W_Freq_Males * W_Freq_Females
	  WD <- W_Freq_Males * D_Freq_Females + D_Freq_Males * W_Freq_Females
	  WR <- W_Freq_Males * R_Freq_Females + R_Freq_Males * W_Freq_Females
	  WB <- W_Freq_Males * B_Freq_Females + B_Freq_Males * W_Freq_Females
	  DD <- D_Freq_Males * D_Freq_Females
	  DR <- D_Freq_Males * R_Freq_Females + R_Freq_Males * D_Freq_Females
	  DB <- D_Freq_Males * B_Freq_Females + B_Freq_Males * D_Freq_Females
	  RR <- R_Freq_Males * R_Freq_Females
	  RB <- R_Freq_Males * B_Freq_Females + B_Freq_Males * R_Freq_Females
	  BB <- B_Freq_Males * B_Freq_Females

	  NormalizingFactor <- (WW + WD + WR + WB + DB + DD + BB) + (DR + RR + RB)*(1-sR/2)

	  WW_XY[i] <- WW/(2*NormalizingFactor) ; WW_XX[i] <- WW/(2*NormalizingFactor)
	  WD_XY[i] <- WD/(2*NormalizingFactor) ; WD_XX[i] <- WD/(2*NormalizingFactor)
	  WR_XY[i] <- WR/(2*NormalizingFactor) ; WR_XX[i] <- WR/(2*NormalizingFactor)
	  WB_XY[i] <- WB/(2*NormalizingFactor) ; WB_XX[i] <- WB/(2*NormalizingFactor)
	  DD_XY[i] <- DD/(2*NormalizingFactor) ; DD_XX[i] <- DD/(2*NormalizingFactor)
	  DR_XY[i] <- DR/(2*NormalizingFactor) ; DR_XX[i] <- DR*(1-sR)/(2*NormalizingFactor)
	  DB_XY[i] <- DB/(2*NormalizingFactor) ; DB_XX[i] <- DB/(2*NormalizingFactor)
	  RR_XY[i] <- RR/(2*NormalizingFactor) ; RR_XX[i] <- RR*(1-sR)/(2*NormalizingFactor)
	  RB_XY[i] <- RB/(2*NormalizingFactor) ; RB_XX[i] <- RB*(1-sR)/(2*NormalizingFactor)
	  BB_XY[i] <- BB/(2*NormalizingFactor) ; BB_XX[i] <- BB/(2*NormalizingFactor)

	  Male_Pos_Pred[i] <- (WD_XY[i] + DD_XY[i] + DR_XY[i] + DB_XY[i]) +
	                      (DD_XX[i] + DB_XX[i] + WD_XX[i])
	  Male_Neg_Pred[i] <- (WR_XY[i] + WB_XY[i] + RR_XY[i] + RB_XY[i] +
	                       BB_XY[i] + WW_XY[i]) + (BB_XX[i])

	  Female_Pos_Pred[i] <- DR_XX[i]
	  Female_Neg_Pred[i] <- WR_XX[i] + WB_XX[i] + RR_XX[i] + RB_XX[i] + WW_XX[i]
	}

	## Multinomial likelihood calculation:

  Male_Pos_Data <- Male_Pos_Data[2:length(Male_Pos_Data)]
  Male_Neg_Data <- Male_Neg_Data[2:length(Male_Neg_Data)]
  Female_Pos_Data <- Female_Pos_Data[2:length(Female_Pos_Data)]
  Female_Neg_Data <- Female_Neg_Data[2:length(Female_Neg_Data)]
  Total_Data <- Total_Data[2:length(Total_Data)]

	Male_Pos_Pred <- Male_Pos_Pred[2:length(Male_Pos_Pred)]
	Male_Neg_Pred <- Male_Neg_Pred[2:length(Male_Neg_Pred)]
	Female_Pos_Pred <- Female_Pos_Pred[2:length(Female_Pos_Pred)]
	Female_Neg_Pred <- Female_Neg_Pred[2:length(Female_Neg_Pred)]

	logLike <- 0

	for (i in 2:length(Total_Data)) {
	  if (is.na(Male_Pos_Data[i])) {
	    logLike <- logLike + Male_Total_Data[i]*log(max(Male_Pos_Pred[i] + Male_Neg_Pred[i],1e-10)) +
	      Female_Total_Data[i]*log(max(Female_Pos_Pred[i] + Female_Neg_Pred[i],1e-10))
 	  }
	  else {
	    logLike <- logLike + Male_Pos_Data[i]*log(max(Male_Pos_Pred[i],1e-10)) +
	      Male_Neg_Data[i]*log(max(Male_Neg_Pred[i],1e-10)) +
	      Female_Pos_Data[i]*log(max(Female_Pos_Pred[i],1e-10)) +
	      Female_Neg_Data[i]*log(max(Female_Neg_Pred[i],1e-10))
	  }
	}

	logLike
}

## Calculate log likelihood for all experiments:

logLike_AllExpts <- function(prHome, prCleave, prR, sR) {

  logLikeAll <- (

    logLike_HomingMod(prHome, prCleave, prR, sR, A_R1_M_Pos, A_R1_M_Neg, A_R1_M_All, A_R1_F_Pos, A_R1_F_Neg, A_R1_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, A_R2_M_Pos, A_R2_M_Neg, A_R2_M_All, A_R2_F_Pos, A_R2_F_Neg, A_R2_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, A_R3_M_Pos, A_R3_M_Neg, A_R3_M_All, A_R3_F_Pos, A_R3_F_Neg, A_R3_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, A_R4_M_Pos, A_R4_M_Neg, A_R4_M_All, A_R4_F_Pos, A_R4_F_Neg, A_R4_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, A_R5_M_Pos, A_R5_M_Neg, A_R5_M_All, A_R5_F_Pos, A_R5_F_Neg, A_R5_F_All) +

    logLike_HomingMod(prHome, prCleave, prR, sR, B_R1_M_Pos, B_R1_M_Neg, B_R1_M_All, B_R1_F_Pos, B_R1_F_Neg, B_R1_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, B_R2_M_Pos, B_R2_M_Neg, B_R2_M_All, B_R2_F_Pos, B_R2_F_Neg, B_R2_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, B_R3_M_Pos, B_R3_M_Neg, B_R3_M_All, B_R3_F_Pos, B_R3_F_Neg, B_R3_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, B_R4_M_Pos, B_R4_M_Neg, B_R4_M_All, B_R4_F_Pos, B_R4_F_Neg, B_R4_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, B_R5_M_Pos, B_R5_M_Neg, B_R5_M_All, B_R5_F_Pos, B_R5_F_Neg, B_R5_F_All) +

    logLike_HomingMod(prHome, prCleave, prR, sR, C_R1_M_Pos, C_R1_M_Neg, C_R1_M_All, C_R1_F_Pos, C_R1_F_Neg, C_R1_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, C_R2_M_Pos, C_R2_M_Neg, C_R2_M_All, C_R2_F_Pos, C_R2_F_Neg, C_R2_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, C_R3_M_Pos, C_R3_M_Neg, C_R3_M_All, C_R3_F_Pos, C_R3_F_Neg, C_R3_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, C_R4_M_Pos, C_R4_M_Neg, C_R4_M_All, C_R4_F_Pos, C_R4_F_Neg, C_R4_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, C_R5_M_Pos, C_R5_M_Neg, C_R5_M_All, C_R5_F_Pos, C_R5_F_Neg, C_R5_F_All) +

    logLike_HomingMod(prHome, prCleave, prR, sR, D_R1_M_Pos, D_R1_M_Neg, D_R1_M_All, D_R1_F_Pos, D_R1_F_Neg, D_R1_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, D_R2_M_Pos, D_R2_M_Neg, D_R2_M_All, D_R2_F_Pos, D_R2_F_Neg, D_R2_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, D_R3_M_Pos, D_R3_M_Neg, D_R3_M_All, D_R3_F_Pos, D_R3_F_Neg, D_R3_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, D_R4_M_Pos, D_R4_M_Neg, D_R4_M_All, D_R4_F_Pos, D_R4_F_Neg, D_R4_F_All) +
    logLike_HomingMod(prHome, prCleave, prR, sR, D_R5_M_Pos, D_R5_M_Neg, D_R5_M_All, D_R5_F_Pos, D_R5_F_Neg, D_R5_F_All) )

  logLikeAll
}

## Prior function:
logPrior <- function(prHome, prCleave, prR, sR) {
  # Prior on eD:
  # logPrior_eD <- dunif(eD, min = 0, max = 1, log = TRUE)
  logPrior_prHome <- dnorm(prHome, mean = 0.56, sd = 0.16, log = TRUE)
  # Prior on rho:
  logPrior_prCleave <- dunif(prCleave, min = 0, max = 1, log = TRUE)
  # Prior on prR:
  # logPrior_prR <- dunif(prR, min = 0, max = 1, log = TRUE)
  logPrior_prR <- dnorm(prR, mean = 0.36, sd = 0.052, log = TRUE)
  # Prior on sR:
  logPrior_sR <- dunif(sR, min = 0, max = 1, log = TRUE)

  return(logPrior_prHome + logPrior_prCleave + logPrior_prR + logPrior_sR)
}

## Posterior function:
logPosterior <- function(theta) {
  ## Extract parameters from theta vector:
  prHome <- theta[["prHome"]]; prCleave <- theta[["prCleave"]]; prR <- theta[["prR"]]; sR <- theta[["sR"]];

  ## Calculate the log prior:
  logPrior <- logPrior(prHome, prCleave, prR, sR)

  ## Calculate the log likelihood of the data:
  logLike <- logLike_AllExpts(prHome, prCleave, prR, sR)

  ## Calculate the posterior:
  logPosterior <- logPrior + logLike
  return(logPosterior)
}

theta <- c(prHome = 0.30, prCleave = 1, prR = 0.33, sR = 0.1)
logPosterior(theta)

#########################################################################
## MCMC analysis:                                                      ##
#########################################################################

## Metropolis-Hastings algorithm for searching parameter space:

mcmcMH <- function(logPosterior, initTheta, proposalSD, numIterations) {

  # Evaluate the log posterior at initTheta:
  logPosteriorThetaCurrent <- logPosterior(initTheta)

  # Initialise variables:
  thetaCurrent <- initTheta
  samples <- initTheta
  accepted <- 0

  # Run the MCMC algorithm for numIterations interations:
  for (i in 1:numIterations) {
    # Draw a new theta from a Gaussian proposal distribution and
    # assign this to a variable called thetaProposed:
    prHome_Proposed <- rtruncnorm(n = 1, a = 0, b = 1, mean = thetaCurrent[["prHome"]],
                                  sd = proposalSD[["prHome"]])
    prCleave_Proposed <- rnorm(n = 1, mean = thetaCurrent[["prCleave"]], sd = proposalSD[["prCleave"]])
    prR_Proposed <- rtruncnorm(n = 1, a = 0, b = 1, mean = thetaCurrent[["prR"]],
                               sd = proposalSD[["prR"]])
    sR_Proposed <- rnorm(n = 1, mean = thetaCurrent[["sR"]],
                         sd = proposalSD[["sR"]])
    thetaProposed <- c(prHome = prHome_Proposed, prCleave = prCleave_Proposed, prR = prR_Proposed,
                       sR = sR_Proposed)

    # Evaluate the log posterior function at the proposed theta value:
    logPosteriorThetaProposed <- logPosterior(thetaProposed)

    # Compute the Metropolis-Hastings (log) acceptance probability:
    logAcceptance <- logPosteriorThetaProposed - logPosteriorThetaCurrent

    # Use a random number to determine if thetaProposed will be accepted:
    randNum <- runif(n = 1, min = 0, max = 1)

    ## If accepted, update the thetaCurrent vector, etc.:
    if (randNum < exp(logAcceptance)) {
      thetaCurrent <- thetaProposed
      logPosteriorThetaCurrent <- logPosteriorThetaProposed
      accepted <- accepted + 1
    }

    # Add the current theta to the vector of samples:
    samples <- c(samples, thetaCurrent)

    # Print the current state of chain and acceptance rate:
    cat("iteration:", i, "chain:", thetaCurrent, "acceptance rate:", accepted / i, "\n")
  }
  return(samples)
}

# Running the MCMC algorithm to vary the parameters e, sHet and sHom:
mcmcTrace <- mcmcMH(logPosterior = logPosterior, # posterior distribution
                    initTheta = c(prHome = 0.56, prCleave = 1,
                                  prR = 0.33, sR = 0), # intial parameter guess
                                  # prR = 0.25, sR = 0), # intial parameter guess
                    proposalSD = c(prHome = 0.01, prCleave = 0,
                                   prR = 0, sR = 0.01), # standard deviations of
                                   # prR = 0.01, sR = 0), # standard deviations of
                    # parameters for Gaussian proposal distribution
                    numIterations = 10000) # number of iterations

trace <- matrix(mcmcTrace, ncol = 4, byrow = T)

# Use the package "coda" to convert the trace into this format:
trace <- mcmc(trace)
plot(trace)
summary(trace)

## Remove the first 1000 iterations to allow for burn-in:
traceBurn <- trace[-(1:1000),]
traceBurn <- mcmc(traceBurn)
plot(traceBurn)
summary(traceBurn)

## Check for autocorrelation:
autocorr.plot(traceBurn)

## Subsamping to account for autocorrelation:
subsample <- 20
traceBurn_eD <- traceBurn[,1]
traceBurn_rho <- traceBurn[,2]
traceBurn_prR <- traceBurn[,3]
traceBurn_sR <- traceBurn[,4]
traceBurnAndThin_eD <- traceBurn_eD[seq(1,length(traceBurn_eD),subsample)]
traceBurnAndThin_rho <- traceBurn_rho[seq(1,length(traceBurn_rho),subsample)]
traceBurnAndThin_prR <- traceBurn_prR[seq(1,length(traceBurn_prR),subsample)]
traceBurnAndThin_sR <- traceBurn_sR[seq(1,length(traceBurn_sR),subsample)]
traceBurnAndThin <- cbind(traceBurnAndThin_eD, traceBurnAndThin_rho, traceBurnAndThin_prR, traceBurnAndThin_sR)
traceBurnAndThin <- mcmc(traceBurnAndThin)
plot(traceBurnAndThin)
summary(traceBurnAndThin)

## Check for autocorrelation:
autocorr.plot(traceBurnAndThin)

## Based on the autocorrelation, we should run the chain for 5*20*100 =
## 10000 iterations beyond the 1000 iterations for burn-in:

mcmcTrace <- mcmcMH(logPosterior = logPosterior, # posterior distribution
                    initTheta = c(prHome = 0.89, prCleave = 0.99,
                                  prR = 0.99, sR = 0), # intial parameter guess
                    # prR = 0.25, sR = 0), # intial parameter guess
                    proposalSD = c(prHome = 0.003, prCleave = 0.0003,
                                   prR = 0.003, sR = 0.003), # standard deviations of
                    # prR = 0.01, sR = 0), # standard deviations of
                    # parameters for Gaussian proposal distribution
                    numIterations = 20000) # number of iterations

trace <- matrix(mcmcTrace, ncol = 4, byrow = T)

# Use the package "coda" to convert the trace into this format:
trace <- mcmc(trace)
plot(trace)
summary(trace)

## Remove the first 1000 iterations to allow for burn-in:
traceBurn <- trace[-(1:1000),]
traceBurn <- mcmc(traceBurn)
plot(traceBurn)
summary(traceBurn)

#########################################################################
## Model selection                                                     ##
#########################################################################

## Calculating the DIC, we have:
logLikeTraceBurn <- rep(0, length(traceBurn[,1]))
for (i in 1:length(traceBurn[,1])) {
  logLikeTraceBurn[i] <- logLike_AllExpts(prHome = traceBurn[[i,1]], prCleave = traceBurn[[i,2]],
                                          prR = traceBurn[[i,3]], sR = traceBurn[[i,4]])
  if (i %% 100 == 0) {
    print(i)
  }
}
pD <- 0.5*var(-2*logLikeTraceBurn) # Effective number of parameters

logLikeMean <- logLike_AllExpts(prHome = mean(traceBurn[,1]),
                                prCleave = mean(traceBurn[,2]),
                                prR = mean(traceBurn[,3]),
                                sR = mean(traceBurn[,4]))

DIC_ModelFitnessCosts <- 2*pD - 2*logLikeMean
DIC_ModelFitnessCosts

#########################################################################
## Now let's plot the results:                                         ##
#########################################################################

## Let's use the following parameter values for the plots:
# prHome_Fitted <- 0.886; prCleave_Fitted <- 0.9999; prR_Fitted <- 0.9943; sR_Fitted <- 0.1268;
prHome_Fitted <- 0.892; prCleave_Fitted <- 1; prR_Fitted <- 0.33; sR_Fitted <- 0;

## Trajectory calculation for model spread of Drosophila homing system
## through a population:

traj_HomingMod <- function(prHome, prCleave, prR, WW_XX_0, WD_XX_0, WW_XY_0, WD_XY_0, numGens, sR) {

  rhoR <- prCleave*(1-prHome)*prR
  rhoB <- prCleave*(1-prHome)*(1-prR)
  eD   <- prCleave*prHome

  WW_XX <- WW_XX_0
  WD_XX <- WD_XX_0
  WR_XX <- 0
  WB_XX <- 0
  DD_XX <- 0
  DR_XX <- 0
  DB_XX <- 0
  RR_XX <- 0
  RB_XX <- 0
  BB_XX <- 0

  WW_XY <- WW_XY_0
  WD_XY <- WD_XY_0
  WR_XY <- 0
  WB_XY <- 0
  DD_XY <- 0
  DR_XY <- 0
  DB_XY <- 0
  RR_XY <- 0
  RB_XY <- 0
  BB_XY <- 0

  Male_Pos_Pred <- (DD_XY + DR_XY + DB_XY + WD_XY) + (DD_XX + DB_XX + WD_XX)
  Male_Neg_Pred <- (WW_XY + WR_XY + WB_XY + RR_XY + RB_XY + BB_XY) + (BB_XX)

  Female_Pos_Pred <- DR_XX
  Female_Neg_Pred <- WR_XX + WB_XX + RR_XX + RB_XX + WW_XX

  Female_Pred <- DR_XX + WR_XX + WB_XX + RR_XX + RB_XX + WW_XX
  D_Allele_Pred <- (0.5*WD_XY + DD_XY + 0.5*DR_XY + 0.5*DB_XY +
                      DD_XX + 0.5*DB_XX + 0.5*WD_XX + 0.5*DR_XX)
  W_Allele_Pred <- (0.5*WD_XY + 0.5*WD_XX + 0.5*WR_XY + 0.5*WB_XY +
                      WW_XY + 0.5*WR_XX + 0.5*WB_XX + WW_XX)
  R_Allele_Pred <- (0.5*DR_XY + 0.5*WR_XY + RR_XY + 0.5*RB_XY +
                      0.5*DR_XX + 0.5*WR_XX + RR_XX + 0.5*RB_XX)
  B_Allele_Pred <- (0.5*DB_XY + 0.5*DB_XX + 0.5*WB_XY + 0.5*RB_XY +
                      BB_XY + BB_XX + 0.5*WB_XX + 0.5*RB_XX)

  for (i in 2:numGens) {

    # Calculate allele frequencies among male parents:
    D_Freq_Males <- (0.5*(1+eD)*WD_XY[i-1] + 1*DD_XY[i-1] + 0.5*DR_XY[i-1] + 0.5*DB_XY[i-1])
    R_Freq_Males <- (0.5*rhoR*WD_XY[i-1] + 0.5*DR_XY[i-1] + 0.5*WR_XY[i-1] + RR_XY[i-1] + 0.5*RB_XY[i-1])
    B_Freq_Males <- (0.5*rhoB*WD_XY[i-1] + 0.5*DB_XY[i-1] + 0.5*WB_XY[i-1] + 0.5*RB_XY[i-1] + BB_XY[i-1])
    W_Freq_Males <- (0.5*(1-eD-rhoR-rhoB)*WD_XY[i-1] + 0.5*WR_XY[i-1] + 0.5*WB_XY[i-1] + WW_XY[i-1])

    # Calculate allele frequencies among female parents:
    D_Freq_Females <- (0.5*DR_XX[i-1])
    R_Freq_Females <- (0.5*DR_XX[i-1] + 0.5*WR_XX[i-1] + RR_XX[i-1] + 0.5*RB_XX[i-1])
    B_Freq_Females <- (0.5*WB_XX[i-1] + 0.5*RB_XX[i-1])
    W_Freq_Females <- (0.5*WR_XX[i-1] + 0.5*WB_XX[i-1] + WW_XX[i-1])

    # Calculate allele frequencies among offspring (unnormalized):
    WW <- W_Freq_Males * W_Freq_Females
    WD <- W_Freq_Males * D_Freq_Females + D_Freq_Males * W_Freq_Females
    WR <- W_Freq_Males * R_Freq_Females + R_Freq_Males * W_Freq_Females
    WB <- W_Freq_Males * B_Freq_Females + B_Freq_Males * W_Freq_Females
    DD <- D_Freq_Males * D_Freq_Females
    DR <- D_Freq_Males * R_Freq_Females + R_Freq_Males * D_Freq_Females
    DB <- D_Freq_Males * B_Freq_Females + B_Freq_Males * D_Freq_Females
    RR <- R_Freq_Males * R_Freq_Females
    RB <- R_Freq_Males * B_Freq_Females + B_Freq_Males * R_Freq_Females
    BB <- B_Freq_Males * B_Freq_Females

    NormalizingFactor <- (WW + WD + WR + WB + DB + DD + BB) + (DR + RR + RB)*(1-sR/2)

    WW_XY[i] <- WW/(2*NormalizingFactor) ; WW_XX[i] <- WW/(2*NormalizingFactor)
    WD_XY[i] <- WD/(2*NormalizingFactor) ; WD_XX[i] <- WD/(2*NormalizingFactor)
    WR_XY[i] <- WR/(2*NormalizingFactor) ; WR_XX[i] <- WR/(2*NormalizingFactor)
    WB_XY[i] <- WB/(2*NormalizingFactor) ; WB_XX[i] <- WB/(2*NormalizingFactor)
    DD_XY[i] <- DD/(2*NormalizingFactor) ; DD_XX[i] <- DD/(2*NormalizingFactor)
    DR_XY[i] <- DR/(2*NormalizingFactor) ; DR_XX[i] <- DR*(1-sR)/(2*NormalizingFactor)
    DB_XY[i] <- DB/(2*NormalizingFactor) ; DB_XX[i] <- DB/(2*NormalizingFactor)
    RR_XY[i] <- RR/(2*NormalizingFactor) ; RR_XX[i] <- RR*(1-sR)/(2*NormalizingFactor)
    RB_XY[i] <- RB/(2*NormalizingFactor) ; RB_XX[i] <- RB*(1-sR)/(2*NormalizingFactor)
    BB_XY[i] <- BB/(2*NormalizingFactor) ; BB_XX[i] <- BB/(2*NormalizingFactor)

    Male_Pos_Pred[i] <- ((WD_XY[i] + DD_XY[i] + DR_XY[i] + DB_XY[i]) +
                          (DD_XX[i] + DB_XX[i] + WD_XX[i]))
    Male_Neg_Pred[i] <- ((WR_XY[i] + WB_XY[i] + RR_XY[i] + RB_XY[i] +
                           BB_XY[i] + WW_XY[i]) + (BB_XX[i]))

    Female_Pos_Pred[i] <- DR_XX[i]
    Female_Neg_Pred[i] <- (WR_XX[i] + WB_XX[i] + RR_XX[i] + RB_XX[i] + WW_XX[i])

    Female_Pred[i] <- DR_XX[i] + WR_XX[i] + WB_XX[i] + RR_XX[i] + RB_XX[i] + WW_XX[i]
    D_Allele_Pred[i] <- (0.5*WD_XY[i] + DD_XY[i] + 0.5*DR_XY[i] + 0.5*DB_XY[i] +
                           DD_XX[i] + 0.5*DB_XX[i] + 0.5*WD_XX[i] + 0.5*DR_XX[i])
    W_Allele_Pred[i] <- (0.5*WD_XY[i] + 0.5*WD_XX[i] + 0.5*WR_XY[i] + 0.5*WB_XY[i] +
                           WW_XY[i] + 0.5*WR_XX[i] + 0.5*WB_XX[i] + WW_XX[i])
    R_Allele_Pred[i] <- (0.5*DR_XY[i] + 0.5*WR_XY[i] + RR_XY[i] + 0.5*RB_XY[i] +
                           0.5*DR_XX[i] + 0.5*WR_XX[i] + RR_XX[i] + 0.5*RB_XX[i])
    B_Allele_Pred[i] <- (0.5*DB_XY[i] + 0.5*DB_XX[i] + 0.5*WB_XY[i] + 0.5*RB_XY[i] +
                           BB_XY[i] + BB_XX[i] + 0.5*WB_XX[i] + 0.5*RB_XX[i])
  }

  cbind(Male_Pos_Pred, Female_Pos_Pred, Female_Pred, D_Allele_Pred,
        W_Allele_Pred, R_Allele_Pred, B_Allele_Pred)
}

## Prepare the observed data for plotting (M & F frequencies together):

A_R1_M_Obs <- A_R1_M_Pos / (A_R1_M_All + A_R1_F_All); A_R1_F_Obs <- A_R1_F_Pos / (A_R1_M_All + A_R1_F_All);
A_R2_M_Obs <- A_R2_M_Pos / (A_R2_M_All + A_R2_F_All); A_R2_F_Obs <- A_R2_F_Pos / (A_R2_M_All + A_R2_F_All);
A_R3_M_Obs <- A_R3_M_Pos / (A_R3_M_All + A_R3_F_All); A_R3_F_Obs <- A_R3_F_Pos / (A_R3_M_All + A_R3_F_All);
A_R4_M_Obs <- A_R4_M_Pos / (A_R4_M_All + A_R4_F_All); A_R4_F_Obs <- A_R4_F_Pos / (A_R4_M_All + A_R4_F_All);
A_R5_M_Obs <- A_R5_M_Pos / (A_R5_M_All + A_R5_F_All); A_R5_F_Obs <- A_R5_F_Pos / (A_R5_M_All + A_R5_F_All);

B_R1_M_Obs <- B_R1_M_Pos / (B_R1_M_All + B_R1_F_All); B_R1_F_Obs <- B_R1_F_Pos / (B_R1_M_All + B_R1_F_All);
B_R2_M_Obs <- B_R2_M_Pos / (B_R2_M_All + B_R2_F_All); B_R2_F_Obs <- B_R2_F_Pos / (B_R2_M_All + B_R2_F_All);
B_R3_M_Obs <- B_R3_M_Pos / (B_R3_M_All + B_R3_F_All); B_R3_F_Obs <- B_R3_F_Pos / (B_R3_M_All + B_R3_F_All);
B_R4_M_Obs <- B_R4_M_Pos / (B_R4_M_All + B_R4_F_All); B_R4_F_Obs <- B_R4_F_Pos / (B_R4_M_All + B_R4_F_All);
B_R5_M_Obs <- B_R5_M_Pos / (B_R5_M_All + B_R5_F_All); B_R5_F_Obs <- B_R5_F_Pos / (B_R5_M_All + B_R5_F_All);

C_R1_M_Obs <- C_R1_M_Pos / (C_R1_M_All + C_R1_F_All); C_R1_F_Obs <- C_R1_F_Pos / (C_R1_M_All + C_R1_F_All);
C_R2_M_Obs <- C_R2_M_Pos / (C_R2_M_All + C_R2_F_All); C_R2_F_Obs <- C_R2_F_Pos / (C_R2_M_All + C_R2_F_All);
C_R3_M_Obs <- C_R3_M_Pos / (C_R3_M_All + C_R3_F_All); C_R3_F_Obs <- C_R3_F_Pos / (C_R3_M_All + C_R3_F_All);
C_R4_M_Obs <- C_R4_M_Pos / (C_R4_M_All + C_R4_F_All); C_R4_F_Obs <- C_R4_F_Pos / (C_R4_M_All + C_R4_F_All);
C_R5_M_Obs <- C_R5_M_Pos / (C_R5_M_All + C_R5_F_All); C_R5_F_Obs <- C_R5_F_Pos / (C_R5_M_All + C_R5_F_All);

D_R1_M_Obs <- D_R1_M_Pos / (D_R1_M_All + D_R1_F_All); D_R1_F_Obs <- D_R1_F_Pos / (D_R1_M_All + D_R1_F_All);
D_R2_M_Obs <- D_R2_M_Pos / (D_R2_M_All + D_R2_F_All); D_R2_F_Obs <- D_R2_F_Pos / (D_R2_M_All + D_R2_F_All);
D_R3_M_Obs <- D_R3_M_Pos / (D_R3_M_All + D_R3_F_All); D_R3_F_Obs <- D_R3_F_Pos / (D_R3_M_All + D_R3_F_All);
D_R4_M_Obs <- D_R4_M_Pos / (D_R4_M_All + D_R4_F_All); D_R4_F_Obs <- D_R4_F_Pos / (D_R4_M_All + D_R4_F_All);
D_R5_M_Obs <- D_R5_M_Pos / (D_R5_M_All + D_R5_F_All); D_R5_F_Obs <- D_R5_F_Pos / (D_R5_M_All + D_R5_F_All);

# Filling in NAs:

A_R1_M_Obs[2] <- (A_R1_M_Obs[1] + A_R1_M_Obs[3])/2 ; A_R1_F_Obs[2] <- (A_R1_F_Obs[1] + A_R1_F_Obs[3])/2
A_R2_M_Obs[2] <- (A_R2_M_Obs[1] + A_R2_M_Obs[3])/2 ; A_R2_F_Obs[2] <- (A_R2_F_Obs[1] + A_R2_F_Obs[3])/2
A_R3_M_Obs[2] <- (A_R3_M_Obs[1] + A_R3_M_Obs[3])/2 ; A_R3_F_Obs[2] <- (A_R3_F_Obs[1] + A_R3_F_Obs[3])/2
A_R4_M_Obs[2] <- (A_R4_M_Obs[1] + A_R4_M_Obs[3])/2 ; A_R4_F_Obs[2] <- (A_R4_F_Obs[1] + A_R4_F_Obs[3])/2

B_R1_M_Obs[2] <- (B_R1_M_Obs[1] + B_R1_M_Obs[3])/2 ; B_R1_F_Obs[2] <- (B_R1_F_Obs[1] + B_R1_F_Obs[3])/2
B_R2_M_Obs[2] <- (B_R2_M_Obs[1] + B_R2_M_Obs[3])/2 ; B_R2_F_Obs[2] <- (B_R2_F_Obs[1] + B_R2_F_Obs[3])/2
B_R3_M_Obs[2] <- (B_R3_M_Obs[1] + B_R3_M_Obs[3])/2 ; B_R3_F_Obs[2] <- (B_R3_F_Obs[1] + B_R3_F_Obs[3])/2
B_R4_M_Obs[2] <- (B_R4_M_Obs[1] + B_R4_M_Obs[3])/2 ; B_R4_F_Obs[2] <- (B_R4_F_Obs[1] + B_R4_F_Obs[3])/2

C_R1_M_Obs[2] <- (C_R1_M_Obs[1] + C_R1_M_Obs[3])/2 ; C_R1_F_Obs[2] <- (C_R1_F_Obs[1] + C_R1_F_Obs[3])/2
C_R2_M_Obs[2] <- (C_R2_M_Obs[1] + C_R2_M_Obs[3])/2 ; C_R2_F_Obs[2] <- (C_R2_F_Obs[1] + C_R2_F_Obs[3])/2
C_R3_M_Obs[2] <- (C_R3_M_Obs[1] + C_R3_M_Obs[3])/2 ; C_R3_F_Obs[2] <- (C_R3_F_Obs[1] + C_R3_F_Obs[3])/2
C_R4_M_Obs[2] <- (C_R4_M_Obs[1] + C_R4_M_Obs[3])/2 ; C_R4_F_Obs[2] <- (C_R4_F_Obs[1] + C_R4_F_Obs[3])/2

D_R1_M_Obs[2] <- (D_R1_M_Obs[1] + D_R1_M_Obs[3])/2 ; D_R1_F_Obs[2] <- (D_R1_F_Obs[1] + D_R1_F_Obs[3])/2
D_R2_M_Obs[2] <- (D_R2_M_Obs[1] + D_R2_M_Obs[3])/2 ; D_R2_F_Obs[2] <- (D_R2_F_Obs[1] + D_R2_F_Obs[3])/2
D_R3_M_Obs[2] <- (D_R3_M_Obs[1] + D_R3_M_Obs[3])/2 ; D_R3_F_Obs[2] <- (D_R3_F_Obs[1] + D_R3_F_Obs[3])/2
D_R4_M_Obs[2] <- (D_R4_M_Obs[1] + D_R4_M_Obs[3])/2 ; D_R4_F_Obs[2] <- (D_R4_F_Obs[1] + D_R4_F_Obs[3])/2

## Calculate the predicted trajectories for plotting:

A_Pred <- traj_HomingMod(prHome = prHome_Fitted, prCleave = prCleave_Fitted, prR = prR_Fitted,
                         WW_XX_0 = 30/(30 + 33), WD_XX_0 = 0/(30 + 33), WW_XY_0 = 30/(30 + 33), WD_XY_0 =  3/(30 + 33), numGens = 9,
                         sR = sR_Fitted)
A_M_Pred <- as.data.frame(A_Pred)$Male_Pos_Pred
A_F_Pred <- as.data.frame(A_Pred)$Female_Pos_Pred

B_Pred <- traj_HomingMod(prHome = prHome_Fitted, prCleave = prCleave_Fitted, prR = prR_Fitted,
                         WW_XX_0 = 30/(30 + 45), WD_XX_0 = 0/(30 + 45), WW_XY_0 = 30/(30 + 45), WD_XY_0 = 15/(30 + 45), numGens = 9,
                         sR = sR_Fitted)
B_M_Pred <- as.data.frame(B_Pred)$Male_Pos_Pred
B_F_Pred <- as.data.frame(B_Pred)$Female_Pos_Pred

C_Pred <- traj_HomingMod(prHome = prHome_Fitted, prCleave = prCleave_Fitted, prR = prR_Fitted,
                         WW_XX_0 = 30/(30 + 60), WD_XX_0 = 0/(30 + 60), WW_XY_0 = 30/(30 + 60), WD_XY_0 = 30/(30 + 60), numGens = 9,
                         sR = sR_Fitted)
C_M_Pred <- as.data.frame(C_Pred)$Male_Pos_Pred
C_F_Pred <- as.data.frame(C_Pred)$Female_Pos_Pred

D_Pred <- traj_HomingMod(prHome = prHome_Fitted, prCleave = prCleave_Fitted, prR = prR_Fitted,
                         WW_XX_0 = 30/(30 + 30), WD_XX_0 = 0/(30 + 30), WW_XY_0 =  0/(30 + 30), WD_XY_0 = 30/(30 + 30), numGens = 9,
                         sR = sR_Fitted)
D_M_Pred <- as.data.frame(D_Pred)$Male_Pos_Pred
D_F_Pred <- as.data.frame(D_Pred)$Female_Pos_Pred

## Collate observed & predicted data into a data frame:

gens <- 1:9

homingResults <- data.frame(gens,
                            A_M_Pred, A_R1_M_Obs, A_R2_M_Obs, A_R3_M_Obs, A_R4_M_Obs, A_R5_M_Obs,
                            A_F_Pred, A_R1_F_Obs, A_R2_F_Obs, A_R3_F_Obs, A_R4_F_Obs, A_R5_F_Obs,
                            B_M_Pred, B_R1_M_Obs, B_R2_M_Obs, B_R3_M_Obs, B_R4_M_Obs, B_R5_M_Obs,
                            B_F_Pred, B_R1_F_Obs, B_R2_F_Obs, B_R3_F_Obs, B_R4_F_Obs, B_R5_F_Obs,
                            C_M_Pred, C_R1_M_Obs, C_R2_M_Obs, C_R3_M_Obs, C_R4_M_Obs, C_R5_M_Obs,
                            C_F_Pred, C_R1_F_Obs, C_R2_F_Obs, C_R3_F_Obs, C_R4_F_Obs, C_R5_F_Obs,
                            D_M_Pred, D_R1_M_Obs, D_R2_M_Obs, D_R3_M_Obs, D_R4_M_Obs, D_R5_M_Obs,
                            D_F_Pred, D_R1_F_Obs, D_R2_F_Obs, D_R3_F_Obs, D_R4_F_Obs, D_R5_F_Obs)

## Plot the results:

## Release scenario A:

pA_M <- ggplot(homingResults, aes(x = gens, y = A_M_Pred, color = Expt)) +
  geom_line(aes(y = A_M_Pred  , col = "0"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = A_R1_M_Obs, col = "1"), size = 1.2) +
  geom_line(aes(y = A_R2_M_Obs, col = "2"), size = 1.2) +
  geom_line(aes(y = A_R3_M_Obs, col = "3"), size = 1.2) +
  geom_line(aes(y = A_R4_M_Obs, col = "4"), size = 1.2) +
  geom_line(aes(y = A_R5_M_Obs, col = "5"), size = 1.2) +
  labs(x = "Generation", y = "Proportion positive") +
  ggtitle("Release A, Males")

pA_F <- ggplot(homingResults, aes(x = gens, y = A_F_Pred, color = Expt)) +
  geom_line(aes(y = A_F_Pred  , col = "0"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = A_R1_F_Obs, col = "1"), size = 1.2) +
  geom_line(aes(y = A_R2_F_Obs, col = "2"), size = 1.2) +
  geom_line(aes(y = A_R3_F_Obs, col = "3"), size = 1.2) +
  geom_line(aes(y = A_R4_F_Obs, col = "4"), size = 1.2) +
  geom_line(aes(y = A_R5_F_Obs, col = "5"), size = 1.2) +
  labs(x = "Generation", y = "Proportion positive") +
  ggtitle("Release A, Females")

## Release scenario B:

pB_M <- ggplot(homingResults, aes(x = gens, y = B_M_Pred, color = Expt)) +
  geom_line(aes(y = B_M_Pred  , col = "0"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = B_R1_M_Obs, col = "1"), size = 1.2) +
  geom_line(aes(y = B_R2_M_Obs, col = "2"), size = 1.2) +
  geom_line(aes(y = B_R3_M_Obs, col = "3"), size = 1.2) +
  geom_line(aes(y = B_R4_M_Obs, col = "4"), size = 1.2) +
  geom_line(aes(y = B_R5_M_Obs, col = "5"), size = 1.2) +
  labs(x = "Generation", y = "Proportion positive") +
  ggtitle("Release B, Males")

pB_F <- ggplot(homingResults, aes(x = gens, y = B_F_Pred, color = Expt)) +
  geom_line(aes(y = B_F_Pred  , col = "0"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = B_R1_F_Obs, col = "1"), size = 1.2) +
  geom_line(aes(y = B_R2_F_Obs, col = "2"), size = 1.2) +
  geom_line(aes(y = B_R3_F_Obs, col = "3"), size = 1.2) +
  geom_line(aes(y = B_R4_F_Obs, col = "4"), size = 1.2) +
  geom_line(aes(y = B_R5_F_Obs, col = "5"), size = 1.2) +
  labs(x = "Generation", y = "Proportion positive") +
  ggtitle("Release B, Females")

## Release scenario C:

pC_M <- ggplot(homingResults, aes(x = gens, y = C_M_Pred, color = Expt)) +
  geom_line(aes(y = C_M_Pred  , col = "0"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = C_R1_M_Obs, col = "1"), size = 1.2) +
  geom_line(aes(y = C_R2_M_Obs, col = "2"), size = 1.2) +
  geom_line(aes(y = C_R3_M_Obs, col = "3"), size = 1.2) +
  geom_line(aes(y = C_R4_M_Obs, col = "4"), size = 1.2) +
  geom_line(aes(y = C_R5_M_Obs, col = "5"), size = 1.2) +
  labs(x = "Generation", y = "Proportion positive") +
  ggtitle("Release C, Males")

pC_F <- ggplot(homingResults, aes(x = gens, y = C_F_Pred, color = Expt)) +
  geom_line(aes(y = C_F_Pred  , col = "0"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = C_R1_F_Obs, col = "1"), size = 1.2) +
  geom_line(aes(y = C_R2_F_Obs, col = "2"), size = 1.2) +
  geom_line(aes(y = C_R3_F_Obs, col = "3"), size = 1.2) +
  geom_line(aes(y = C_R4_F_Obs, col = "4"), size = 1.2) +
  geom_line(aes(y = C_R5_F_Obs, col = "5"), size = 1.2) +
  labs(x = "Generation", y = "Proportion positive") +
  ggtitle("Release C, Females")

## Release scenario D:

pD_M <- ggplot(homingResults, aes(x = gens, y = D_M_Pred, color = Expt)) +
  geom_line(aes(y = D_M_Pred  , col = "0"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = D_R1_M_Obs, col = "1"), size = 1.2) +
  geom_line(aes(y = D_R2_M_Obs, col = "2"), size = 1.2) +
  geom_line(aes(y = D_R3_M_Obs, col = "3"), size = 1.2) +
  geom_line(aes(y = D_R4_M_Obs, col = "4"), size = 1.2) +
  geom_line(aes(y = D_R5_M_Obs, col = "5"), size = 1.2) +
  labs(x = "Generation", y = "Proportion positive") +
  ggtitle("Release D, Males")

pD_F <- ggplot(homingResults, aes(x = gens, y = D_F_Pred, color = Expt)) +
  geom_line(aes(y = D_F_Pred  , col = "0"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = D_R1_F_Obs, col = "1"), size = 1.2) +
  geom_line(aes(y = D_R2_F_Obs, col = "2"), size = 1.2) +
  geom_line(aes(y = D_R3_F_Obs, col = "3"), size = 1.2) +
  geom_line(aes(y = D_R4_F_Obs, col = "4"), size = 1.2) +
  geom_line(aes(y = D_R5_F_Obs, col = "5"), size = 1.2) +
  labs(x = "Generation", y = "Proportion positive") +
  ggtitle("Release D, Females")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(pA_F, pB_F, pC_F, pD_F, pA_M, pB_M, pC_M, pD_M, cols=2)

#########################################################################
## Now let's plot a figure that shows the averaged results for each    ##
## release alongside model predictions:                                ##
#########################################################################

A_Av_M_Obs <- (A_R1_M_Obs + A_R2_M_Obs + A_R3_M_Obs + A_R4_M_Obs + A_R5_M_Obs) / 5;
A_Av_F_Obs <- (A_R1_F_Obs + A_R2_F_Obs + A_R3_F_Obs + A_R4_F_Obs + A_R5_F_Obs) / 5;

B_Av_M_Obs <- (B_R1_M_Obs + B_R2_M_Obs + B_R3_M_Obs + B_R4_M_Obs + B_R5_M_Obs) / 5;
B_Av_F_Obs <- (B_R1_F_Obs + B_R2_F_Obs + B_R3_F_Obs + B_R4_F_Obs + B_R5_F_Obs) / 5;

C_Av_M_Obs <- (C_R1_M_Obs + C_R2_M_Obs + C_R3_M_Obs + C_R4_M_Obs + C_R5_M_Obs) / 5;
C_Av_F_Obs <- (C_R1_F_Obs + C_R2_F_Obs + C_R3_F_Obs + C_R4_F_Obs + C_R5_F_Obs) / 5;

D_Av_M_Obs <- (D_R1_M_Obs + D_R2_M_Obs + D_R3_M_Obs + D_R4_M_Obs + D_R5_M_Obs) / 5;
D_Av_F_Obs <- (D_R1_F_Obs + D_R2_F_Obs + D_R3_F_Obs + D_R4_F_Obs + D_R5_F_Obs) / 5;

## Collate observed & predicted data into a data frame:

homingResultsFigure <- data.frame(gens,
                                  A_M_Pred, A_Av_M_Obs, A_F_Pred, A_Av_F_Obs,
                                  B_M_Pred, B_Av_M_Obs, B_F_Pred, B_Av_F_Obs,
                                  C_M_Pred, C_Av_M_Obs, C_F_Pred, C_Av_F_Obs,
                                  D_M_Pred, D_Av_M_Obs, D_F_Pred, D_Av_F_Obs)

## Plot the results:

## All release scenarios, males:

gens <- 1:9

Male_Plot <- ggplot(homingResultsFigure, aes(x = gens, y = A_M_Pred, color = Expt)) +
  geom_line(aes(y = A_M_Pred  , col = "1"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = A_Av_M_Obs, col = "1"), size = 1.2) +
  geom_line(aes(y = B_M_Pred, col = "2"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = B_Av_M_Obs, col = "2"), size = 1.2) +
  geom_line(aes(y = C_M_Pred, col = "3"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = C_Av_M_Obs, col = "3"), size = 1.2) +
  geom_line(aes(y = D_M_Pred, col = "4"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = D_Av_M_Obs, col = "4"), size = 1.2) +
  labs(x = "Generation", y = "Proportion positive") +
  ggtitle("All release scenarios, males")

## All release scenarios, females:

Female_Plot <- ggplot(homingResultsFigure, aes(x = gens, y = A_F_Pred, color = Expt)) +
  geom_line(aes(y = A_F_Pred  , col = "1"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = A_Av_F_Obs, col = "1"), size = 1.2) +
  geom_line(aes(y = B_F_Pred, col = "2"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = B_Av_F_Obs, col = "2"), size = 1.2) +
  geom_line(aes(y = C_F_Pred, col = "3"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = C_Av_F_Obs, col = "3"), size = 1.2) +
  geom_line(aes(y = D_F_Pred, col = "4"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = D_Av_F_Obs, col = "4"), size = 1.2) +
  labs(x = "Generation", y = "Proportion positive") +
  ggtitle("All release scenarios, females")

#####################################################################
## Model predictions for fitted model parameters:                  ##
#####################################################################

numGens <- 150
D_Pred <- traj_HomingMod(prHome = prHome_Fitted, prCleave = prCleave_Fitted, prR = prR_Fitted,
                         WW_XX_0 = 30/(30 + 45), WD_XX_0 = 0/(30 + 45), WW_XY_0 = 30/(30 + 45), WD_XY_0 = 15/(30 + 45), numGens = numGens,
                         sR = sR_Fitted)
D_Female_Pred <- as.data.frame(D_Pred)$Female_Pred
D_D_Allele_Pred <- as.data.frame(D_Pred)$D_Allele_Pred
D_W_Allele_Pred <- as.data.frame(D_Pred)$W_Allele_Pred
D_R_Allele_Pred <- as.data.frame(D_Pred)$R_Allele_Pred
D_B_Allele_Pred <- as.data.frame(D_Pred)$B_Allele_Pred

## Model predictions for ideal model parameters:

I_Pred <- traj_HomingMod(prHome = 0.98, prCleave = 1, prR = 0.33,
                         WW_XX_0 = 30/(30 + 45), WD_XX_0 = 0/(30 + 45), WW_XY_0 = 30/(30 + 45), WD_XY_0 = 15/(30 + 45), numGens = numGens,
                         sR = 0)
I_Female_Pred <- as.data.frame(I_Pred)$Female_Pred
I_D_Allele_Pred <- as.data.frame(I_Pred)$D_Allele_Pred
I_W_Allele_Pred <- as.data.frame(I_Pred)$W_Allele_Pred
I_R_Allele_Pred <- as.data.frame(I_Pred)$R_Allele_Pred
I_B_Allele_Pred <- as.data.frame(I_Pred)$B_Allele_Pred

gens <- 1:numGens

DrosophilaModelPred <- data.frame(gens,
                                  D_Female_Pred, D_D_Allele_Pred, D_W_Allele_Pred,
                                  D_R_Allele_Pred, D_B_Allele_Pred,
                                  I_Female_Pred, I_D_Allele_Pred, I_W_Allele_Pred,
                                  I_R_Allele_Pred, I_B_Allele_Pred)

## Model predictions for

Fitted_Params_Plot <- ggplot(DrosophilaModelPred, aes(x = gens, y = D_Female_Pred, color = Allele/Sex)) +
  geom_line(aes(y = D_Female_Pred, col = "Female"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = D_D_Allele_Pred, col = "D allele"), size = 1.2) +
  geom_line(aes(y = D_W_Allele_Pred, col = "W allele"), size = 1.2) +
  geom_line(aes(y = D_R_Allele_Pred, col = "R allele"), size = 1.2) +
  geom_line(aes(y = D_B_Allele_Pred, col = "B allele"), size = 1.2) +
  labs(x = "Generation", y = "Proportion") +
  ggtitle("Release 2, fitted parameters")

Ideal_Params_Plot <- ggplot(DrosophilaModelPred, aes(x = gens, y = I_Female_Pred, color = Allele/Sex)) +
  geom_line(aes(y = I_Female_Pred, col = "Female"), linetype = "dashed", size = 1.2) +
  geom_line(aes(y = I_D_Allele_Pred, col = "D allele"), size = 1.2) +
  geom_line(aes(y = I_W_Allele_Pred, col = "W allele"), size = 1.2) +
  geom_line(aes(y = I_R_Allele_Pred, col = "R allele"), size = 1.2) +
  geom_line(aes(y = I_B_Allele_Pred, col = "B allele"), size = 1.2) +
  labs(x = "Generation", y = "Proportion") +
  ggtitle("Release 2, ideal parameters")

multiplot(Male_Plot, Fitted_Params_Plot, Female_Plot, Ideal_Params_Plot, cols=2)

write.csv(homingResultsFigure,file="modelFit.csv",sep = ",")
write.csv(DrosophilaModelPred,file="modelFit2.csv",sep = ",")
