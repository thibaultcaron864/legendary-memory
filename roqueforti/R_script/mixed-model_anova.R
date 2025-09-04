###################################################################
###  ch16.r                                                     ###
###                                                             ###
###  This is an R script for producing the examples in          ###
###  Chapter 16 of                                              ###
###  Snijders, Tom A.B., and Bosker, Roel J.                    ###
###  Multilevel Analysis:                                       ###
###  An Introduction to Basic and Advanced Multilevel Modeling, ###
###  second edition.                                            ### 
###  London etc.: Sage Publishers, 2012                         ### 
###                                                             ###
###  Thanks to Ruth Ripley for helpful guidance.                ###
###                                                             ###
###  version March 14, 2017                                     ### 
###                                                             ###
###################################################################

# Read data
mlbook_b <- read.table("mlbook2_b.dat", header=TRUE)
# what are the names
names(mlbook_b)

# Note that school means can be calculated by functions such as
# schmeans <- aggregate(mlbook_red, by = list(mlbook_red$schoolnr), mean)
# But this is not necessary here because the school means
# are already in the data set.

# First we have to reorganize the multivariate data 
# with an extra lowest level indicating the responses,
# as described on page 284.
# Copy the data, with extra columns:
# the response variable, and dummy variables whichpart indicating
# whether this is the language or the arithmetic response.

part1 <- cbind(mlbook_b, response = mlbook_b$langPOST, whichpart = factor(1))
part2 <- cbind(mlbook_b, response = mlbook_b$aritPOST, whichpart = factor(2))
mlbook_disaggregated <- rbind(part1, part2)
mlbook_disaggregated <- 
  mlbook_disaggregated[order(mlbook_disaggregated$schoolnr,
                             mlbook_disaggregated$pupilNR_new),]
names(mlbook_disaggregated)
rownames(mlbook_disaggregated) <- 1:dim(mlbook_disaggregated)[1]
# A short look at what has been produced:
mlbook_disaggregated[1:10,]
# Omit rows with missings in response variable:
mlbook_disaggregated <- mlbook_disaggregated[!is.na(mlbook_disaggregated$response),]
# How much is left:
dim(mlbook_disaggregated)

# Attach nlme
library(nlme)

# The next command yields Table 16.1.
# Having "- 1" as part of a formula drops the intercept.
# For weights=varIdent, see Pinheiro & Bates (2000), page 208-209.
# For corr=corSymm, see Pinheiro & Bates (2000), page 234-235.

mlb161 <- lme(response ~ - 1 + whichpart, random = ~ -1 + whichpart|schoolnr, 
              weights=varIdent(form=~1|whichpart),
              corr=corSymm(form=~as.numeric(whichpart)|schoolnr/pupilNR_new),
              data=mlbook_disaggregated, method="ML",
              control = list(maxIter=500, msMaxIter=500, tolerance=1e-8,
                             niterEM=250))
summary(mlb161)

# Some notes about the extraction of estimates of parameters in the random part:
# The covariance matrix at the school level and the residual variance for
# the last category (arithmetic) of "whichpart":
VarCorr(mlb161)
# The factors for the various levels of the varIdent structure are obtained as follows:
mlb161$modelStruct$varStruct
# These numbers give the proportionality constants between the standard deviations.
# The internal coefficients are the logarithms, and
exp(coef(mlb161$modelStruct$varStruct))
# gives the same result.
# Except for rounding and numerical precision issues,
# this is the same as the ratio of standard deviations
sqrt(62.87/32.12)
# as given according to Table 16.1.

# The next command yields (almost) Table 16.2.
mlb162 <- lme(response ~ - 1 + whichpart + whichpart:(IQ_verb*ses + sch_iqv*sch_ses), 
              random = ~ -1 + whichpart|schoolnr, 
              weights=varIdent(form=~1|whichpart),
              corr=corSymm(form=~as.numeric(whichpart)|schoolnr/pupilNR_new),
              data=mlbook_disaggregated, method="ML",
              control = list(maxIter=500, msMaxIter=500, tolerance=1e-8,
                             niterEM=250))
summary(mlb162)

# The results differ slightly from those in the book,
# mainly the parameter estimates for mean_iqv (in the 2nd decimal).
# Substantively this is not important.
# Perhaps this is due to different approaches to rounding at some stage
# for the calculation and storage of the school means.
