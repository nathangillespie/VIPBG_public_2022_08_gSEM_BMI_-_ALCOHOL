library(GenomicSEM)
library(lavaan)

# Read in S and V matrices
# Autoregression models
# Factor analyses 
# 1 factor
# 2 factors
# 3 factors
# Model comparisons



# Read in S and V matrices


require(OpenMx)
# Covariance matrices created here:
# /Volumes/ngillespie/Documents/work/projects/2020/2020 2. genomicSEM_bmi/4. Alcohol_GWAS_results_40_PCs_separate_MFs/Alcohol_GWAS_results_40_PCs_separate_MFs.R

load("/Volumes/ngillespie/Documents/work/projects/2020/2020 2. genomicSEM_bmi/4. Alcohol_GWAS_results_40_PCs_separate_MFs/results/UKB23704_AlcCons_BMI_6ageTranche_sexStrat_20210201_m.RData")
S <- NA
V <- NA
W <- NA
S <- LDSCoutput$S
S <- as.matrix(Matrix::nearPD(S)$mat)
rownames(S) <- colnames(S)
cov2cor(S)
V <- LDSCoutput$V
W 					<- chol2inv(chol(LDSCoutput$V))
#dataDWLS 			<- mxData(S, numObs = 2, means = NA, type = "acov", acov=diag(diag(V)), fullWeight=W)

round(S,4)
 #        ALC1_m ALC2_m ALC3_m ALC4_m  ALC5_m  ALC6_m  BMI1_m  BMI2_m  BMI3_m BMI4_m BMI5_m BMI6_m
 # ALC1_m 0.0997 0.0888 0.0798 0.0769  0.0729  0.0767  0.0208  0.0134  0.0305 0.0335 0.0332 0.0492
 # ALC2_m 0.0888 0.1130 0.0865 0.0779  0.0658  0.0893  0.0102  0.0141  0.0486 0.0403 0.0211 0.0265
 # ALC3_m 0.0798 0.0865 0.0926 0.0768  0.0840  0.0789  0.0027  0.0085  0.0143 0.0276 0.0208 0.0138
 # ALC4_m 0.0769 0.0779 0.0768 0.0881  0.0671  0.0763  0.0174  0.0321  0.0366 0.0416 0.0322 0.0382
 # ALC5_m 0.0729 0.0658 0.0840 0.0671  0.0852  0.0673 -0.0026 -0.0055 -0.0103 0.0080 0.0048 0.0051
 # ALC6_m 0.0767 0.0893 0.0789 0.0763  0.0673  0.0819 -0.0107  0.0046  0.0150 0.0194 0.0004 0.0103
 # BMI1_m 0.0208 0.0102 0.0027 0.0174 -0.0026 -0.0107  0.3054  0.2143  0.2609 0.2448 0.2105 0.2833
 # BMI2_m 0.0134 0.0141 0.0085 0.0321 -0.0055  0.0046  0.2143  0.2015  0.2079 0.2208 0.1893 0.2338
 # BMI3_m 0.0305 0.0486 0.0143 0.0366 -0.0103  0.0150  0.2609  0.2079  0.2668 0.2413 0.2031 0.2547
 # BMI4_m 0.0335 0.0403 0.0276 0.0416  0.0080  0.0194  0.2448  0.2208  0.2413 0.2530 0.2165 0.2647
 # BMI5_m 0.0332 0.0211 0.0208 0.0322  0.0048  0.0004  0.2105  0.1893  0.2031 0.2165 0.2220 0.2309
 # BMI6_m 0.0492 0.0265 0.0138 0.0382  0.0051  0.0103  0.2833  0.2338  0.2547 0.2647 0.2309 0.3154

 # SNP heritabilities
 h2 <- t(diag2vec(round(S,4)))
 rownames(h2) <- "h2"
 h2
 #    ALC1_m ALC2_m ALC3_m ALC4_m ALC5_m ALC6_m BMI1_m BMI2_m BMI3_m BMI4_m BMI5_m BMI6_m
 # h2 0.0997  0.113 0.0926 0.0881 0.0852 0.0819 0.3054 0.2015 0.2668  0.253  0.222 0.3154

 # Average genome-wdie SNP heritabilities 
 alc_h2 <- rowMeans(as.data.frame(h2)[1, c(1:6)], na.rm=TRUE) 	# 0.09
 bmi_h2 <- rowMeans(as.data.frame(h2)[1, c(7:12)], na.rm=TRUE)	# 0.26

 # Genome-wide SNP covariances 
 cov <- t(diag2vec(round(S[7:12,1:6],4)))
 rownames(cov) <- "covs"
 colnames(cov) <- c("cov1","cov2","cov3","cov4","cov5","cov6"); 
 #        cov1   cov2   cov3   cov4   cov5   cov6
 # covs 0.0208 0.0141 0.0143 0.0416 0.0048 0.0103
 
 # Average genome-wide SNP covariance
 rowMeans(as.data.frame(cov)[1, c(1:6)], na.rm=TRUE) 	
 # 0.01765 



round(cov2cor(S),4)
 #        ALC1_m ALC2_m ALC3_m ALC4_m  ALC5_m  ALC6_m  BMI1_m  BMI2_m  BMI3_m BMI4_m BMI5_m BMI6_m
 # ALC1_m 1.0000 0.8362 0.8301 0.8198  0.7912  0.8488  0.1193  0.0943  0.1871 0.2108 0.2228 0.2774
 # ALC2_m 0.8362 1.0000 0.8453 0.7807  0.6702  0.9280  0.0550  0.0932  0.2801 0.2381 0.1334 0.1406
 # ALC3_m 0.8301 0.8453 1.0000 0.8504  0.9451  0.9060  0.0159  0.0623  0.0909 0.1805 0.1447 0.0804
 # ALC4_m 0.8198 0.7807 0.8504 1.0000  0.7740  0.8983  0.1058  0.2408  0.2385 0.2785 0.2305 0.2292
 # ALC5_m 0.7912 0.6702 0.9451 0.7740  1.0000  0.8058 -0.0160 -0.0418 -0.0681 0.0544 0.0352 0.0308
 # ALC6_m 0.8488 0.9280 0.9060 0.8983  0.8058  1.0000 -0.0678  0.0361  0.1012 0.1347 0.0033 0.0640
 # BMI1_m 0.1193 0.0550 0.0159 0.1058 -0.0160 -0.0678  1.0000  0.8640  0.9141 0.8805 0.8084 0.9128
 # BMI2_m 0.0943 0.0932 0.0623 0.2408 -0.0418  0.0361  0.8640  1.0000  0.8966 0.9779 0.8950 0.9276
 # BMI3_m 0.1871 0.2801 0.0909 0.2385 -0.0681  0.1012  0.9141  0.8966  1.0000 0.9288 0.8346 0.8781
 # BMI4_m 0.2108 0.2381 0.1805 0.2785  0.0544  0.1347  0.8805  0.9779  0.9288 1.0000 0.9134 0.9370
 # BMI5_m 0.2228 0.1334 0.1447 0.2305  0.0352  0.0033  0.8084  0.8950  0.8346 0.9134 1.0000 0.8726
 # BMI6_m 0.2774 0.1406 0.0804 0.2292  0.0308  0.0640  0.9128  0.9276  0.8781 0.9370 0.8726 1.0000
 mean(diag2vec(round(cov2cor(S),4)[7:12,1:6])) #= 0.1135167


# 2-factor CFA - uncorrelated
m1 <- '
F1 =~ 1*ALC1_m + ALC2_m + ALC3_m + ALC4_m + ALC5_m + ALC6_m
F2 =~ 1*BMI1_m + BMI2_m + BMI3_m + BMI4_m + BMI5_m + BMI6_m
F1 ~~ start(0.11)*F1
F2 ~~ start(0.12)*F2
BMI1_m ~~ BMI1_m
BMI2_m ~~ BMI2_m
BMI3_m ~~ BMI3_m
BMI4_m ~~ BMI4_m
BMI5_m ~~ BMI5_m
BMI6_m ~~ BMI6_m
ALC1_m ~~ ALC1_m
ALC2_m ~~ ALC2_m
ALC3_m ~~ ALC3_m
ALC4_m ~~ ALC4_m
ALC5_m ~~ ALC5_m
ALC6_m ~~ ALC6_m
'  
model_m1 <- lavaan(
	model=m1,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(model_m1)
model_m1@optim$fx


# 2-factor correlated EFA
m2 <- '
 F1 =~ 1*ALC1_m + ALC2_m + ALC3_m + ALC4_m + ALC5_m + ALC6_m
 F2 =~ 1*BMI1_m + BMI2_m + BMI3_m + BMI4_m + BMI5_m + BMI6_m
 F1 ~~ start(0.11)*F1
 F2 ~~ start(0.12)*F2
 F1 ~~ start(0.11)*F2 + corF1F2*F2
 BMI1_m ~~ BMI1_m
 BMI2_m ~~ BMI2_m
 BMI3_m ~~ BMI3_m
 BMI4_m ~~ BMI4_m
 BMI5_m ~~ BMI5_m
 BMI6_m ~~ BMI6_m
 ALC1_m ~~ ALC1_m
 ALC2_m ~~ ALC2_m
 ALC3_m ~~ ALC3_m
 ALC4_m ~~ ALC4_m
 ALC5_m ~~ ALC5_m
 ALC6_m ~~ ALC6_m
 '
model_m2 <- lavaan(
    model=m2,
    sample.cov=S,
    sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(model_m2)
model_m2@optim$fx

standardizedSolution(model_m2, type = "std.all", se = TRUE, zstat = TRUE, 
                     pvalue = TRUE, ci = TRUE, level = 0.95, cov.std = TRUE)



# Model 3
# 2-factor correlated EFA, cross-sectional causality from ALC to BMI
m3 <- '
F1 =~ 1*ALC1_m + ALC2_m + ALC3_m + ALC4_m + ALC5_m + ALC6_m
F2 =~ 1*BMI1_m + BMI2_m + BMI3_m + BMI4_m + BMI5_m + BMI6_m
F1 ~~ start(0.11)*F1
F2 ~~ start(0.12)*F2
F1 ~~ start(0.11)*F2 + corF1F2*F2
 BMI1_m ~~ BMI1_m
 BMI2_m ~~ BMI2_m
 BMI3_m ~~ BMI3_m
 BMI4_m ~~ BMI4_m
 BMI5_m ~~ BMI5_m
 BMI6_m ~~ BMI6_m
 ALC1_m ~~ ALC1_m
 ALC2_m ~~ ALC2_m
 ALC3_m ~~ ALC3_m
 ALC4_m ~~ ALC4_m
 ALC5_m ~~ ALC5_m
 ALC6_m ~~ ALC6_m
 BMI1_m ~ b1*ALC1_m + ALC1_m
 BMI2_m ~ b2*ALC2_m + ALC2_m
 BMI3_m ~ b3*ALC3_m + ALC3_m
 BMI4_m ~ b4*ALC4_m + ALC4_m
 BMI5_m ~ b5*ALC5_m + ALC5_m
 BMI6_m ~ b6*ALC6_m + ALC6_m
 '  
model_m3 <- lavaan(
	model=m3,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(model_m3)
model_m3@optim$fx




# Model 4
# 2-factor correlated EFA, cross-sectional & cross-lagged causality from ALC to BMI
m3 <- '
F1 =~ 1*ALC1_m + ALC2_m + ALC3_m + ALC4_m + ALC5_m + ALC6_m
F2 =~ 1*BMI1_m + BMI2_m + BMI3_m + BMI4_m + BMI5_m + BMI6_m
F1 ~~ start(0.11)*F1
F2 ~~ start(0.12)*F2
F1 ~~ start(0.11)*F2 + corF1F2*F2
 BMI1_m ~~ BMI1_m
 BMI2_m ~~ BMI2_m
 BMI3_m ~~ BMI3_m
 BMI4_m ~~ BMI4_m
 BMI5_m ~~ BMI5_m
 BMI6_m ~~ BMI6_m
 ALC1_m ~~ ALC1_m
 ALC2_m ~~ ALC2_m
 ALC3_m ~~ ALC3_m
 ALC4_m ~~ ALC4_m
 ALC5_m ~~ ALC5_m
 ALC6_m ~~ ALC6_m
 BMI1_m ~ b1*ALC1_m + ALC1_m
 BMI2_m ~ b2*ALC2_m + ALC2_m
 BMI3_m ~ b3*ALC3_m + ALC3_m
 BMI4_m ~ b4*ALC4_m + ALC4_m
 BMI5_m ~ b5*ALC5_m + ALC5_m
 BMI6_m ~ b6*ALC6_m + ALC6_m
# Cross-lagged pathways
 BMI2_m ~ b7*ALC1_m + ALC1_m
 BMI3_m ~ b8*ALC2_m + ALC2_m
 BMI4_m ~ b9*ALC3_m + ALC3_m
 BMI5_m ~ b10*ALC4_m + ALC4_m
 BMI6_m ~ b11*ALC5_m + ALC5_m
 '  
model_m3 <- lavaan(
	model=m3,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(model_m3)
model_m3@optim$fx




# Model 5
# 2-factor EFA, correlated & cross-sectional causality from BMI to ALC
m5 <- '
F1 =~ 1*ALC1_m + ALC2_m + ALC3_m + ALC4_m + ALC5_m + ALC6_m
F2 =~ 1*BMI1_m + BMI2_m + BMI3_m + BMI4_m + BMI5_m + BMI6_m
F1 ~~ start(0.11)*F1
F2 ~~ start(0.12)*F2
F1 ~~ start(0.11)*F2 + corF1F2*F2
 BMI1_m ~~ BMI1_m
 BMI2_m ~~ BMI2_m
 BMI3_m ~~ BMI3_m
 BMI4_m ~~ BMI4_m
 BMI5_m ~~ BMI5_m
 BMI6_m ~~ BMI6_m
 ALC1_m ~~ ALC1_m
 ALC2_m ~~ ALC2_m
 ALC3_m ~~ ALC3_m
 ALC4_m ~~ ALC4_m
 ALC5_m ~~ ALC5_m
 ALC6_m ~~ ALC6_m
 ALC1_m ~ b1*BMI1_m + BMI1_m
 ALC2_m ~ b2*BMI2_m + BMI2_m
 ALC3_m ~ b3*BMI3_m + BMI3_m
 ALC4_m ~ b4*BMI4_m + BMI4_m
 ALC5_m ~ b5*BMI5_m + BMI5_m
 ALC6_m ~ b6*BMI6_m + BMI6_m
 '  
model_m5 <- lavaan(
	model=m5,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(model_m5)
model_m5@optim$fx


standardizedSolution(model_m5, type = "std.all", se = TRUE, zstat = TRUE, 
                     pvalue = TRUE, ci = TRUE, level = 0.95, cov.std = TRUE)



			         lhs op    rhs   label est.std    se      z pvalue ci.lower ci.upper
			   1      F1 =~ ALC1_m           0.880 0.182  4.824  0.000    0.523    1.238
			   2      F1 =~ ALC2_m           0.854 0.125  6.840  0.000    0.610    1.099
			   3      F1 =~ ALC3_m           0.979 0.134  7.320  0.000    0.717    1.241
			   4      F1 =~ ALC4_m           0.881 0.126  7.010  0.000    0.634    1.127
			   5      F1 =~ ALC5_m           0.867 0.105  8.233  0.000    0.661    1.073
			   6      F1 =~ ALC6_m           0.973 0.114  8.519  0.000    0.749    1.197
			   7      F2 =~ BMI1_m           0.922 0.051 18.115  0.000    0.823    1.022
			   8      F2 =~ BMI2_m           0.979 0.055 17.684  0.000    0.870    1.087
			   9      F2 =~ BMI3_m           0.935 0.039 23.868  0.000    0.858    1.012
			   10     F2 =~ BMI4_m           0.988 0.035 28.477  0.000    0.920    1.056
			   11     F2 =~ BMI5_m           0.909 0.031 29.140  0.000    0.848    0.970
			   12     F2 =~ BMI6_m           0.956 0.028 34.494  0.000    0.901    1.010
			   13     F1 ~~     F1           1.000 0.000     NA     NA    1.000    1.000
			   14     F2 ~~     F2           1.000 0.000     NA     NA    1.000    1.000
			   15     F1 ~~     F2 corF1F2   0.027 0.292  0.092  0.927   -0.545    0.598
			   
			   
			   16 BMI1_m ~~ BMI1_m           0.149 0.094  1.588  0.112   -0.035    0.333
			   17 BMI2_m ~~ BMI2_m           0.042 0.108  0.384  0.701   -0.171    0.254
			   18 BMI3_m ~~ BMI3_m           0.126 0.073  1.722  0.085   -0.017    0.270
			   19 BMI4_m ~~ BMI4_m           0.023 0.069  0.340  0.734   -0.111    0.158
			   20 BMI5_m ~~ BMI5_m           0.173 0.057  3.051  0.002    0.062    0.284
			   21 BMI6_m ~~ BMI6_m           0.086 0.053  1.633  0.103   -0.017    0.190
			   22 ALC1_m ~~ ALC1_m           0.176 0.328  0.536  0.592   -0.468    0.820
			   23 ALC2_m ~~ ALC2_m           0.240 0.219  1.099  0.272   -0.188    0.669
			   24 ALC3_m ~~ ALC3_m           0.030 0.262  0.116  0.908   -0.482    0.543
			   25 ALC4_m ~~ ALC4_m           0.166 0.230  0.720  0.471   -0.285    0.616
			   26 ALC5_m ~~ ALC5_m           0.249 0.182  1.364  0.173   -0.109    0.606
			   27 ALC6_m ~~ ALC6_m           0.050 0.223  0.224  0.823   -0.388    0.488
			   
			   28 ALC1_m  ~ BMI1_m      b1   0.201 0.289  0.696  0.486   -0.365    0.768
			   29 ALC2_m  ~ BMI2_m      b2   0.152 0.262  0.579  0.563   -0.362    0.665
			   30 ALC3_m  ~ BMI3_m      b3   0.087 0.300  0.291  0.771   -0.500    0.675
			   31 ALC4_m  ~ BMI4_m      b4   0.221 0.265  0.831  0.406   -0.299    0.741
			   32 ALC5_m  ~ BMI5_m      b5  -0.016 0.276 -0.059  0.953   -0.557    0.525
			   33 ALC6_m  ~ BMI6_m      b6   0.032 0.291  0.112  0.911   -0.537    0.602


# Model 6
# 2-factor EFA, correlated & cross-sectional & cross-lagged causality from BMI to ALC
m6 <- '
F1 =~ 1*ALC1_m + ALC2_m + ALC3_m + ALC4_m + ALC5_m + ALC6_m
F2 =~ 1*BMI1_m + BMI2_m + BMI3_m + BMI4_m + BMI5_m + BMI6_m
F1 ~~ start(0.31)*F1
F2 ~~ start(0.31)*F2
 F1 ~~ start(0.11)*F2 + corF1F2*F2
 BMI1_m ~~ BMI1_m
 BMI2_m ~~ BMI2_m
 BMI3_m ~~ BMI3_m
 BMI4_m ~~ BMI4_m
 BMI5_m ~~ BMI5_m
 BMI6_m ~~ BMI6_m
 ALC1_m ~~ ALC1_m
 ALC2_m ~~ ALC2_m
 ALC3_m ~~ ALC3_m
 ALC4_m ~~ ALC4_m
 ALC5_m ~~ ALC5_m
 ALC6_m ~~ ALC6_m
 ALC1_m ~ b1*BMI1_m + BMI1_m
 ALC2_m ~ b2*BMI2_m + BMI2_m
 ALC3_m ~ b3*BMI3_m + BMI3_m
 ALC4_m ~ b4*BMI4_m + BMI4_m
 ALC5_m ~ b5*BMI5_m + BMI5_m
 ALC6_m ~ b6*BMI6_m + BMI6_m
 ALC2_m ~ b7*BMI1_m + BMI1_m
 ALC3_m ~ b8*BMI2_m + BMI2_m
 ALC4_m ~ b9*BMI3_m + BMI3_m
 ALC5_m ~ b10*BMI4_m + BMI4_m
 ALC6_m ~ b11*BMI5_m + BMI5_m
 '  
model_m6 <- lavaan(
	model=m6,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(model_m6)
model_m6@optim$fx


standardizedSolution(model_m6, type = "std.all", se = TRUE, zstat = TRUE, 
                     pvalue = TRUE, ci = TRUE, level = 0.95, cov.std = TRUE)




parameterEstimates(model_m6)
standardizedSolution(model_m6)
fitMeasures(model_m6)
fitMeasures(model_m6, "chisq")
lavInspect(model_m6)
lavInspect(model_m6, what = "list")
#AIC(model_m6)
#BIC(model_m6)

#Pseudo-AICs:
model_m1_aic 	<- model_m1@Fit@test$scaled.shifted$stat 	  + 2*model_m1@Fit@npar
model_m2_aic 	<- model_m2@Fit@test$scaled.shifted$stat 	  + 2*model_m2@Fit@npar
model_m3_aic 	<- model_m3@Fit@test$scaled.shifted$stat 	  + 2*model_m3@Fit@npar
model_m3_aic 	<- model_m3@Fit@test$scaled.shifted$stat 	  + 2*model_m3@Fit@npar
model_m5_aic 	<- model_m5@Fit@test$scaled.shifted$stat 	  + 2*model_m5@Fit@npar
model_m6_aic 	<- model_m6@Fit@test$scaled.shifted$stat 	  + 2*model_m6@Fit@npar

AICs <- rbind(
model_m1_aic, 	
model_m2_aic ,	
model_m3_aic ,	
model_m3_aic ,	
model_m5_aic ,	
model_m6_aic )

fits <- rbind(
fitmeasures(model_m1)[c("npar","chisq","srmr","cfi","tli")],
fitmeasures(model_m2)[c("npar","chisq","srmr","cfi","tli")],
fitmeasures(model_m3)[c("npar","chisq","srmr","cfi","tli")],
fitmeasures(model_m3)[c("npar","chisq","srmr","cfi","tli")],
fitmeasures(model_m5)[c("npar","chisq","srmr","cfi","tli")],
fitmeasures(model_m6)[c("npar","chisq","srmr","cfi","tli")])

results <- cbind(fits[,1:2],AICs,fits[,3:5])
colnames(results) <- c("NP","Chisq","Pseudo-AIC","SRMR","CFI","TLI")
rownames(results) <- c("M1","M2","M3","M4","M5","M6")
results
write.csv(noquote(results) ,file="/Users/ngillespie/Desktop/temp.csv")

   NP     Chisq Pseudo-AIC       SRMR       CFI       TLI
M1 24 203.62383  136.28479 0.10673639 0.9842043 0.9806942
M2 25  88.82218  100.22561 0.07376391 0.9962183 0.9952907
M3 36  52.53878  104.63903 0.05931686 0.9988874 0.9982517
M4 36  52.53878  104.63903 0.05931686 0.9988874 0.9982517
M5 31  34.04005   87.13642 0.04801438 1.0000000 1.0019213
M6 36  28.21210   93.09329 0.04292542 1.0000000 1.0022873