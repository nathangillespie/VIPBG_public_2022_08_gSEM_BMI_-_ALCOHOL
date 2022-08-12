install.packages("GenomicSEM")
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


require(umx)
# Covariance matrices created here:
# /Volumes/ngillespie/Documents/work/projects/2020/2020 2. genomicSEM_bmi/4. Alcohol_GWAS_results_40_PCs_separate_MFs/Alcohol_GWAS_results_40_PCs_separate_MFs.R

load("/Volumes/ngillespie/Documents/work/projects/2020/2020 2. genomicSEM_bmi/4. Alcohol_GWAS_results_40_PCs_separate_MFs/results/UKB23704_AlcCons_BMI_6ageTranche_sexStrat_20210201_f.RData")
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
 #        ALC1_f ALC2_f  ALC3_f  ALC4_f ALC5_f  ALC6_f  BMI1_f  BMI2_f BMI3_f BMI4_f BMI5_f  BMI6_f
 # ALC1_f 0.0648 0.0697  0.0530  0.0541 0.0554  0.0606  0.0170  0.0248 0.0304 0.0187 0.0213  0.0118
 # ALC2_f 0.0697 0.0906  0.0591  0.0531 0.0529  0.0527  0.0265  0.0352 0.0427 0.0393 0.0195  0.0151
 # ALC3_f 0.0530 0.0591  0.0666  0.0486 0.0524  0.0452  0.0185  0.0046 0.0192 0.0209 0.0108 -0.0001
 # ALC4_f 0.0541 0.0531  0.0486  0.0786 0.0669  0.0495  0.0045 -0.0070 0.0075 0.0006 0.0038 -0.0072
 # ALC5_f 0.0554 0.0529  0.0524  0.0669 0.0707  0.0639  0.0046  0.0082 0.0163 0.0217 0.0155  0.0107
 # ALC6_f 0.0606 0.0527  0.0452  0.0495 0.0639  0.0805 -0.0048  0.0213 0.0189 0.0144 0.0203  0.0165
 # BMI1_f 0.0170 0.0265  0.0185  0.0045 0.0046 -0.0048  0.2788  0.2680 0.2688 0.2617 0.2456  0.2281
 # BMI2_f 0.0248 0.0352  0.0046 -0.0070 0.0082  0.0213  0.2680  0.2999 0.2850 0.2822 0.2532  0.2488
 # BMI3_f 0.0304 0.0427  0.0192  0.0075 0.0163  0.0189  0.2688  0.2850 0.2792 0.2742 0.2417  0.2323
 # BMI4_f 0.0187 0.0393  0.0209  0.0006 0.0217  0.0144  0.2617  0.2822 0.2742 0.3062 0.2418  0.2472
 # BMI5_f 0.0213 0.0195  0.0108  0.0038 0.0155  0.0203  0.2456  0.2532 0.2417 0.2418 0.2508  0.2378
 # BMI6_f 0.0118 0.0151 -0.0001 -0.0072 0.0107  0.0165  0.2281  0.2488 0.2323 0.2472 0.2378  0.2361
 h2 <- t(diag2vec(round(S,4)))
 rownames(h2) <- "h2"
 h2
 #   ALC1_f ALC2_f ALC3_f ALC4_f ALC5_f ALC6_f BMI1_f BMI2_f BMI3_f BMI4_f BMI5_f BMI6_f
 # h2 0.0648 0.0906 0.0666 0.0786 0.0707 0.0805 0.2788 0.2999 0.2792 0.3062 0.2508 0.2361
 
alc_h2 <- rowMeans(as.data.frame(h2)[1, c(1:6)], na.rm=TRUE) 	# 0.0753 
bmi_h2 <- rowMeans(as.data.frame(h2)[1, c(7:12)], na.rm=TRUE)	# 0.2751667


round(cov2cor(S),4)
 #        ALC1_f ALC2_f  ALC3_f  ALC4_f ALC5_f  ALC6_f  BMI1_f  BMI2_f BMI3_f BMI4_f BMI5_f  BMI6_f
 # ALC1_f 1.0000 0.9102  0.8068  0.7582 0.8176  0.8393  0.1262  0.1776 0.2263 0.1327 0.1674  0.0958
 # ALC2_f 0.9102 1.0000  0.7609  0.6289 0.6609  0.6174  0.1668  0.2138 0.2686 0.2357 0.1295  0.1034
 # ALC3_f 0.8068 0.7609  1.0000  0.6711 0.7636  0.6176  0.1356  0.0329 0.1409 0.1465 0.0834 -0.0006
 # ALC4_f 0.7582 0.6289  0.6711  1.0000 0.8971  0.6231  0.0302 -0.0455 0.0509 0.0040 0.0273 -0.0528
 # ALC5_f 0.8176 0.6609  0.7636  0.8971 1.0000  0.8463  0.0329  0.0561 0.1161 0.1477 0.1165  0.0827
 # ALC6_f 0.8393 0.6174  0.6176  0.6231 0.8463  1.0000 -0.0318  0.1374 0.1264 0.0917 0.1427  0.1199
 # BMI1_f 0.1262 0.1668  0.1356  0.0302 0.0329 -0.0318  1.0000  0.9268 0.9635 0.8957 0.9285  0.8891
 # BMI2_f 0.1776 0.2138  0.0329 -0.0455 0.0561  0.1374  0.9268  1.0000 0.9850 0.9312 0.9230  0.9352
 # BMI3_f 0.2263 0.2686  0.1409  0.0509 0.1161  0.1264  0.9635  0.9850 1.0000 0.9380 0.9134  0.9050
 # BMI4_f 0.1327 0.2357  0.1465  0.0040 0.1477  0.0917  0.8957  0.9312 0.9380 1.0000 0.8724  0.9197
 # BMI5_f 0.1674 0.1295  0.0834  0.0273 0.1165  0.1427  0.9285  0.9230 0.9134 0.8724 1.0000  0.9773
 # BMI6_f 0.0958 0.1034 -0.0006 -0.0528 0.0827  0.1199  0.8891  0.9352 0.9050 0.9197 0.9773  1.0000

 mean(diag2vec(round(cov2cor(S),4)[7:12,1:6])) # = 0.1202167
 

# 2-factor CFA - uncorrelated
m1 <- '
F1 =~ 1*ALC1_f + ALC2_f + ALC3_f + ALC4_f + ALC5_f + ALC6_f
F2 =~ 1*BMI1_f + BMI2_f + BMI3_f + BMI4_f + BMI5_f + BMI6_f
F1 ~~ start(0.31)*F1
F2 ~~ start(0.31)*F2
BMI1_f ~~ BMI1_f
BMI2_f ~~ BMI2_f
BMI3_f ~~ BMI3_f
BMI4_f ~~ BMI4_f
BMI5_f ~~ BMI5_f
BMI6_f ~~ BMI6_f
ALC1_f ~~ ALC1_f
ALC2_f ~~ ALC2_f
ALC3_f ~~ ALC3_f
ALC4_f ~~ ALC4_f
ALC5_f ~~ ALC5_f
ALC6_f ~~ ALC6_f
'  
model_f1 <- lavaan(
	model=m1,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(model_f1)
model_f1@optim$fx


# 2-factor correlated EFA
m2 <- '
F1 =~ 1*ALC1_f + ALC2_f + ALC3_f + ALC4_f + ALC5_f + ALC6_f
F2 =~ 1*BMI1_f + BMI2_f + BMI3_f + BMI4_f + BMI5_f + BMI6_f
F1 ~~ start(0.31)*F1
F2 ~~ start(0.31)*F2
 F1 ~~ start(0.11)*F2 + corF1F2*F2
 BMI1_f ~~ BMI1_f
 BMI2_f ~~ BMI2_f
 BMI3_f ~~ BMI3_f
 BMI4_f ~~ BMI4_f
 BMI5_f ~~ BMI5_f
 BMI6_f ~~ BMI6_f
 ALC1_f ~~ ALC1_f
 ALC2_f ~~ ALC2_f
 ALC3_f ~~ ALC3_f
 ALC4_f ~~ ALC4_f
 ALC5_f ~~ ALC5_f
 ALC6_f ~~ ALC6_f
 '
model_f2 <- lavaan(
    model=m2,
    sample.cov=S,
    sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(model_f2)
model_f2@optim$fx
# 14.62668


standardizedSolution(model_f2, type = "std.all", se = TRUE, zstat = TRUE, 
                     pvalue = TRUE, ci = TRUE, level = 0.95, cov.std = TRUE)


# Model 3
# 2-factor correlated EFA, cross-sectional causality from ALC to BMI
m3 <- '
F1 =~ 1*ALC1_f + ALC2_f + ALC3_f + ALC4_f + ALC5_f + ALC6_f
F2 =~ 1*BMI1_f + BMI2_f + BMI3_f + BMI4_f + BMI5_f + BMI6_f
F1 ~~ start(0.31)*F1
F2 ~~ start(0.31)*F2
 F1 ~~ start(0.11)*F2 + corF1F2*F2
 BMI1_f ~~ BMI1_f
 BMI2_f ~~ BMI2_f
 BMI3_f ~~ BMI3_f
 BMI4_f ~~ BMI4_f
 BMI5_f ~~ BMI5_f
 BMI6_f ~~ BMI6_f
 ALC1_f ~~ ALC1_f
 ALC2_f ~~ ALC2_f
 ALC3_f ~~ ALC3_f
 ALC4_f ~~ ALC4_f
 ALC5_f ~~ ALC5_f
 ALC6_f ~~ ALC6_f
 BMI1_f ~ b1*ALC1_f + ALC1_f
 BMI2_f ~ b2*ALC2_f + ALC2_f
 BMI3_f ~ b3*ALC3_f + ALC3_f
 BMI4_f ~ b4*ALC4_f + ALC4_f
 BMI5_f ~ b5*ALC5_f + ALC5_f
 BMI6_f ~ b6*ALC6_f + ALC6_f
 '  
model_f3 <- lavaan(
	model=m3,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(model_f3)
model_f3@optim$fx




# Model 4
# 2-factor correlated EFA, cross-sectional & cross-lagged causality from ALC to BMI
m4 <- '
F1 =~ 1*ALC1_f + ALC2_f + ALC3_f + ALC4_f + ALC5_f + ALC6_f
F2 =~ 1*BMI1_f + BMI2_f + BMI3_f + BMI4_f + BMI5_f + BMI6_f
F1 ~~ start(0.31)*F1
F2 ~~ start(0.31)*F2
 F1 ~~ start(0.11)*F2 + corF1F2*F2
 BMI1_f ~~ BMI1_f
 BMI2_f ~~ BMI2_f
 BMI3_f ~~ BMI3_f
 BMI4_f ~~ BMI4_f
 BMI5_f ~~ BMI5_f
 BMI6_f ~~ BMI6_f
 ALC1_f ~~ ALC1_f
 ALC2_f ~~ ALC2_f
 ALC3_f ~~ ALC3_f
 ALC4_f ~~ ALC4_f
 ALC5_f ~~ ALC5_f
 ALC6_f ~~ ALC6_f
 BMI1_f ~ b1*ALC1_f + ALC1_f
 BMI2_f ~ b2*ALC2_f + ALC2_f
 BMI3_f ~ b3*ALC3_f + ALC3_f
 BMI4_f ~ b4*ALC4_f + ALC4_f
 BMI5_f ~ b5*ALC5_f + ALC5_f
 BMI6_f ~ b6*ALC6_f + ALC6_f
# Cross-lagged pathways
 BMI2_f ~ b7*ALC1_f + ALC1_f
 BMI3_f ~ b8*ALC2_f + ALC2_f
 BMI4_f ~ b9*ALC3_f + ALC3_f
 BMI5_f ~ b10*ALC4_f + ALC4_f
 BMI6_f ~ b11*ALC5_f + ALC5_f
 '  
model_f4 <- lavaan(
	model=m4,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(model_f4)
model_f4@optim$fx




# Model 5
# 2-factor EFA, correlated & cross-sectional causality from BMI to ALC
m5 <- '
F1 =~ 1*ALC1_f + ALC2_f + ALC3_f + ALC4_f + ALC5_f + ALC6_f
F2 =~ 1*BMI1_f + BMI2_f + BMI3_f + BMI4_f + BMI5_f + BMI6_f
F1 ~~ start(0.31)*F1
F2 ~~ start(0.31)*F2
 F1 ~~ start(0.11)*F2 + corF1F2*F2
 BMI1_f ~~ BMI1_f
 BMI2_f ~~ BMI2_f
 BMI3_f ~~ BMI3_f
 BMI4_f ~~ BMI4_f
 BMI5_f ~~ BMI5_f
 BMI6_f ~~ BMI6_f
 ALC1_f ~~ ALC1_f
 ALC2_f ~~ ALC2_f
 ALC3_f ~~ ALC3_f
 ALC4_f ~~ ALC4_f
 ALC5_f ~~ ALC5_f
 ALC6_f ~~ ALC6_f
 ALC1_f ~ b1*BMI1_f + BMI1_f
 ALC2_f ~ b2*BMI2_f + BMI2_f
 ALC3_f ~ b3*BMI3_f + BMI3_f
 ALC4_f ~ b4*BMI4_f + BMI4_f
 ALC5_f ~ b5*BMI5_f + BMI5_f
 ALC6_f ~ b6*BMI6_f + BMI6_f
 '  
model_f5 <- lavaan(
	model=m5,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(model_f5)
model_f5@optim$fx

standardizedSolution(model_f5, type = "std.all", se = TRUE, zstat = TRUE, 
                     pvalue = TRUE, ci = TRUE, level = 0.95, cov.std = TRUE)



			         lhs op    rhs   label est.std    se      z pvalue ci.lower ci.upper
			   1      F1 =~ ALC1_f           0.955 0.261  3.661  0.000    0.444    1.466
			   2      F1 =~ ALC2_f           0.769 0.119  6.435  0.000    0.534    1.003
			   3      F1 =~ ALC3_f           0.808 0.141  5.727  0.000    0.531    1.084
			   4      F1 =~ ALC4_f           0.858 0.124  6.928  0.000    0.615    1.100
			   5      F1 =~ ALC5_f           0.971 0.112  8.709  0.000    0.753    1.190
			   6      F1 =~ ALC6_f           0.802 0.137  5.868  0.000    0.534    1.070
			   
			   7      F2 =~ BMI1_f           0.959 0.046 20.822  0.000    0.869    1.050
			   8      F2 =~ BMI2_f           0.979 0.035 27.822  0.000    0.910    1.048
			   9      F2 =~ BMI3_f           0.980 0.030 33.078  0.000    0.922    1.039
			   10     F2 =~ BMI4_f           0.941 0.025 37.550  0.000    0.892    0.990
			   11     F2 =~ BMI5_f           0.953 0.024 39.679  0.000    0.906    1.000
			   12     F2 =~ BMI6_f           0.967 0.030 32.118  0.000    0.908    1.026
			   
			   13     F1 ~~     F1           1.000 0.000     NA     NA    1.000    1.000
			   14     F2 ~~     F2           1.000 0.000     NA     NA    1.000    1.000
			   15     F1 ~~     F2 corF1F2  -0.069 0.423 -0.163  0.870   -0.897    0.759
			   
			   16 BMI1_f ~~ BMI1_f           0.080 0.088  0.901  0.368   -0.094    0.253
			   17 BMI2_f ~~ BMI2_f           0.041 0.069  0.589  0.556   -0.095    0.176
			   18 BMI3_f ~~ BMI3_f           0.039 0.058  0.665  0.506   -0.075    0.153
			   19 BMI4_f ~~ BMI4_f           0.114 0.047  2.412  0.016    0.021    0.206
			   20 BMI5_f ~~ BMI5_f           0.092 0.046  2.015  0.044    0.003    0.182
			   21 BMI6_f ~~ BMI6_f           0.066 0.058  1.127  0.260   -0.048    0.180
			   22 ALC1_f ~~ ALC1_f           0.060 0.511  0.117  0.907   -0.942    1.062
			   23 ALC2_f ~~ ALC2_f           0.371 0.186  2.000  0.045    0.007    0.735
			   24 ALC3_f ~~ ALC3_f           0.341 0.228  1.494  0.135   -0.106    0.787
			   25 ALC4_f ~~ ALC4_f           0.268 0.208  1.286  0.199   -0.140    0.675
			   26 ALC5_f ~~ ALC5_f           0.048 0.213  0.228  0.820   -0.369    0.465
			   27 ALC6_f ~~ ALC6_f           0.345 0.218  1.586  0.113   -0.081    0.772
			   
			   28 ALC1_f  ~ BMI1_f      b1   0.242 0.428  0.564  0.572   -0.598    1.082
			   29 ALC2_f  ~ BMI2_f      b2   0.254 0.341  0.744  0.457   -0.415    0.922
			   30 ALC3_f  ~ BMI3_f      b3   0.155 0.348  0.446  0.656   -0.526    0.836
			   31 ALC4_f  ~ BMI4_f      b4   0.065 0.374  0.174  0.862   -0.669    0.799
			   32 ALC5_f  ~ BMI5_f      b5   0.175 0.426  0.410  0.682   -0.660    1.009
			   33 ALC6_f  ~ BMI6_f      b6   0.173 0.356  0.487  0.626   -0.525    0.872


# Model 6
# 2-factor EFA, correlated & cross-sectional & cross-lagged causality from BMI to ALC
m6 <- '
F1 =~ 1*ALC1_f + ALC2_f + ALC3_f + ALC4_f + ALC5_f + ALC6_f
F2 =~ 1*BMI1_f + BMI2_f + BMI3_f + BMI4_f + BMI5_f + BMI6_f
F1 ~~ start(0.31)*F1
F2 ~~ start(0.31)*F2
 F1 ~~ start(0.11)*F2 + corF1F2*F2
 BMI1_f ~~ BMI1_f
 BMI2_f ~~ BMI2_f
 BMI3_f ~~ BMI3_f
 BMI4_f ~~ BMI4_f
 BMI5_f ~~ BMI5_f
 BMI6_f ~~ BMI6_f
 ALC1_f ~~ ALC1_f
 ALC2_f ~~ ALC2_f
 ALC3_f ~~ ALC3_f
 ALC4_f ~~ ALC4_f
 ALC5_f ~~ ALC5_f
 ALC6_f ~~ ALC6_f
 ALC1_f ~ b1*BMI1_f + BMI1_f
 ALC2_f ~ b2*BMI2_f + BMI2_f
 ALC3_f ~ b3*BMI3_f + BMI3_f
 ALC4_f ~ b4*BMI4_f + BMI4_f
 ALC5_f ~ b5*BMI5_f + BMI5_f
 ALC6_f ~ b6*BMI6_f + BMI6_f
 ALC2_f ~ b7*BMI1_f + BMI1_f
 ALC3_f ~ b8*BMI2_f + BMI2_f
 ALC4_f ~ b9*BMI3_f + BMI3_f
 ALC5_f ~ b10*BMI4_f + BMI4_f
 ALC6_f ~ b11*BMI5_f + BMI5_f
 '  
model_f6 <- lavaan(
	model=m6,
	sample.cov=S,
	sample.nobs=2,
	NACOV=V,
	WLS.V=W,
	optim.method="BFGS",
	estimator="WLSMV"
)
summary(model_f6)
model_f6@optim$fx

parameterEstimates(model_f6)
standardizedSolution(model_f6)
fitMeasures(model_f5)
fitMeasures(model_f6, "chisq")
lavInspect(model_f6)
lavInspect(model_f6, what = "list")
#AIC(model_f6)
#BIC(model_f6)

#Pseudo-AICs:
model_f1_aic 	<- model_f1@Fit@test$scaled.shifted$stat 	  + 2*model_f1@Fit@npar
model_f2_aic 	<- model_f2@Fit@test$scaled.shifted$stat 	  + 2*model_f2@Fit@npar
model_f3_aic 	<- model_f3@Fit@test$scaled.shifted$stat 	  + 2*model_f3@Fit@npar
model_f4_aic 	<- model_f4@Fit@test$scaled.shifted$stat 	  + 2*model_f4@Fit@npar
model_f5_aic 	<- model_f5@Fit@test$scaled.shifted$stat 	  + 2*model_f5@Fit@npar
model_f6_aic 	<- model_f6@Fit@test$scaled.shifted$stat 	  + 2*model_f6@Fit@npar

AICs <- rbind(
model_f1_aic, 	
model_f2_aic ,	
model_f3_aic ,	
model_f4_aic ,	
model_f5_aic ,	
model_f6_aic )

fits <- rbind(
fitmeasures(model_f1)[c("npar","chisq","srmr","cfi","tli")],
fitmeasures(model_f2)[c("npar","chisq","srmr","cfi","tli")],
fitmeasures(model_f3)[c("npar","chisq","srmr","cfi","tli")],
fitmeasures(model_f4)[c("npar","chisq","srmr","cfi","tli")],
fitmeasures(model_f5)[c("npar","chisq","srmr","cfi","tli")],
fitmeasures(model_f6)[c("npar","chisq","srmr","cfi","tli")])

results <- cbind(fits[,1],AICs,fits[,2:5])
colnames(results) <- c("npar", "Chisq","Pseudo-AIC","srmr","cfi","tli")
rownames(results) <- c("M1","M2","M3","M4","M5","M6")
results
          npar     Chisq Pseudo-AIC       srmr       cfi       tli
Model_f1    24 119.71707  147.23354 0.09450189 0.9931160 0.9915862
Model_f2    25  86.92903   58.50671 0.06188796 0.9995934 0.9994937
Model_f3    31  94.84963   51.86373 0.06853065 0.9996409 0.9994957
Model_f3b   36  95.71399   35.26884 0.06544781 1.0000000 1.0007810
Model_4     31  85.64536   31.30600 0.04876278 1.0000000 1.0016272
Model_4b    36  93.20841   28.40128 0.04614188 1.0000000 1.0015778

write.csv(noquote(results) ,file="/Users/ngillespie/Desktop/temp.csv")



