# Alchol GWAS results 

# Unzip gwas & concatenate
# Prepare summary statistics in R
# Munge summary statistics files

# Genetic covariance matrices
# Alcohol - women
# Alcohol - men
# Alcohol - men + women
# Alcohol & bmi - men
# Alcohol & bmi - women


Print the columns
1 The rsID of the SNP.
2 An A1 allele column, with A1 indicating the effect allele.
3 An A2 allele column, with A2 indicating the non-effect allele.
4 Either a logistic or continuous regression effect.



# Alcohol

    
zcat UKB23704_AlcCons_6ageTranche_sexStrat_20210201_chr20.gz | head -1  | tr ' ' '\n' | awk '{print NR,$1}'	

 # 1 chr
 # 2 rsid
 # 3 pos
 # 4 a_0
 # 5 a_1
 # 6 af
 # 7 info
 # 8 alc_t6_female_a_beta
 # 9 alc_t6_female_a_se
 # 10 alc_t6_female_a_t
 # 11 alc_t6_female_a-log10p
 # 12 alc_t6_female_b_beta
 # 13 alc_t6_female_b_se
 # 14 alc_t6_female_b_t
 # 15 alc_t6_female_b-log10p
 # 16 alc_t6_female_c_beta
 # 17 alc_t6_female_c_se
 # 18 alc_t6_female_c_t
 # 19 alc_t6_female_c-log10p
 # 20 alc_t6_female_d_beta
 # 21 alc_t6_female_d_se
 # 22 alc_t6_female_d_t
 # 23 alc_t6_female_d-log10p
 # 24 alc_t6_female_e_beta
 # 25 alc_t6_female_e_se
 # 26 alc_t6_female_e_t
 # 27 alc_t6_female_e-log10p
 # 28 alc_t6_female_f_beta
 # 29 alc_t6_female_f_se
 # 30 alc_t6_female_f_t
 # 31 alc_t6_female_f-log10p
 #  32 alc_t6_male_a_beta
 # 33 alc_t6_male_a_se
 # 34 alc_t6_male_a_t
 # 35 alc_t6_male_a-log10p
 # 36 alc_t6_male_b_beta
 # 37 alc_t6_male_b_se
 # 38 alc_t6_male_b_t
 # 39 alc_t6_male_b-log10p
 # 40 alc_t6_male_c_beta
 # 41 alc_t6_male_c_se
 # 42 alc_t6_male_c_t
 # 43 alc_t6_male_c-log10p
 # 44 alc_t6_male_d_beta
 # 45 alc_t6_male_d_se
 # 46 alc_t6_male_d_t
 # 47 alc_t6_male_d-log10p
 # 48 alc_t6_male_e_beta
 # 49 alc_t6_male_e_se
 # 50 alc_t6_male_e_t
 # 51 alc_t6_male_e-log10p
 # 52 alc_t6_male_f_beta
 # 53 alc_t6_male_f_se
 # 54 alc_t6_male_f_t
 # 55 alc_t6_male_f-log10p

# Unzip gwas & concatenate

for f in UKB23704_AlcCons_6ageTranche_sexStrat_20210201_chr*.gz ; do gunzip -c "$f" > /***/"${f%.*}" ; done

#> [ngillespie@light out]$ for f in UKB23704_AlcCons_6ageTranche_sexStrat_20210201_chr*.gz ; do gunzip -c "$f" > /***/"${f%.*}" ; done
#gzip: UKB23704_AlcCons_6ageTranche_sexStrat_20210201_chr19.gz: unexpected end of file

# Copy over chr19 & use zcat to decompress
gunzip -c /***/UKB23704_AlcCons_6ageTranche_sexStrat_20210201_chr19.gz > /***/UKB23704_AlcCons_6ageTranche_sexStrat_20210201_chr19


# Concatenate
cd /***/
cat UKB23704_AlcCons_6ageTranche_sexStrat_20210201_chr* > UKB23704_AlcCons_6ageTranche_sexStrat_20210201_chr_allPCs_chr1to22

# wc -l UKB23704_AlcCons_6ageTranche_sexStrat_20210201_chr_allPCs_chr1to22

head -1 UKB23704_AlcCons_6ageTranche_sexStrat_20210201_chr_allPCs_chr1to22 | tr ' ' '\n' | awk '{print NR,$1}'	

awk '{print $2,$5,$4,$8,$11}'  *allPCs_chr1to22 > alc_t6_female_a_beta
awk '{print $2,$5,$4,$12,$15}' *allPCs_chr1to22 > alc_t6_female_b_beta
awk '{print $2,$5,$4,$16,$19}' *allPCs_chr1to22 > alc_t6_female_c_beta
awk '{print $2,$5,$4,$20,$23}' *allPCs_chr1to22 > alc_t6_female_d_beta
awk '{print $2,$5,$4,$24,$27}' *allPCs_chr1to22 > alc_t6_female_e_beta
awk '{print $2,$5,$4,$28,$31}' *allPCs_chr1to22 > alc_t6_female_f_beta
awk '{print $2,$5,$4,$32,$35}' *allPCs_chr1to22 > alc_t6_male_a_beta
awk '{print $2,$5,$4,$36,$39}' *allPCs_chr1to22 > alc_t6_male_b_beta
awk '{print $2,$5,$4,$40,$43}' *allPCs_chr1to22 > alc_t6_male_c_beta
awk '{print $2,$5,$4,$44,$47}' *allPCs_chr1to22 > alc_t6_male_d_beta
awk '{print $2,$5,$4,$48,$51}' *allPCs_chr1to22 > alc_t6_male_e_beta
awk '{print $2,$5,$4,$52,$55}' *allPCs_chr1to22 > alc_t6_male_f_beta

# wc -l alc_t6_female_a_beta 	= 16159843
# wc -l alc_t6_male_f_beta 		= 16159843

# GenomicSEM URL 
 # Read in summary statistics using fread from the data.table R package
 # https://github.com/MichelNivard/GenomicSEM/wiki
 # https://github.com/MichelNivard/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects



# Install GenomicSEM 

 install.packages("devtools")
 library(devtools)
 install_github("MichelNivard/GenomicSEM")
 require(GenomicSEM)  
  


# Prepare summary statistics in R

 > nano prepare.R
 > qRh11 prepare.R
 #/home/jpritikin/bin/qR --queue openmx --ncpu 6 /home/ngillespie/GenomicSEM_bmi/temp.R

 setwd("/home/***/")
 require(data.table)
 ln_alc_t6_female_a		<- fread("/***/alc_t6_female_a_beta", data.table=FALSE)
 ln_alc_t6_female_b		<- fread("/***/alc_t6_female_b_beta", data.table=FALSE)
 ln_alc_t6_female_c		<- fread("/***/alc_t6_female_c_beta", data.table=FALSE)
 ln_alc_t6_female_d		<- fread("/***/alc_t6_female_d_beta", data.table=FALSE) 
 ln_alc_t6_female_e		<- fread("/***/alc_t6_female_e_beta", data.table=FALSE)  
 ln_alc_t6_female_f 	<- fread("/***/alc_t6_female_f_beta", data.table=FALSE) 
 ln_alc_t6_male_a		<- fread("/***/alc_t6_male_a_beta", data.table=FALSE)
 ln_alc_t6_male_b		<- fread("/***/alc_t6_male_b_beta", data.table=FALSE)
 ln_alc_t6_male_c		<- fread("/***/alc_t6_male_c_beta", data.table=FALSE)
 ln_alc_t6_male_d		<- fread("/***/alc_t6_male_d_beta", data.table=FALSE) 
 ln_alc_t6_male_e		<- fread("/***/alc_t6_male_e_beta", data.table=FALSE)  
 ln_alc_t6_male_f 		<- fread("/***/alc_t6_male_f_beta", data.table=FALSE)  
     
 names(ln_alc_t6_female_a) 	<- c("SNP","A1","A2","OR","P")
 names(ln_alc_t6_female_b) 	<- c("SNP","A1","A2","OR","P") 
 names(ln_alc_t6_female_c) 	<- c("SNP","A1","A2","OR","P") 
 names(ln_alc_t6_female_d) 	<- c("SNP","A1","A2","OR","P") 
 names(ln_alc_t6_female_e) 	<- c("SNP","A1","A2","OR","P")
 names(ln_alc_t6_female_f) 	<- c("SNP","A1","A2","OR","P")
 names(ln_alc_t6_male_a  ) 	<- c("SNP","A1","A2","OR","P")
 names(ln_alc_t6_male_b  ) 	<- c("SNP","A1","A2","OR","P") 
 names(ln_alc_t6_male_c  ) 	<- c("SNP","A1","A2","OR","P") 
 names(ln_alc_t6_male_d  ) 	<- c("SNP","A1","A2","OR","P") 
 names(ln_alc_t6_male_e  ) 	<- c("SNP","A1","A2","OR","P")
 names(ln_alc_t6_male_f  ) 	<- c("SNP","A1","A2","OR","P")

# Extract only rows with rsIDs    
 ln_alc_t6_female_a			<- ln_alc_t6_female_a[grep("rs", ln_alc_t6_female_a$SNP),]
 ln_alc_t6_female_b			<- ln_alc_t6_female_b[grep("rs", ln_alc_t6_female_b$SNP),]
 ln_alc_t6_female_c			<- ln_alc_t6_female_c[grep("rs", ln_alc_t6_female_c$SNP),]
 ln_alc_t6_female_d			<- ln_alc_t6_female_d[grep("rs", ln_alc_t6_female_d$SNP),]
 ln_alc_t6_female_e			<- ln_alc_t6_female_e[grep("rs", ln_alc_t6_female_e$SNP),]
 ln_alc_t6_female_f			<- ln_alc_t6_female_f[grep("rs", ln_alc_t6_female_f$SNP),]   
 ln_alc_t6_male_a			<- ln_alc_t6_male_a[grep("rs", ln_alc_t6_male_a$SNP),]
 ln_alc_t6_male_b			<- ln_alc_t6_male_b[grep("rs", ln_alc_t6_male_b$SNP),]
 ln_alc_t6_male_c			<- ln_alc_t6_male_c[grep("rs", ln_alc_t6_male_c$SNP),]
 ln_alc_t6_male_d			<- ln_alc_t6_male_d[grep("rs", ln_alc_t6_male_d$SNP),]
 ln_alc_t6_male_e			<- ln_alc_t6_male_e[grep("rs", ln_alc_t6_male_e$SNP),]
 ln_alc_t6_male_f			<- ln_alc_t6_male_f[grep("rs", ln_alc_t6_male_f$SNP),]     

 ln_alc_t6_female_a$SNP 	<- sapply(strsplit(ln_alc_t6_female_a$SNP, ":"), `[`, 1) 
 ln_alc_t6_female_b$SNP 	<- sapply(strsplit(ln_alc_t6_female_b$SNP, ":"), `[`, 1)  
 ln_alc_t6_female_c$SNP 	<- sapply(strsplit(ln_alc_t6_female_c$SNP, ":"), `[`, 1) 
 ln_alc_t6_female_d$SNP 	<- sapply(strsplit(ln_alc_t6_female_d$SNP, ":"), `[`, 1)  
 ln_alc_t6_female_e$SNP 	<- sapply(strsplit(ln_alc_t6_female_e$SNP, ":"), `[`, 1) 
 ln_alc_t6_female_f$SNP 	<- sapply(strsplit(ln_alc_t6_female_f$SNP, ":"), `[`, 1) 
 ln_alc_t6_male_a$SNP 		<- sapply(strsplit(ln_alc_t6_male_a$SNP, ":"), `[`, 1) 
 ln_alc_t6_male_b$SNP 		<- sapply(strsplit(ln_alc_t6_male_b$SNP, ":"), `[`, 1)  
 ln_alc_t6_male_c$SNP 		<- sapply(strsplit(ln_alc_t6_male_c$SNP, ":"), `[`, 1) 
 ln_alc_t6_male_d$SNP 		<- sapply(strsplit(ln_alc_t6_male_d$SNP, ":"), `[`, 1)  
 ln_alc_t6_male_e$SNP 		<- sapply(strsplit(ln_alc_t6_male_e$SNP, ":"), `[`, 1) 
 ln_alc_t6_male_f$SNP 		<- sapply(strsplit(ln_alc_t6_male_f$SNP, ":"), `[`, 1) 
 
 cols.num 			<- c("OR","P")
 ln_alc_t6_female_a[cols.num] <- sapply(ln_alc_t6_female_a[cols.num],as.numeric); str(ln_alc_t6_female_a) ; 
 ln_alc_t6_female_b[cols.num] <- sapply(ln_alc_t6_female_b[cols.num],as.numeric); str(ln_alc_t6_female_b) ; 
 ln_alc_t6_female_c[cols.num] <- sapply(ln_alc_t6_female_c[cols.num],as.numeric); str(ln_alc_t6_female_c) ; 
 ln_alc_t6_female_d[cols.num] <- sapply(ln_alc_t6_female_d[cols.num],as.numeric); str(ln_alc_t6_female_d) ; 
 ln_alc_t6_female_e[cols.num] <- sapply(ln_alc_t6_female_e[cols.num],as.numeric); str(ln_alc_t6_female_e) ; 
 ln_alc_t6_female_f[cols.num] <- sapply(ln_alc_t6_female_f[cols.num],as.numeric); str(ln_alc_t6_female_f) ; 
 ln_alc_t6_male_a[cols.num]   <- sapply(ln_alc_t6_male_a[cols.num],as.numeric); str(ln_alc_t6_male_a) ; 
 ln_alc_t6_male_b[cols.num]   <- sapply(ln_alc_t6_male_b[cols.num],as.numeric); str(ln_alc_t6_male_b) ; 
 ln_alc_t6_male_c[cols.num]   <- sapply(ln_alc_t6_male_c[cols.num],as.numeric); str(ln_alc_t6_male_c) ; 
 ln_alc_t6_male_d[cols.num]   <- sapply(ln_alc_t6_male_d[cols.num],as.numeric); str(ln_alc_t6_male_d) ; 
 ln_alc_t6_male_e[cols.num]   <- sapply(ln_alc_t6_male_e[cols.num],as.numeric); str(ln_alc_t6_male_e) ; 
 ln_alc_t6_male_f[cols.num]   <- sapply(ln_alc_t6_male_f[cols.num],as.numeric); str(ln_alc_t6_male_f) ; 

 ln_alc_t6_female_a$P <- 10^(-1*ln_alc_t6_female_a$P)
 ln_alc_t6_female_b$P <- 10^(-1*ln_alc_t6_female_b$P)
 ln_alc_t6_female_c$P <- 10^(-1*ln_alc_t6_female_c$P)
 ln_alc_t6_female_d$P <- 10^(-1*ln_alc_t6_female_d$P)  
 ln_alc_t6_female_e$P <- 10^(-1*ln_alc_t6_female_e$P)
 ln_alc_t6_female_f$P <- 10^(-1*ln_alc_t6_female_f$P)
 ln_alc_t6_male_a$P   <- 10^(-1*ln_alc_t6_male_a$P)
 ln_alc_t6_male_b$P   <- 10^(-1*ln_alc_t6_male_b$P)
 ln_alc_t6_male_c$P   <- 10^(-1*ln_alc_t6_male_c$P)
 ln_alc_t6_male_d$P   <- 10^(-1*ln_alc_t6_male_d$P)  
 ln_alc_t6_male_e$P   <- 10^(-1*ln_alc_t6_male_e$P)
 ln_alc_t6_male_f$P   <- 10^(-1*ln_alc_t6_male_f$P)

 #output the result in a .txt file called scz_withRS.txt#
 write.table(ln_alc_t6_female_a, file = "ln_alc_t6_female_a_withRS.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
 write.table(ln_alc_t6_female_b, file = "ln_alc_t6_female_b_withRS.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE) 
 write.table(ln_alc_t6_female_c, file = "ln_alc_t6_female_c_withRS.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
 write.table(ln_alc_t6_female_d, file = "ln_alc_t6_female_d_withRS.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
 write.table(ln_alc_t6_female_e, file = "ln_alc_t6_female_e_withRS.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
 write.table(ln_alc_t6_female_f, file = "ln_alc_t6_female_f_withRS.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
 write.table(ln_alc_t6_male_a  , file   = "ln_alc_t6_male_a_withRS.txt"  , sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
 write.table(ln_alc_t6_male_b  , file   = "ln_alc_t6_male_b_withRS.txt"  , sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE) 
 write.table(ln_alc_t6_male_c  , file   = "ln_alc_t6_male_c_withRS.txt"  , sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
 write.table(ln_alc_t6_male_d  , file   = "ln_alc_t6_male_d_withRS.txt"  , sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
 write.table(ln_alc_t6_male_e  , file   = "ln_alc_t6_male_e_withRS.txt"  , sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)
 write.table(ln_alc_t6_male_f  , file   = "ln_alc_t6_male_f_withRS.txt"  , sep = "\t", quote = FALSE, row.names = FALSE, col.names=TRUE)

 save.image("/***/GenomicSEM_alc.RData")    
 
 less prepare.Rout
 
 
 
 
 
# Munge summary statistics files
 
  # https://github.com/MichelNivard/GenomicSEM/wiki/6.-GenomicSEM-and-OpenMx
   # note that we include the argument names (e.g., files, hm3) for
   # completeness in the code below, but that this is not necessary to run the code.
   
 nano munge_alcohol.R
 qRh11 munge_alcohol.R

  library(GenomicSEM)
  setwd("/***/")
    munge(files=c(
 	   "ln_alc_t6_female_a_withRS.txt",
 	   "ln_alc_t6_female_b_withRS.txt",
 	   "ln_alc_t6_female_c_withRS.txt",
 	   "ln_alc_t6_female_d_withRS.txt",
 	   "ln_alc_t6_female_e_withRS.txt",
 	   "ln_alc_t6_female_f_withRS.txt",
 	   "ln_alc_t6_male_a_withRS.txt",  
	   "ln_alc_t6_male_b_withRS.txt",  
	   "ln_alc_t6_male_c_withRS.txt",  
	   "ln_alc_t6_male_d_withRS.txt",  
	   "ln_alc_t6_male_e_withRS.txt",  
	   "ln_alc_t6_male_f_withRS.txt"),
    hm3="/***/w_hm3.snplist",
    trait.names=c(
	 	"alc_a_f","alc_b_f","alc_c_f","alc_d_f","alc_e_f","alc_f_f",
 	   	"alc_a_m","alc_b_m","alc_c_m","alc_d_m","alc_e_m","alc_f_m"), 
    N=c(13694, 19124, 23081, 26934, 35835, 24932,
		13272, 17017, 20062, 24801, 35720, 30312), 
 	info.filter = 0.9, 
 	maf.filter = 0.01)
	
	
	


# Genetic covariance matrices


# Alcohol - women

alc_corrs_f.R
qRh11 alc_corrs_f.R

  library(GenomicSEM)
  setwd("/***/")
  traits 			<- c("alc_a_f.sumstats.gz",
  	 					 "alc_b_f.sumstats.gz",
  	 					 "alc_c_f.sumstats.gz",
  	 					 "alc_d_f.sumstats.gz",
  	 					 "alc_e_f.sumstats.gz",
  	 					 "alc_f_f.sumstats.gz")
  sample.prev 		<- c(NA,NA,NA,NA,NA,NA)
  population.prev 	<- c(NA,NA,NA,NA,NA,NA)
  ld 				<- "/***/eur_w_ld_chr"
  wld 				<- "/***/eur_w_ld_chr"
  trait.names		<- c("ALC1_f","ALC2_f","ALC3_f","ALC4_f","ALC5_f","ALC6_f")
  LDSCoutput 		<- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
  save(LDSCoutput, file="UKB23704_AlcCons_6ageTranche_sexStrat_20210201_f.RData")
  cov2cor(LDSCoutput$S)
  OpenMx::diag2vec(LDSCoutput$S)  
  
  load("UKB23704_AlcCons_6ageTranche_sexStrat_20210201_f.RData")
  S2 <- as.matrix(Matrix::nearPD(LDSCoutput$S)$mat)
  round(cov2cor(S2),4)
  
        ALC1_f ALC2_f ALC3_f ALC4_f ALC5_f ALC6_f
 ALC1_f 1.0000 0.9406 0.9675 0.8018 0.8681 0.8549
 ALC2_f 0.9406 1.0000 0.8742 0.6627 0.6702 0.6349
 ALC3_f 0.9675 0.8742 1.0000 0.6815 0.8288 0.9126
 ALC4_f 0.8018 0.6627 0.6815 1.0000 0.9435 0.7436
 ALC5_f 0.8681 0.6702 0.8288 0.9435 1.0000 0.9220
 ALC6_f 0.8549 0.6349 0.9126 0.7436 0.9220 1.0000
 
       ALC1_m ALC2_m ALC3_m ALC4_m ALC5_m ALC6_m
ALC1_m 1.0000 0.9396 0.9317 0.9669 0.8512 0.9821
ALC2_m 0.9396 1.0000 0.9472 0.8279 0.8031 0.9872
ALC3_m 0.9317 0.9472 1.0000 0.8809 0.9507 0.9520
ALC4_m 0.9669 0.8279 0.8809 1.0000 0.8727 0.9045
ALC5_m 0.8512 0.8031 0.9507 0.8727 1.0000 0.8334
ALC6_m 0.9821 0.9872 0.9520 0.9045 0.8334 1.000
 
 
 

# Alcohol - men
 
alc_corrs_m.R
qRh11 alc_corrs_m.R

  library(GenomicSEM)
  setwd("/***/")
  traits 			<- c("alc_a_m.sumstats.gz",
  	 				     "alc_b_m.sumstats.gz",
  	 				     "alc_c_m.sumstats.gz",
  	 				     "alc_d_m.sumstats.gz",
  	 				     "alc_e_m.sumstats.gz",
  	 				     "alc_f_m.sumstats.gz")
  sample.prev 		<- c(NA,NA,NA,NA,NA,NA)
  population.prev 	<- c(NA,NA,NA,NA,NA,NA)
  ld 				<- "/***/eur_w_ld_chr"
  wld 				<- "/***/eur_w_ld_chr"
  trait.names		<- c("ALC1_m","ALC2_m","ALC3_m","ALC4_m","ALC5_m","ALC6_m")
  LDSCoutput 		<- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
  save(LDSCoutput, file="UKB23704_AlcCons_6ageTranche_sexStrat_20210201_m.RData")
  cov2cor(LDSCoutput$S)
  
  load("UKB23704_AlcCons_6ageTranche_sexStrat_20210201_m.RData")
  S2 <- as.matrix(Matrix::nearPD(LDSCoutput$S)$mat)
  round(cov2cor(S2),4)
  
    OpenMx::diag2vec(LDSCoutput$S)  
  
         ALC1_m ALC2_m ALC3_m ALC4_m ALC5_m ALC6_m
  ALC1_m 1.0000 0.9396 0.9317 0.9669 0.8512 0.9821
  ALC2_m 0.9396 1.0000 0.9472 0.8279 0.8031 0.9872
  ALC3_m 0.9317 0.9472 1.0000 0.8809 0.9507 0.9520
  ALC4_m 0.9669 0.8279 0.8809 1.0000 0.8727 0.9045
  ALC5_m 0.8512 0.8031 0.9507 0.8727 1.0000 0.8334
  ALC6_m 0.9821 0.9872 0.9520 0.9045 0.8334 1.0000
  
  
  
  
# Alcohol - men + women
 
 alc_corrs_mf.R
 qRh11 alc_corrs_mf.R

   library(GenomicSEM)
   setwd("/***/")
   traits 			<- c("alc_a_m.sumstats.gz",
	                     "alc_b_m.sumstats.gz",
	                     "alc_c_m.sumstats.gz",
	                     "alc_d_m.sumstats.gz",
	                     "alc_e_m.sumstats.gz",
	                     "alc_f_m.sumstats.gz",
	   				     "alc_a_f.sumstats.gz",
   	 					 "alc_b_f.sumstats.gz",
   	 					 "alc_c_f.sumstats.gz",
   	 					 "alc_d_f.sumstats.gz",
   	 					 "alc_e_f.sumstats.gz",
   	 					 "alc_f_f.sumstats.gz")
   sample.prev 		<- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
   population.prev 	<- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
   ld 				<- "/***/eur_w_ld_chr"
   wld 				<- "/***/eur_w_ld_chr"
   trait.names		<- c("ALC1_m","ALC2_m","ALC3_m","ALC4_m","ALC5_m","ALC6_m",
   						 "ALC1_f","ALC2_f","ALC3_f","ALC4_f","ALC5_f","ALC6_f")
   LDSCoutput 		<- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
   save(LDSCoutput, file="UKB23704_AlcCons_6ageTranche_sexStrat_20210201_mf.RData")
   cov2cor(LDSCoutput$S)
   
   load("***/UKB23704_AlcCons_6ageTranche_sexStrat_20210201_mf.RData")
   S2 <- as.matrix(Matrix::nearPD(LDSCoutput$S)$mat)
   round(cov2cor(S2),4)
   
        ALC1_m ALC2_m ALC3_m ALC4_m ALC5_m ALC6_m  ALC1_f ALC2_f ALC3_f ALC4_f ALC5_f ALC6_f
  ALC1_m 1.0000 0.8043 0.7882 0.8184 0.7451 0.8611  0.6956 0.4153 0.5613 0.7007 0.5911 0.3435
  ALC2_m 0.8043 1.0000 0.8492 0.8020 0.7888 0.9158  0.8024 0.7411 0.7802 0.6896 0.8073 0.6240
  ALC3_m 0.7882 0.8492 1.0000 0.7244 0.8597 0.9210  0.9151 0.6109 0.7685 0.5569 0.6664 0.7305
  ALC4_m 0.8184 0.8020 0.7244 1.0000 0.8855 0.7929  0.5957 0.4611 0.7995 0.8028 0.8371 0.6253
  ALC5_m 0.7451 0.7888 0.8597 0.8855 1.0000 0.7882  0.6858 0.4858 0.8491 0.6462 0.7924 0.7186
  ALC6_m 0.8611 0.9158 0.9210 0.7929 0.7882 1.0000  0.8331 0.6034 0.6860 0.6128 0.6539 0.6488
  
  ALC1_f 0.6956 0.8024 0.9151 0.5957 0.6858 0.8331  1.0000 0.8047 0.7561 0.6438 0.6642 0.7243
  ALC2_f 0.4153 0.7411 0.6109 0.4611 0.4858 0.6034  0.8047 1.0000 0.7768 0.5986 0.6543 0.6049
  ALC3_f 0.5613 0.7802 0.7685 0.7995 0.8491 0.6860  0.7561 0.7768 1.0000 0.6562 0.8284 0.7635
  ALC4_f 0.7007 0.6896 0.5569 0.8028 0.6462 0.6128  0.6438 0.5986 0.6562 1.0000 0.8742 0.5981
  ALC5_f 0.5911 0.8073 0.6664 0.8371 0.7924 0.6539  0.6642 0.6543 0.8284 0.8742 1.0000 0.7834
  ALC6_f 0.3435 0.6240 0.7305 0.6253 0.7186 0.6488  0.7243 0.6049 0.7635 0.5981 0.7834 1.0000
 
 
 
 
 
# Alcohol & bmi - men
  
  alc_bmi_corrs_m.R
  qRh11 alc_bmi_corrs_m.R
  less alc_bmi_corrs_m.Rout

    library(GenomicSEM)
    setwd("/***/")
    traits 			<- c(    "alc_a_m.sumstats.gz",
    	 				     "alc_b_m.sumstats.gz",
    	 				     "alc_c_m.sumstats.gz",
    	 				     "alc_d_m.sumstats.gz",
    	 				     "alc_e_m.sumstats.gz",
    	 				     "alc_f_m.sumstats.gz",
							 "/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_a_m.sumstats.gz",
							 "/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_b_m.sumstats.gz",
							 "/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_c_m.sumstats.gz",
							 "/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_d_m.sumstats.gz",
							 "/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_e_m.sumstats.gz",
							 "/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_f_m.sumstats.gz")
    sample.prev 		<- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
    population.prev 	<- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
    ld 				<- "/***/eur_w_ld_chr"
    wld 				<- "/***/eur_w_ld_chr"
    trait.names		<- c("ALC1_m","ALC2_m","ALC3_m","ALC4_m","ALC5_m","ALC6_m",
						 "BMI1_m","BMI2_m","BMI3_m","BMI4_m","BMI5_m","BMI6_m")
    LDSCoutput 		<- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
    save(LDSCoutput, file="UKB23704_AlcCons_BMI_6ageTranche_sexStrat_20210201_m.RData")
    cov2cor(LDSCoutput$S)
  
    load("UKB23704_AlcCons_BMI_6ageTranche_sexStrat_20210201_m.RData")
    S2 <- as.matrix(Matrix::nearPD(LDSCoutput$S)$mat)
    round(cov2cor(S2),4)
	
	
    load("UKB23704_AlcCons_BMI_6ageTranche_sexStrat_20210201_m.RData")
    S2 <- as.matrix(Matrix::nearPD(LDSCoutput$S)$mat)
    round(cov2cor(S2),4)
	
        ALC1_m ALC2_m ALC3_m ALC4_m  ALC5_m  ALC6_m  BMI1_m  BMI2_m  BMI3_m BMI4_m BMI5_m BMI6_m
 ALC1_m 1.0000 0.8362 0.8301 0.8198  0.7912  0.8488  0.1193  0.0943  0.1871 0.2108 0.2228 0.2774
 ALC2_m 0.8362 1.0000 0.8453 0.7807  0.6702  0.9280  0.0550  0.0932  0.2801 0.2381 0.1334 0.1406
 ALC3_m 0.8301 0.8453 1.0000 0.8504  0.9451  0.9060  0.0159  0.0623  0.0909 0.1805 0.1447 0.0804
 ALC4_m 0.8198 0.7807 0.8504 1.0000  0.7740  0.8983  0.1058  0.2408  0.2385 0.2785 0.2305 0.2292
 ALC5_m 0.7912 0.6702 0.9451 0.7740  1.0000  0.8058 -0.0160 -0.0418 -0.0681 0.0544 0.0352 0.0308
 ALC6_m 0.8488 0.9280 0.9060 0.8983  0.8058  1.0000 -0.0678  0.0361  0.1012 0.1347 0.0033 0.0640
 BMI1_m 0.1193 0.0550 0.0159 0.1058 -0.0160 -0.0678  1.0000  0.8640  0.9141 0.8805 0.8084 0.9128
 BMI2_m 0.0943 0.0932 0.0623 0.2408 -0.0418  0.0361  0.8640  1.0000  0.8966 0.9779 0.8950 0.9276
 BMI3_m 0.1871 0.2801 0.0909 0.2385 -0.0681  0.1012  0.9141  0.8966  1.0000 0.9288 0.8346 0.8781
 BMI4_m 0.2108 0.2381 0.1805 0.2785  0.0544  0.1347  0.8805  0.9779  0.9288 1.0000 0.9134 0.9370
 BMI5_m 0.2228 0.1334 0.1447 0.2305  0.0352  0.0033  0.8084  0.8950  0.8346 0.9134 1.0000 0.8726
 BMI6_m 0.2774 0.1406 0.0804 0.2292  0.0308  0.0640  0.9128  0.9276  0.8781 0.9370 0.8726 1.0000
	
	
	

# Alcohol & bmi - women
  
 nano alc_bmi_corrs_f.R
 qRh11 alc_bmi_corrs_f.R

 install.packages("GenomicSEM")
  library(GenomicSEM)
  setwd("/***/")
  traits 			<- c("alc_a_f.sumstats.gz",
  	 				     "alc_b_f.sumstats.gz",
  	 				     "alc_c_f.sumstats.gz",
  	 				     "alc_d_f.sumstats.gz",
  	 				     "alc_e_f.sumstats.gz",
  	 				     "alc_f_f.sumstats.gz",
						 "/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_a_f.sumstats.gz",
						 "/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_b_f.sumstats.gz",
						 "/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_c_f.sumstats.gz",
						 "/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_d_f.sumstats.gz",
						 "/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_e_f.sumstats.gz",
						 "/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_f_f.sumstats.gz")
  sample.prev 		<- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  population.prev 	<- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  ld 				<- "/***/eur_w_ld_chr"
  wld 				<- "/***/eur_w_ld_chr"
  trait.names		<- c("ALC1_f","ALC2_f","ALC3_f","ALC4_f","ALC5_f","ALC6_f",
						 "BMI1_f","BMI2_f","BMI3_f","BMI4_f","BMI5_f","BMI6_f")
  LDSCoutput 		<- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
  save(LDSCoutput, file="UKB23704_AlcCons_BMI_6ageTranche_sexStrat_20210201_f.RData")
  cov2cor(LDSCoutput$S)
 
    load("/Volumes/ngillespie/Documents/work/projects/2020/2020 2. genomicSEM_bmi/4. Alcohol_GWAS_results_40_PCs_separate_MFs/results/UKB23704_AlcCons_BMI_6ageTranche_sexStrat_20210201_f.RData")
    S2 <- as.matrix(Matrix::nearPD(LDSCoutput$S)$mat)
    round(cov2cor(S2),4)
	
       ALC1_f ALC2_f  ALC3_f  ALC4_f ALC5_f  ALC6_f  BMI1_f  BMI2_f BMI3_f BMI4_f BMI5_f  BMI6_f
ALC1_f 1.0000 0.9102  0.8068  0.7582 0.8176  0.8393  0.1262  0.1776 0.2263 0.1327 0.1674  0.0958
ALC2_f 0.9102 1.0000  0.7609  0.6289 0.6609  0.6174  0.1668  0.2138 0.2686 0.2357 0.1295  0.1034
ALC3_f 0.8068 0.7609  1.0000  0.6711 0.7636  0.6176  0.1356  0.0329 0.1409 0.1465 0.0834 -0.0006
ALC4_f 0.7582 0.6289  0.6711  1.0000 0.8971  0.6231  0.0302 -0.0455 0.0509 0.0040 0.0273 -0.0528
ALC5_f 0.8176 0.6609  0.7636  0.8971 1.0000  0.8463  0.0329  0.0561 0.1161 0.1477 0.1165  0.0827
ALC6_f 0.8393 0.6174  0.6176  0.6231 0.8463  1.0000 -0.0318  0.1374 0.1264 0.0917 0.1427  0.1199
BMI1_f 0.1262 0.1668  0.1356  0.0302 0.0329 -0.0318  1.0000  0.9268 0.9635 0.8957 0.9285  0.8891
BMI2_f 0.1776 0.2138  0.0329 -0.0455 0.0561  0.1374  0.9268  1.0000 0.9850 0.9312 0.9230  0.9352
BMI3_f 0.2263 0.2686  0.1409  0.0509 0.1161  0.1264  0.9635  0.9850 1.0000 0.9380 0.9134  0.9050
BMI4_f 0.1327 0.2357  0.1465  0.0040 0.1477  0.0917  0.8957  0.9312 0.9380 1.0000 0.8724  0.9197
BMI5_f 0.1674 0.1295  0.0834  0.0273 0.1165  0.1427  0.9285  0.9230 0.9134 0.8724 1.0000  0.9773
BMI6_f 0.0958 0.1034 -0.0006 -0.0528 0.0827  0.1199  0.8891  0.9352 0.9050 0.9197 0.9773  1.0000
	
	
	
#alcohol
cd /***/

alc_a_f.sumstats.gz
alc_a_m.sumstats.gz
alc_b_f.sumstats.gz
alc_b_m.sumstats.gz
alc_c_f.sumstats.gz
alc_c_m.sumstats.gz
alc_d_f.sumstats.gz
alc_d_m.sumstats.gz
alc_e_f.sumstats.gz
alc_e_m.sumstats.gz
alc_f_f.sumstats.gz
alc_f_m.sumstats.gz

#bmi
cd /home/ngillespie/GenomicSEM_bmi_40PCs/

bmi_a_f.sumstats.gz
bmi_a_m.sumstats.gz
bmi_b_f.sumstats.gz
bmi_b_m.sumstats.gz
bmi_c_f.sumstats.gz
bmi_c_m.sumstats.gz
bmi_d_f.sumstats.gz
bmi_d_m.sumstats.gz
bmi_e_f.sumstats.gz
bmi_e_m.sumstats.gz
bmi_f_f.sumstats.gz
bmi_f_m.sumstats.gz


# GWAS by subtraction
 # https://rpubs.com/MichelNivard/565885
 
# Prepare the data 
 (AlcBmi_v1.R)
 library(GenomicSEM)
 setwd("/home/ngillespie/GenomicSEM_bmi_alc/") 
 traits 			<- c("/***/alc_a_f.sumstats.gz","/home/ngillespie/GenomicSEM_bmi_40PCs/bmi_a_f.sumstats.gz")
 sample.prev 		<- c(NA,NA)
 population.prev	<- c(NA,NA)
 ld					<- "/***/eur_w_ld_chr"
 wld				<- "/***/eur_w_ld_chr"
 trait.names		<-c("ALC", "BMI")
 LDSCoutput 		<- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
 save(LDSCoutput, file="LDSCoutputAlcBmi.RData")
 cov2cor(LDSCoutput$S)
 #             ALC        BMI
 # [1,] 1.00000000 0.08173671
 # [2,] 0.08173671 1.00000000



# Run SEM model without SNPs

library(lavaan)
load(file="/home/ngillespie/GenomicSEM_bmi_alc/LDSCoutputAlcBmi.RData")
model <-'A=~NA*BMI + start(0.4)*ALC
        B=~NA*BMI
# Latent variances & covariance		
		B ~~ 1*B
		A ~~ 1*A
		A ~~ 0*B
# Phenotypic variances & covariance		
		ALC ~~ 0*BMI
		ALC ~~ 0*ALC
		BMI ~~ 0*BMI'
output	<-usermodel(LDSCoutput,estimation="DWLS",model=model)
save(output, file="/home/ngillespie/GenomicSEM_bmi_alc/Modeloutput.Rdata" )

load("/home/ngillespie/GenomicSEM_bmi_alc/Modeloutput.Rdata")
  # $modelfit
  #    chisq df p_chisq AIC CFI         SRMR
  # df    NA  0      NA  NA   1 1.956288e-09
  # 
  # $results
  #   lhs op rhs Unstand_Est         Unstand_SE STD_Genotype    STD_Genotype_SE
  # 1   A =~ BMI  0.04153037   0.09595010764072   0.08173672  0.188841233839844
  # 2   A =~ ALC  0.20881978 0.0747025564479532   1.00000000  0.357736975586839
  # 3   B =~ BMI  0.50639922 0.0278216372344786   0.99665396 0.0547562941013872
  # 4   B ~~   B  1.00000000                      1.00000000                   
  # 5   A ~~   A  1.00000000                      1.00000000                   
  # 6   A ~~   B  0.00000000                      0.00000000                   
  # 7 BMI ~~ ALC  0.00000000                      0.00000000                   
  # 8 ALC ~~ ALC  0.00000000                      0.00000000                   
  # 9 BMI ~~ BMI  0.00000000                      0.00000000                   
  #      STD_All
  # 1 0.08173672
  # 2 1.00000000
  # 3 0.99665396
  # 4 1.00000000
  # 5 1.00000000
  # 6 0.00000000
  # 7 0.00000000
  # 8 0.00000000
  # 9 0.00000000






# Check if GWAS sumstat p-values are standard error of logistic beta

 # https://github.com/GenomicSEM/GenomicSEM/wiki/2.-Important-resources-and-key-information
 # Use decision tree if analyzing binary data

  head /***/UKB23704_AlcCons_6ageTranche_sexStrat_20210201_chr_allPCs_chr1to22
  
  awk '{print $2,$5,$4,$8,$11,$9}' /***/UKB23704_AlcCons_6ageTranche_sexStrat_20210201_chr_allPCs_chr1to22 | head -100  > /***/alc_100.txt
  head /***/alc_100.txt
  df <- read.table("/***/alc_100.txt",header=T)
  head(df)
  names(df) <- c("SNP","A1","A2","OR","P", "SE")
   #             SNP A1 A2         OR        P       SE
   # 1 1:692794_CA_C  C CA -0.0101230 0.313080 0.014541
   # 2    rs12238997  G  A -0.0163500 0.627130 0.013796
   # 3   rs371890604  C  G -0.0089355 0.249410 0.015453
   # 4   rs149887893  C  T -0.0065897 0.098181 0.025703
   # 5    rs12184267  T  C -0.0041274 0.065760 0.023315
   # 6    rs12184277  G  A -0.0037832 0.060165 0.023230
  df$OR <-  as.numeric(df$OR)
  df$SE <-  as.numeric(df$SE)
  df$temp <- 2*pnorm(abs(log(df$OR))/df$SE,lower.tail=FALSE)
  head(df)
  
  2*pnorm(abs(log(0.56167000))/0.440810,lower.tail=FALSE) = 0.1906723
  
  	# If df$temp = df$P, then it's likely to be the SE of a logistic beta. In many cases it is the standard error of a logistic beta even when the effect column says i 
  	# says it's an OR. 
	
	# alcohol and BMI are continuous.


# Clean the SNPs for GWAS

(CleanSNPs.R)
 library(GenomicSEM)
 library(lavaan)
 setwd("/home/ngillespie/GenomicSEM_bmi_alc")
 files 			= c("/***/ln_alc_t6_female_a_withRS.txt", "/home/ngillespie/GenomicSEM_bmi_40PCs/ln_bmi_t6_female_a_withRS.txt")
 ref 			= "/home/ngillespie/GenomicSEM_bmi_alc/reference.1000G.maf.0.005.txt"
 trait.names 	= c("ALC","BMI")
 se.logit 		= c(F,F)
 info.filter 	= 0.6
 maf.filter 	= 0.01
 p_sumstats		<- sumstats(files, ref, trait.names, se.logit, info.filter, maf.filter, OLS=c(T,T),linprob=NULL, prop=NULL, N=c(13694,18110))
 save(p_sumstats, file="/home/ngillespie/GenomicSEM_bmi_alc/Sumstats.RData")
load("/home/ngillespie/GenomicSEM_bmi_alc/Sumstats.RData")

	# > head(p_sumstats)
	#           SNP CHR     BP       MAF A1 A2    beta.ALC     se.ALC    beta.BMI
	# 1 rs371890604   1 707522 0.0994036  G  C 0.011677710 0.02019546  0.01958800
	# 2  rs12184267   1 715265 0.0437376  C  T 0.005230457 0.02954635  0.05430986
	# 3  rs12184277   1 715367 0.0437376  A  G 0.004811765 0.02954635  0.05264064
	# 4 rs116801199   1 720381 0.0447316  G  T 0.005098111 0.02923142  0.05496437
	# 5  rs12565286   1 721290 0.0457256  G  C 0.004797786 0.02892701  0.05260029
	# 6   rs2977670   1 723891 0.0516900  G  C 0.007782306 0.02729238 -0.05638639
	#       se.BMI
	# 1 0.01756142
	# 2 0.02569270
	# 3 0.02569270
	# 4 0.02541885
	# 5 0.02515414
	# 6 0.02373271
	# > dim(p_sumstats)
	# [1] 7875957      10

# Run the GWAS-by-subtraction
 # We are now ready to run the GenomicSEM GWAS-by-subtraction. Additional to the model ran above, we introduce the SNP variable. The SNP is regressed on the latent variables C and NC concurrently ( e.g. C~SNP and N~CNP).
 # model with SNP
model<-'C=~NA*EA +start(0.2)*EA + start(0.5)*CP
         NC=~NA*EA +start(0.2)*EA
         
         C~SNP
         NC~SNP

         NC~~1*NC
         C~~1*C
         C~~0*NC

         CP ~~ 0*EA
         CP~~0*CP
         EA~~0*EA
         SNP~~SNP'
outputGWAS<-userGWAS(covstruc=LDSCoutput,SNPs=p_sumstats,estimation="DWLS",model=model,sub =c("C~SNP","NC~SNP"))# printwarn = FALSE

