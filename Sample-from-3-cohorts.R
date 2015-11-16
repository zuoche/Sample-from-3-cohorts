setwd("H:/plink-1.07-dos/")

eira_icfam <- read.table("ICHIP//cleandata_4.fam")
eira_hlafam <- read.table("SE-E.fam")
eira_key <- read.table("key_ICHIP.txt")
tmp <- merge(eira_key, eira_hlafam, by="V2")
names(eira_icfam)[2] <- "iid"
names(tmp)[4] <- "iid"
tmp <- merge(tmp, eira_icfam, by="iid") #EIRA:IC+HLA(found in key)=4618

test <- eira_icfam[,2] %in% eira_key[,4]
summary(test)   #4815/4830 found in key
test <- eira_hlafam[,2] %in% eira_key[,2]
summary(test)#4702/4687 found in key

#Convert sample id (4815 ids can be converted)
eira_icnew <- eira_icfam
names(eira_key)[1:2] <- c("fid","iid")
names(eira_icnew)[1:2] <- c("fid","iid")
key_fid <- sapply(eira_key$fid, as.character)
key_iid <- sapply(eira_key$iid, as.character)
fid <- sapply(eira_icnew$fid, as.character)
iid <- sapply(eira_icnew$iid, as.character)
index <- match(iid, eira_key[,4], nomatch=0)
#fid[index!=0] <- key_fid[index]
iid[index!=0] <- key_iid[index]
#eira_icnew <- cbind(eira_icnew[,1], iid, eira_icnew[,3:6])
#Using same fid as HLA .bim
hla_fid <- sapply(eira_hlafam[,1], as.character)
index <- match(iid, eira_hlafam[,2], nomatch=0)
fid[index!=0] <- hla_fid[index]
eira_icnew <- cbind(fid,iid,eira_icnew[,3:6])
write.table(eira_icnew, "eira_ichip_newID.fam", quote=F, row.names=F, col.names=F, sep="\t")

#Select individuals with ab data
wtccc_icfam <- read.table("Experimental data-selected/iChip_RACI_PhaseII_UK_QCgp.fam")
wtccc_hlafam <- read.table("UK.fam")
narac_icfam <- read.table("Experimental data-selected/iChip_RACI_PhaseII_US_QCgp.fam")
narac_hlafam <- read.table("US.fam")

pheno <- read.table("Merged_pheno1-35_case-ctrl.txt")
test <- pheno[,2] %in% wtccc_hlafam[,2] #10470/12300 with ab data
test <- pheno[,2] %in% narac_hlafam[,2] #4042/4325 with ab data
test <- pheno[,2] %in% tmp[,2] #4325/4618 with ab data
cohort <- pheno[which(pheno[,2] %in% wtccc_hlafam[,2] | pheno[,2] %in% narac_hlafam[,2] | pheno[,2] %in% tmp[,2]),1:2]
cov_eira <- rep(NA, 18837)
cov_wtccc <- rep(NA, 18837)
cov_narac <- rep(NA, 18837)
for (i in 1:18837) cov_eira[i] <- ifelse(cohort[i,2] %in% tmp[,2], 2, 1)
for (i in 1:18837) cov_wtccc[i] <- ifelse(cohort[i,2] %in% wtccc_hlafam[,2], 2, 1)
for (i in 1:18837) cov_narac[i] <- ifelse(cohort[i,2] %in% narac_hlafam[,2], 2, 1)
cohort <- cbind(cohort, cov_eira, cov_wtccc, cov_narac)
write.table(cohort, "3Cohorts_dummy_covar.txt", quote=F, row.names=F, col.names=F, sep="\t")


#Convert ichip- marker names to rs- names
#Memory space not enough!
wtccc_bim <- read.table("Experimental data-selected/iChip_RACI_PhaseII_UK_QCgp.bim")
narac_bim <- read.table("Experimental data-selected/iChip_RACI_PhaseII_US_QCgp.bim")
eira_bim <- read.table("ICHIP//cleandata_4.bim")
ic_key <- read.table("ichip.hg18.hg19.dbsnpID.chr1-26.txt/ichip.hg18.hg19.dbsnpID.chr1-26.txt", header=T)
wtccc_ic <- sapply(wtccc_bim[,2], as.character)
