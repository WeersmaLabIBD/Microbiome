
#####################################################################################################################
#####################################################gender check####################################################
#####################################################################################################################

lld_gender=read.table(file = "LLD_gender_info.txt",header = F,sep = "\t")
lld_plink=read.table(file = "LLD_checksex.sexcheck",header = T,as.is=T)
colnames(lld_gender)=c("IID","gender")
gender_check=merge(lld_plink,lld_gender,all = T,by="IID")

male=subset(gender_check,gender_check$gender=="Male")
female=subset(gender_check,gender_check$gender=="Female")
no_sex=subset(gender_check,gender_check$SNPSEX==0)
male_problem=subset(male,male$F<0.4)
femal_problem=subset(female,female$F>0.7)
lld_problematic=rbind(male_problem,femal_problem)

write.table(lld_problematic,file = "LLD_sex_problmatic.txt",quote = F,row.names = F,sep = "\t")

jpeg(res=300, width=10, height=5, units="in")
par(mfcol=c(1,2), pty="m")
hist(male[,6], main="Male", ylab="Freq", xlab="Chr X inbreeding estimate", breaks =100)
hist(female[,6], main="Female", ylab="Freq", xlab="Chr X inbreeding estimate", breaks=100)
dev.off()

IBD_gender=read.table(file = "IBD_gender_info.txt",header = F,sep = "\t")
IBD_plink=read.table(file = "IBD_checksex.sexcheck",header = T,as.is=T)
colnames(IBD_gender)=c("IID","gender")
gender_check_ibd=merge(IBD_plink,IBD_gender,all = T,by="IID")

male_ibd=subset(gender_check_ibd,gender_check_ibd$gender=="Male")
female_ibd=subset(gender_check_ibd,gender_check_ibd$gender=="Female")
no_sex_ibd=subset(gender_check_ibd,gender_check_ibd$SNPSEX==0)
male_ibd_problem=subset(male_ibd,male_ibd$F<0.4)
femal_ibd_problem=subset(female_ibd,female_ibd$F>0.7)
ibd_problematic=rbind(male_ibd_problem,femal_ibd_problem)

write.table(ibd_problematic,file = "IBD_sex_problematic.txt",quote = F,row.names = F)

jpeg(res=300, width=10, height=5, units="in")
par(mfcol=c(2,2), pty="m")
hist(male[,6], main="Male", ylab="Freq", xlab="Chr X inbreeding estimate", breaks =100)
hist(female[,6], main="Female", ylab="Freq", xlab="Chr X inbreeding estimate", breaks=100)
hist(male_ibd[,6], main="Male", ylab="Freq", xlab="Chr X inbreeding estimate", breaks =100)
hist(female_ibd[,6], main="Female", ylab="Freq", xlab="Chr X inbreeding estimate", breaks=100)

dev.off()
