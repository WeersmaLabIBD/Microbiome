mIBD=melt(IBD, id.vars = "SID")
subIBD=subset(mIBD, mIBD$value=="User")
subIBD$value=NULL
data3 <- dcast(subIBD, ...~variable)
rownames(data3)=data3$SID
data3$SID=NULL
data3$united <- apply(data3, 1, function(x) paste(x[!is.na(x) & x != "No"], collapse = "+ "))
#counts_IBD=data3
counts_LLD=data3
counts_MIBS=data3
