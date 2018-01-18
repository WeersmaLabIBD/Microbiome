# 16S mice data follow up analysis

# Creator: Shixian
# Year: 2018

# This script is to analyze the outfiles from qiiem2, including:  
#                                                taxonomy files: level-2.csv,level-3.csv,level-4.csv,level-5.csv,level-6.csv
#                                                  shannon file: alpha-diversity.tsv
#                                                 rarefied outs: rarefied_table.tsv
#                                                 metadata file: Metadata_16S_de2.tsv
#
# Prepare all the files above, then try to run the following commands.


#----------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------taxonomy--------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

# download each level taxa from qiime2 result $taxa-bar-plots.qzv

pylum_data=read.table(file = "level-2.csv",sep = ",",header = T)
class_data=read.table(file = "level-3.csv",sep = ",",header=T)
order_data=read.table(file = "level-4.csv",sep = ",",header = T)
family_data=read.table(file = "level-5.csv",sep = ",",header = T)
genus_data=read.table(file = "level-6.csv",sep = ",",header = T)

# open metadata_de2(already removed 2 sample, 83 samples left)

meta_info=read.table(file = "Metadata_16S_de2.tsv",sep = ";",header = T)

all_taxa=merge(pylum_data[,1:which(colnames(pylum_data)=="BarcodeSequence")-1],class_data[,1:which(colnames(class_data)=="BarcodeSequence")-1],all = T,by="index")
all_taxa=merge(all_taxa,order_data[,1:which(colnames(order_data)=="BarcodeSequence")-1],all = T,by="index")
all_taxa=merge(all_taxa,family_data[,1:which(colnames(family_data)=="BarcodeSequence")-1],all = T,by="index")
all_taxa=merge(all_taxa,genus_data[,1:which(colnames(genus_data)=="BarcodeSequence")-1],all = T,by="index")
colnames(all_taxa)[1]="SampleID"

# remove 2 samples from all_taxa

all=merge(meta_info,all_taxa,all = T,by="SampleID")
all=na.omit(all)
rownames(all)=all$SampleID
all=all[order(all[,which(colnames(all)=="Group")],decreasing = T),]
taxa=all[,-1:-10]

# check non_zero sample frequency --------------------------------------------------------------------------------------------

group_age=all[grep("ACTR0",all$Group),-1:-10]
group_normal=all[grep("YCTR0",all$Group),-1:-10]
group_ppi_four=all[grep("YOME4",all$Group),-1:-10]
group_ppi_eight=all[grep("YOME8",all$Group),-1:-10]
group_ppi_twielve=all[grep("YOME12",all$Group),-1:-10]
group_non_four=all[grep("YCTR4",all$Group),-1:-10]
group_non_eight=all[grep("YCTR8",all$Group),-1:-10]
group_non_twielve=all[grep("YCTR12",all$Group),-1:-10]

groups=c("ACTR0","YCTR0","YOME4","YOME8","YOME12","YCTR4","YCTR8","YCTR12")

aa=c()
for(i in groups){
  
  mm=all[grep(i,all$Group),-1:-10]
  nn=matrix(nrow = 1,ncol = ncol(mm))
  rownames(nn)="non_zero number"
  colnames(nn)=colnames(mm)
  
  for(f in 1:ncol(mm)){
    
    sum=sum(mm[,f]!=0)
    nn[1,f]=sum
    nn=as.data.frame(nn)
  }
  
  y=paste(i,"non_zero",sep = "_")
  aa=append(aa,y)
  
  # give the data frame to string, this function is very very important
  
  assign(y,nn)
  
}

library(ggplot2)
library(gridExtra)

# draw frequency plot for 8 groups

figure=c()
for(i in aa){
  
  mm=as.data.frame(t(get(i)))
  y=paste(i,"figure",sep = "_")
  figure=append(figure,y)
  
  assign(y,
  ggplot(data=mm,aes(`non_zero number`))+
    geom_histogram(color='white',fill='grey80',binwidth = 1)+
    ylab(label = 'taxa frequency')+
    theme_classic()+labs(title = i)+
    theme(plot.title = element_text(size = 10, vjust = 1))+
    scale_x_continuous(breaks=seq(0, 20, 1))
  )
}

# this is only a function to put multiple plots together

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

png("non_zero_frequency.png",width = 1400,height = 1400)

multiplot(ACTR0_non_zero_figure,YCTR0_non_zero_figure,YOME4_non_zero_figure,YOME8_non_zero_figure,YOME12_non_zero_figure,YCTR4_non_zero_figure,YCTR8_non_zero_figure,YCTR12_non_zero_figure,cols=2) 

dev.off()

# filter taxa -----------------------------------------------------------------------------------------------------------
# filter taxa in terms of its presence within each group, filtered means surely exist
# here i choose <20% which is about 2 samples within each group as absence, due to the frequency histogram

groups=c("ACTR0","YCTR0","YOME4","YOME8","YOME12","YCTR4","YCTR8","YCTR12")

cutoff=0.2
for(i in groups){
  
  mm=taxa[grep(i,all$Group),]
  mm = mm[,(colSums(mm > 0) >= cutoff*nrow(mm))]
  mm[mm == 0] = NA
  
  y=paste(i,"filtered.tsv",sep = "_")
  write.table(mm,file = y,sep="\t",quote = F)
  
}

# caculate taxa number ------------------------------------------------------------------------------------------------------
# re-count taxa number after filtering

count_taxa=read.table(file = "ACTR0_filtered.tsv",sep = "\t",header = T)
count_taxa=as.data.frame(t(count_taxa))
count_taxa$taxaID=rownames(count_taxa)
groups=c("YCTR0","YOME4","YOME8","YOME12","YCTR4","YCTR8","YCTR12")

for(i in groups){
  
  y=paste(i,"filtered.tsv",sep = "_")
  mm=read.table(file = y,sep = "\t",header = T)
  mm=as.data.frame(t(mm))
  mm$taxaID=rownames(mm)
  count_taxa=merge(count_taxa,mm,all=T,by="taxaID")
  
}

rownames(count_taxa)=count_taxa[,1]
count_taxa=count_taxa[,-1]
all_filtered=as.data.frame(t(count_taxa))

write.table(count_taxa,file = "total_taxa.tsv",sep = "\t")

count_genus=all_filtered[,grep("g__",colnames(all_filtered),invert = F)]
print(paste("genus number is: ", ncol(count_genus)))

count_family=all_filtered[,grep("f__",colnames(all_filtered),invert = F)]
count_family=count_family[,grep("g__",colnames(count_family),invert = T)]
print(paste("family number is: ", ncol(count_family)))

count_order=all_filtered[,grep("o__",colnames(all_filtered),invert = F)]
count_order=count_order[,grep("f__",colnames(count_order),invert = T)]
print(paste("order number is: ", ncol(count_order)))

count_class=all_filtered[,grep("c__",colnames(all_filtered),invert = F)]
count_class=count_class[,grep("o__",colnames(count_class),invert = T)]
print(paste("class number is: ", ncol(count_class)))

count_pylum=all_filtered[,grep("p__",colnames(all_filtered),invert = F)]
count_pylum=count_pylum[,grep("c__",colnames(count_pylum),invert = T)]
print(paste("phylum number is: ", ncol(count_pylum)))

# calculate relative aubndance -----------------------------------------------------------------------------------------------

groups=c("ACTR0","YCTR0","YOME4","YOME8","YOME12","YCTR4","YCTR8","YCTR12")

# genus relative for each group

for(i in groups){
  ii=paste(i,"filtered.tsv",sep = "_")
  aa=read.table(file = ii,sep = "\t",header = T)
  genus=aa[,grep("g__",colnames(aa))]
  genus[is.na(genus)]=0
  bb=matrix(nrow = nrow(aa)+1,ncol = ncol(genus))
  colnames(bb)=colnames(genus)
  
  rownames(bb)=c(rownames(aa),"average")
  
  for(m in 1:nrow(genus)){
    
    sum=sum(genus[m,])
    
    for(f in 1:ncol(genus)){
      
      bb[m,f]=genus[m,f]/sum
      
    }
    
  }
  
  for(q in 1:ncol(bb)){
    
    colmean=mean(bb[1:nrow(bb)-1,q])
    bb[nrow(aa)+1,q]=colmean
    
  }
  
  y=paste(i,"genus_relative.tsv",sep = "_")
  write.table(bb,y,sep = "\t")
  
}

# family relative for each group

for(i in groups){
  
  ii=paste(i,"filtered.tsv",sep = "_")
  aa=read.table(file = ii,sep = "\t",header = T)
  family=aa[,grep("f__",colnames(aa))]
  family=family[,grep("g__",colnames(family),invert=T)]
  family[is.na(family)]=0
  bb=matrix(nrow = nrow(aa)+1,ncol = ncol(family))
  colnames(bb)=colnames(family)
  
  rownames(bb)=c(rownames(aa),"average")
  
  for(m in 1:nrow(family)){
    
    sum=sum(family[m,])
    
    for(f in 1:ncol(family)){
      
      bb[m,f]=family[m,f]/sum
      
    }
    
  }
  
  for(q in 1:ncol(bb)){
    
    colmean=mean(bb[1:nrow(bb)-1,q])
    bb[nrow(aa)+1,q]=colmean
    
  }
  
  y=paste(i,"family_relative.tsv",sep = "_")
  write.table(bb,y,sep = "\t")
  
}

# order relative for each group

for(i in groups){
  
  ii=paste(i,"filtered.tsv",sep = "_")
  aa=read.table(file = ii,sep = "\t",header = T)
  order=aa[,grep("o__",colnames(aa))]
  order=order[,grep("f__",colnames(order),invert=T)]
  order[is.na(order)]=0
  bb=matrix(nrow = nrow(aa)+1,ncol = ncol(order))
  colnames(bb)=colnames(order)
  
  rownames(bb)=c(rownames(aa),"average")
  
  for(m in 1:nrow(order)){
    
    sum=sum(order[m,])
    
    for(f in 1:ncol(order)){
      
      bb[m,f]=order[m,f]/sum
      
    }
    
  }
  
  for(q in 1:ncol(bb)){
    
    colmean=mean(bb[1:nrow(bb)-1,q])
    bb[nrow(aa)+1,q]=colmean
    
  }
  
  y=paste(i,"order_relative.tsv",sep = "_")
  write.table(bb,y,sep = "\t")
  
}

# class relative for each group

for(i in groups){
  
  ii=paste(i,"filtered.tsv",sep = "_")
  aa=read.table(file = ii,sep = "\t",header = T)
  class=aa[,grep("c__",colnames(aa))]
  class=class[,grep("o__",colnames(class),invert=T)]
  class[is.na(class)]=0
  bb=matrix(nrow = nrow(aa)+1,ncol = ncol(class))
  colnames(bb)=colnames(class)
  
  rownames(bb)=c(rownames(aa),"average")
  
  for(m in 1:nrow(class)){
    
    sum=sum(class[m,])
    
    for(f in 1:ncol(class)){
      
      bb[m,f]=class[m,f]/sum
      
    }
    
  }
  
  for(q in 1:ncol(bb)){
    
    colmean=mean(bb[1:nrow(bb)-1,q])
    bb[nrow(aa)+1,q]=colmean
    
  }
  
  y=paste(i,"class_relative.tsv",sep = "_")
  write.table(bb,y,sep = "\t")
  
}

# phylum relative for each group

for(i in groups){
  
  ii=paste(i,"filtered.tsv",sep = "_")
  aa=read.table(file = ii,sep = "\t",header = T)
  phylum=aa[,grep("p__",colnames(aa))]
  phylum=phylum[,grep("c__",colnames(phylum),invert=T)]
  phylum[is.na(phylum)]=0
  bb=matrix(nrow = nrow(aa)+1,ncol = ncol(phylum))
  colnames(bb)=colnames(phylum)
  
  rownames(bb)=c(rownames(aa),"average")
  
  for(m in 1:nrow(phylum)){
    
    sum=sum(phylum[m,])
    
    for(f in 1:ncol(phylum)){
      
      bb[m,f]=phylum[m,f]/sum
      
    }
    
  }
  
  for(q in 1:ncol(bb)){
    
    colmean=mean(bb[1:nrow(bb)-1,q])
    bb[nrow(aa)+1,q]=colmean
    
  }
  
  y=paste(i,"phylum_relative.tsv",sep = "_")
  write.table(bb,y,sep = "\t")
  
}

# total relative abundace

groups=c("ACTR0","YCTR0","YOME4","YOME8","YOME12","YCTR4","YCTR8","YCTR12")
taxa=c("phylum","class","order","family","genus")

for(i in groups){
  
  name=read.table("Metadata_16S_de2.tsv",sep = ";",header = T)
  name=name[grep(i,name$Group),]
  y=data.frame(name$SampleID)
  colnames(y)="SampleID"
  for(n in taxa){
    x=paste(i,n,"relative.tsv",sep = "_")
    m=read.table(file = x,header = T,sep = "\t")
    m$SampleID=rownames(m)
    y=merge(y,m,all=T,by="SampleID")
    
  }
  
  z=paste(i,"taxa",sep = "_")
  assign(z,y)
  write.table(get(z),file = z, sep = "\t")
  
}

aa=read.table("total_taxa.tsv",sep = "\t",header = T)
bb=data.frame(rownames(aa))
colnames(bb)="taxaID"
for(i in groups){
  
  y=paste(i,"taxa",sep = "_")
  m=read.table(file = y)
  m=m[-1,]
  rownames(m)=m$SampleID
  m=m[,-1]
  m=as.data.frame(t(m))
  m$taxaID=rownames(m)
  bb=merge(bb,m,all=T,by="taxaID")
}
bb[is.na(bb)]=0

write.table(bb,file = "total_relative_abundance.tsv",sep = "\t")



# devide data to binary and quantitative classification-----------------------------------------------------------------------
# select continus taxa data

groups=c("YCTR0","YOME4","YOME8","YOME12","YCTR4","YCTR8","YCTR12")

continus_taxa=read.table(file = "ACTR0_filtered.tsv",sep = "\t",header = T)
continus_taxa=as.data.frame(t(continus_taxa))
continus_taxa$taxaID=rownames(continus_taxa)
for(i in groups){
  
  ii=paste(i,"filtered.tsv",sep = "_")
  aa=read.table(file = ii,sep = "\t",header = T)
  bb=as.data.frame(t(aa))
  bb$taxaID=rownames(bb)
  continus_taxa=merge(continus_taxa,bb,by="taxaID")
  
}

# select binary taxa data
# between age and normal

filter_age=read.table(file = "ACTR0_filtered.tsv",sep = "\t",header = T)
filter_age=as.data.frame(t(filter_age))
filter_age$taxaID=rownames(filter_age)
filter_normal=read.table(file = "YCTR0_filtered.tsv",sep = "\t",header =T)
filter_normal=as.data.frame(t(filter_normal))
filter_normal$taxaID=rownames(filter_normal)

binary_taxa_age=merge(filter_normal,filter_age,all=T,by="taxaID")

for(i in continus_taxa$taxaID){
  
  n=grep(i, binary_taxa_age)
  binary_taxa_age=binary_taxa_age[-n,]
  
}

rownames(binary_taxa_age)=binary_taxa_age[,1]
binary_taxa_age=binary_taxa_age[,-1]
binary_taxa_age=as.data.frame(t(binary_taxa_age))

# select binary taxa data
# between ppi and non_ppi

groups=c("YOME8","YOME12","YCTR4","YCTR8","YCTR12")

filter_YOME4=read.table(file = "YOME4_filtered.tsv",sep = "\t",header =T)
filter_YOME4=as.data.frame(t(filter_YOME4))
filter_YOME4$taxaID=rownames(filter_YOME4)

binary_taxa_ppi=filter_YOME4
for(i in groups){
  
  ii=paste(i,"filtered.tsv",sep = "_")
  aa=read.table(file = ii,sep = "\t",header = T)
  bb=as.data.frame(t(aa))
  bb$taxaID=rownames(bb)
  binary_taxa_ppi=merge(binary_taxa_ppi,bb,by="taxaID")
  
}

for(i in continus_taxa$taxaID){
  
  n=grep(i, binary_taxa_ppi)
  binary_taxa_ppi=binary_taxa_ppi[-n,]
  
}
rownames(binary_taxa_ppi)=binary_taxa_ppi[,1]
binary_taxa_ppi=binary_taxa_ppi[,-1]
binary_taxa_ppi=as.data.frame(t(binary_taxa_ppi))


# composition------------------------------------------------------------------------------------------------------------------

# only mark taxa is present more than 1% of samples within each group

# age vs normal

groups=c("phylum","class","order","family","genus")

for(i in groups){
  
  x=paste("ACTR0",i,"relative.tsv",sep = "_")
  y=paste("YCTR0",i,"relative.tsv",sep = "_")
  p=paste("age_normal", i,".png",sep = "_")
  
  png(p,width=700,height=700)
  
  age=as.data.frame(t(read.table(file = x)))
  normal=as.data.frame(t(read.table(file = y)))
  
  tmp_age=data.frame(age[age[,ncol(age)]>0.01,][,ncol(age)],row.names = rownames(age[age[,ncol(age)]>0.01,]))
  colnames(tmp_age)="average"
  tmp_age$group="age"
  tmp_age$taxaID=rownames(tmp_age)
  tmp_normal=data.frame(normal[normal[,ncol(normal)]>0.01,][,ncol(normal)],row.names = rownames(normal[normal[,ncol(normal)]>0.01,]))
  colnames(tmp_normal)="average"
  tmp_normal$group="normal"
  tmp_normal$taxaID=rownames(tmp_normal)
  
  minority_age=c(1-sum(tmp_age$average),"age","less than 0.01")
  minority_normal=c(1-sum(tmp_normal$average),"normal","less than 0.01")
  
  tmp_age=rbind(tmp_age,minority_age)
  tmp_normal=rbind(tmp_normal,minority_normal)
  tmp=rbind(tmp_normal,tmp_age)
  
  tmp$percentage=as.numeric(tmp$average)

  Group=factor(tmp$group,c("normal","age"))

  print(ggplot(data = tmp) + 
    geom_bar( aes(Group,percentage,fill=taxaID),stat= 'identity', position = 'stack')+
    theme_classic())
  
  dev.off()
  
}


# ppi vs non_ppi

groups=c("phylum","class","order","family","genus")

for(i in groups){
  
  x=paste("YOME4",i,"relative.tsv",sep = "_")
  y=paste("YOME8",i,"relative.tsv",sep = "_")
  w=paste("YOME12",i,"relative.tsv",sep = "_")
  e=paste("YCTR4",i,"relative.tsv",sep = "_")
  r=paste("YCTR8",i,"relative.tsv",sep = "_")
  t=paste("YCTR12",i,"relative.tsv",sep = "_")
  p=paste("ppi_ppi_non", i,".png",sep = "_")
  
  png(p,width=700,height=700)
  
  ppi_4=as.data.frame(t(read.table(file = x)))
  ppi_8=as.data.frame(t(read.table(file = y)))
  ppi_12=as.data.frame(t(read.table(file = w)))
  non_4=as.data.frame(t(read.table(file = e)))
  non_8=as.data.frame(t(read.table(file = r)))
  non_12=as.data.frame(t(read.table(file = t)))
  
  mm=c("ppi_4","ppi_8","ppi_12","non_4","non_8","non_12")
  
  for(k in mm){
    
    nn=paste("tmp",k,sep = "_")
    qq=get(k)
    pp=data.frame(qq[qq[,ncol(qq)]>0.01,][,ncol(qq)],row.names = rownames(qq[qq[,ncol(qq)]>0.01,]))
    assign(nn,pp)
    ww=get(nn)
    colnames(ww)="average"
    ww$group=k
    ww$taxaID=rownames(ww)
    hh=c(1-sum(ww$average),k,"less than 0.01")
    ww=rbind(ww,hh)
    assign(nn,ww)
  }
  
  tmp=rbind(tmp_ppi_4,tmp_ppi_8,tmp_ppi_12,tmp_non_4,tmp_non_8,tmp_non_12)
  
  tmp$percentage=as.numeric(tmp$average)
  
  Group=factor(tmp$group,c("ppi_4","ppi_8","ppi_12","non_4","non_8","non_12"))
  
  print(ggplot(data = tmp) + 
          geom_bar( aes(Group,percentage,fill=taxaID),stat= 'identity', position = 'stack')+
          theme_classic())
  
  dev.off()
  
}

# ------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------shannon diversity---------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------

library(ggplot2)
meta_data=read.table(file = "Metadata_16S_de2.tsv", sep = ";", header = T)
colnames(meta_data)[1]="ID"

shannon_input=read.table(file = "alpha-diversity.tsv",sep = "\t", header = T)
colnames(shannon_input)[1]="ID"

shannon=merge(shannon_input,meta_data,all = T,by="ID")
shannon=na.omit(shannon)

# total view shannon ------------------------------------------------------------------------------------------------------------
# box plot

group=factor(shannon$Group,c("ACTR0","YCTR0","YCTR4","YCTR8","YCTR12","YOME4","YOME8","YOME12"))

png("total_shannon_boxplot.png",width=700,height=700)

ggplot (shannon, aes(group,shannon)) + 
  geom_boxplot(aes(fill = PPI)) + 
  scale_fill_manual(values =c("green","pink"))+
  theme_classic() +
  xlab(label = "Groups") + ylab(label = "shannon") +
  theme(legend.text = element_text(size = 8, face = "bold"),legend.position = c(0.08,0.9))+
  scale_x_discrete(breaks= c("ACTR0","YCTR0","YCTR4","YCTR8","YCTR12","YOME4","YOME8","YOME12"),
  labels = c("age","normal","non_ppi \n 4 weeks","non_ppi \n 8 weeks","non_ppi \n 12 weeks","ppi \n 4 weeks","ppi \n 8 weeks","ppi \n 12 weeks"))+
  theme(axis.text = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  theme(axis.title = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text( size=14)) +
  theme(legend.text = element_text( size=10))+
  guides(fill=FALSE)

dev.off()

# line dot plot

png("total_shannon_lineplot.png",width=700,height=700)

shannon$factor_1=paste(shannon$Mice,shannon$PPI,sep = "_")
shannon$factor_2=c(1:83)
shannon$factor_2[grep("_Yes", shannon$factor_1)]=shannon$factor_1[grep("_Yes", shannon$factor_1)]
shannon$factor_3=c(1:83)
shannon$factor_3[grep("_No", shannon$factor_1)]=shannon$factor_1[grep("_No", shannon$factor_1)]
shannon=as.data.frame(shannon)

# use factor re-order $Timepoint
shannon$timepoint=shannon$Timepoint
shannon$timepoint[1:7]="age"
group=factor(shannon$timepoint,c("age","0","4","8","12"))

ggplot (shannon, aes(group,shannon)) + 
  geom_point(aes(colour = PPI),size=4) + 
  theme_classic() +
  xlab(label = "") + ylab(label = "shannon") +
  geom_line(aes(group=factor_2),linetype="dashed",alpha=.6,colour = "red")+
  geom_line(aes(group=factor_3),linetype="dotdash",alpha=1,colour = "green")+
  theme(legend.text = element_text(size = 8, face = "bold"),legend.position = c(0.08,0.9))+
  scale_color_hue(direction = 1, h.start=150)+
  scale_x_discrete(breaks= c("age","0","4","8","12"),labels = c("age","normal","4 weeks","8 weeks","12 weeks"))+
  theme(axis.text = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  theme(axis.title = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text( size=14)) +
  theme(legend.text = element_text( size=10))

dev.off()

# age vs normal------------------------------------------------------------------------------------------------------------------------

# wilcoxon test for shannon

age=shannon[grep("ACTR0",shannon$Group),]
normal=shannon[grep("YCTR0",shannon$Group),]
  
tmp=rbind(age,normal)
wilcox.test(age$shannon,normal$shannon,exact = T, correct = FALSE)

# re-draw plot with p value

png("age_normal_shannon.png",width=700,height=700)

group=factor(tmp$Group,c("YCTR0","ACTR0"))

ggplot (tmp, aes(group,shannon)) + 
  geom_boxplot(aes(fill = group)) + 
  scale_fill_manual(values =c("green","pink"))+
  theme_classic() +
  xlab(label = "Groups") + ylab(label = "shannon") +
  theme(legend.text = element_text(size = 8, face = "bold"),legend.position = c(0.08,0.9))+
  scale_x_discrete(breaks= c("YCTR0","ACTR0"),
                   labels = c("normal","age"))+
  annotate("text", x = 1.5, y = 2, label = "p value=0.003",size=4)+
  theme(axis.text = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  theme(axis.title = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text( size=14)) +
  theme(legend.text = element_text( size=10))+
  guides(fill=FALSE)

dev.off()

# ppi vs non_ppi-----------------------------------------------------------------------------------------------------------------
options(scipen=4)

# 12 weeks
ppi=shannon[grep("YOME12",shannon$Group),]
ppi_non=shannon[grep("YCTR12",shannon$Group),]
a=wilcox.test(ppi_non$shannon,ppi$shannon,exact = T, correct = FALSE)
print(paste("12 week, ppi vs non_ppi : p_value is", a$p.value))

# 8 weeks
ppi=shannon[grep("YOME8",shannon$Group),]
ppi_non=shannon[grep("YCTR8",shannon$Group),]
a=wilcox.test(ppi_non$shannon,ppi$shannon,exact = T, correct = FALSE)
print(paste("8 week, ppi vs non_ppi : p_value is", a$p.value))

# 4 weeks
ppi=shannon[grep("YOME4",shannon$Group),]
ppi_non=shannon[grep("YCTR4",shannon$Group),]
a=wilcox.test(ppi_non$shannon,ppi$shannon,exact = T, correct = FALSE)
print(paste("8 week, ppi vs non_ppi : p_value is", a$p.value))

# re draw box plot with p_value

png("ppi_shannon.png",width=700,height=700)

tmp=shannon[-1:-27,]
group=factor(tmp$Timepoint,c("4","8","12"))
tmp$hline=2.8

ggplot (tmp, aes(group,shannon)) + 
  geom_boxplot(aes(fill = PPI)) + 
  scale_fill_manual(values =c("green","pink"))+
  theme_classic() +
  xlab(label = "Groups") + ylab(label = "shannon") +
  theme(legend.text = element_text(size = 8, face = "bold"),legend.position = c(0.08,0.9))+
  scale_x_discrete(breaks= c("4","8","12"),labels = c("4 weeks","8 weeks","12 weeks"))+
  geom_errorbar(aes(y=hline, ymax=hline, ymin=hline), linetype="dashed")+
  annotate("text", x = 1, y = 3, label = "NS",size=4)+
  annotate("text", x = 2, y = 3, label = "p value=0.00004",size=4)+
  annotate("text", x = 3, y = 3, label = "NS",size=4)+
  theme(axis.text = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  theme(axis.title = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  theme(legend.title = element_text( size=14)) +
  theme(legend.text = element_text( size=10))

dev.off()

# --------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------beta diversity---------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------

library(vegan)
library(ggplot2)
library(gganimate)
library(animation)

# rarefied feature table to calculate PCoA

rare_data=read.table(file ="rarefied_table.tsv",header = T,sep = "\t" )
meta_data=read.table(file = "Metadata_16S_de2.tsv", sep = ";", header = T)

rownames(meta_data)=meta_data$SampleID
rownames(rare_data)=rare_data$"OTU.ID"

rare_data=rare_data[,-1]
rare_data=t(rare_data)

rare_data=merge(meta_data,rare_data,all = T, by="row.names")
rownames(rare_data)=rare_data$SampleID

# total view ----------------------------------------------------------------------------------------------------------------

pure_data=rare_data[,-1:-11]

beta_diversity=vegdist(pure_data,method = "bray")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=5))
table=merge(pca_analysis,meta_data,all = T,by="row.names")

png("total_beta.png",width=700,height=700)

ggplot (table,aes(V1,V2)) +
  geom_point(aes(color=Group,shape=PPI),size=3) + theme_bw() + xlab(label = "PC1") + ylab(label = "PC2")+
  theme_bw()

dev.off()

# age vs normal ----------------------------------------------------------------------------------------------------------------

age=rare_data[grep("ACTR0",rare_data$Group),]
normal=rare_data[grep("YCTR0",rare_data$Group),]
tmp=rbind(normal,age)
tmp_data=tmp[,-1:-11]

beta_diversity=vegdist(tmp_data,method = "bray")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=5))
table=merge(pca_analysis,meta_data,all = T,by="row.names")
table=na.omit(table)

png("aging_normal_beta.png",width=700,height=700)

ggplot (table,aes(V1,V2)) +
  geom_point(aes(color=Group),size=3) + theme_bw() + xlab(label = "PC1") + ylab(label = "PC2")+
  scale_colour_discrete(breaks = c("ACTR0","YCTR0"), labels = c("aging","normal"))+
  theme_bw()

dev.off()

# ppi vs non_ppi --------------------------------------------------------------------------------------------------------------

tmp_data=rare_data[grep("ACTR0",rare_data$Group,invert=T),]
tmp_data=tmp_data[grep("YCTR0",tmp_data$Group,invert=T),]

beta_diversity=vegdist(tmp_data[,-1:-11],method = "bray")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=5))
table=merge(pca_analysis,meta_data,all = T,by="row.names")
table=na.omit(table)

Time_point=factor(table$Timepoint,c("four","eight","twielve"))

table=as.data.frame(table)

p <- ggplot(table, aes(V1, V2, color = Mice, size=3,frame = Timepoint,shape=PPI)) +
  geom_point() +theme_bw() + xlab(label = "PC1") + ylab(label = "PC2")+
  facet_wrap(~PPI, scales = "fixed")+
  theme(legend.title = element_text( size=14)) +
  theme(legend.text = element_text( size=12))+
  ggtitle(" Week: ")+
  theme(plot.title = element_text(size = 20, face = "bold",vjust = 0.5, hjust = 0.5))+
  theme(axis.text = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  theme(axis.title = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  guides(size=FALSE)
gganimate(p,interval = 1.5,"ppi_seperate.gif")

p <- ggplot(table, aes(V1, V2, color = PPI, size=3,frame = Timepoint)) +
  geom_point() +theme_bw() + xlab(label = "PC1") + ylab(label = "PC2")+
  facet_wrap(~PPI, scales = "fixed")+
  theme(legend.title = element_text( size=18,face = "bold")) +
  theme(legend.text = element_text( size=14))+
  ggtitle(" Week: ")+
  theme(plot.title = element_text(size = 20, face = "bold",vjust = 0.5, hjust = 0.5))+
  theme(axis.text = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  theme(axis.title = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  scale_color_manual(values =c("green","red"))+
  guides(size=FALSE)
gganimate(p,interval = 1.5,"ppi_seperate_b_h.gif")

p <- ggplot(table, aes(V1, V2, color = Mice, size=3,frame = Timepoint,shape=PPI)) +
  geom_point() +theme_bw() + xlab(label = "PC1") + ylab(label = "PC2")+
  theme(legend.title = element_text( size=14)) +
  theme(legend.text = element_text( size=12))+
  ggtitle(" Week: ")+
  theme(plot.title = element_text(size = 20, face = "bold",vjust = 0.5, hjust = 0.5))+
  theme(axis.text = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  theme(axis.title = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  guides(size=FALSE)
gganimate(p,interval = 1.5,"ppi_all_colorful.gif")

p <- ggplot(table, aes(V1, V2, color = PPI, size=3,frame = Timepoint)) +
  geom_point() +theme_bw() + xlab(label = "PC1") + ylab(label = "PC2")+
  theme(legend.title = element_text( size=18,face = "bold")) +
  theme(legend.text = element_text( size=14))+
  ggtitle(" Week: ")+
  theme(plot.title = element_text(size = 20, face = "bold",vjust = 0.5, hjust = 0.5))+
  theme(axis.text = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  theme(axis.title = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  scale_color_manual(values =c("green","red"))+
  guides(size=FALSE)
gganimate(p,interval = 1.5,"ppi_all_b_h.gif")

Time_point=factor(table$Timepoint,c("4","8","12"))

colnames(table)[2:3]=c("xis","yis")
table_1=reshape(table, v.names="yis", idvar = "Mice", timevar = "Timepoint", direction = "wide")
table_2=reshape(table, v.names="xis", idvar = "Mice", timevar = "Timepoint", direction = "wide")
n=ncol(table_1)
table_segment=merge(table_2[,c(9,11,12,13,15,n,n-1,n-2)],table_1[,c(9,11,12,13,15,n,n-1,n-2)],all=T,by="Mice")

png("ppi_beta.png",width=700,height=700)

ggplot (table,aes(xis,yis)) +
  geom_point(aes(color=Time_point,shape=PPI),size=4) + theme_bw() + xlab(label = "PC1") + ylab(label = "PC2")+
  scale_colour_brewer(palette = "Reds")+
  theme_bw()+
  geom_segment(data = table_segment,aes(x = xis.4, y = yis.4, xend = xis.8, yend = yis.8),alpha=0.5,linetype="dotted",       
               arrow = arrow(length=unit(0.25,"cm"),angle = 15, type = "closed"))+
  geom_segment(data = table_segment,aes(x = xis.8, y = yis.8, xend = xis.12, yend = yis.12),alpha=0.5,linetype="dotted",         
               arrow = arrow(length=unit(0.25,"cm"),angle = 15, type = "closed"))+
  guides(linetype=FALSE)

dev.off()









