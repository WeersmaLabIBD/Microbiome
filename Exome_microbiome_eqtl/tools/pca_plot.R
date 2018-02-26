
pdf('LLD_1000G.pdf')

pca <- read.table('LLD_1000G.pca', header = TRUE)
pop <- read.table('kgp.pop', header = TRUE)
lld=read.table("LLD_ID_list.txt",header = F)
lld$FID=0
lld$pop="LLD"
colnames(lld)=c("IID","FID","POP")
id=rbind(pop,lld)
pca=merge(id,pca)

p <- list()
p[[1]] <- ggplot(pca, aes(x=PC1, y=PC2))
p[[2]] <- ggplot(pca, aes(x=PC1, y=PC3))
p[[3]] <- ggplot(pca, aes(x=PC2, y=PC3))
p[[4]] <- ggplot(pca, aes(x=PC1, y=PC4))
p[[5]] <- ggplot(pca, aes(x=PC2, y=PC4))
p[[6]] <- ggplot(pca, aes(x=PC3, y=PC4))
for (i in 1:6) {
  print(p[[i]] + geom_point( aes(color=POP),size=1) +
          theme_bw())
}

dev.off()

pdf('IBD_1000G.pdf')

pca <- read.table('IBD_1000G.pca', header = TRUE)
pop <- read.table('kgp.pop', header = TRUE)
ibd=read.table("IBD_ID_list.txt",header = F)
ibd$FID=0
ibd$pop="IBD"
colnames(ibd)=c("IID","FID","POP")
id=rbind(pop,ibd)
pca=merge(id,pca)

p <- list()
p[[1]] <- ggplot(pca, aes(x=PC1, y=PC2))
p[[2]] <- ggplot(pca, aes(x=PC1, y=PC3))
p[[3]] <- ggplot(pca, aes(x=PC2, y=PC3))
p[[4]] <- ggplot(pca, aes(x=PC1, y=PC4))
p[[5]] <- ggplot(pca, aes(x=PC2, y=PC4))
p[[6]] <- ggplot(pca, aes(x=PC3, y=PC4))
for (i in 1:6) {
  print(p[[i]] + geom_point( aes(color=POP),size=1) +
          theme_bw())
}

dev.off()