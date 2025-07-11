pdf("E:\\Pinus SNP paper\\gea_sample\\imiss.pdf",width = 12,height = 8)
par(mfrow=c(2,3))
for(p in sp){
  d <- imiss[which(imiss$sp == p),5]
  d_new <- 1-d
  hist(d_new,main = NULL,xlab = p,ylab = "Number of samples",col = "#70a1d7",border = "black", cex.lab = 2,cex.axis = 1.5)
  abline(v = mean(d_new),col = "#f47c7c",lty = 2,lwd = 2)
}
dev.off()
lmiss <- read.table("E:\\Pinus SNP paper\\gea_sample\\gea.lmiss",header = T,stringsAsFactors = F)
pdf("E:\\Pinus SNP paper\\gea_sample\\lmiss.pdf",width = 4,height = 4)
hist((1-lmiss$F_MISS),main = NULL,xlab = p,ylab = "Number of samples",col = "#70a1d7",border = "black", cex.lab = 1.5,cex.axis = 1.5)
abline(v = mean((1-lmiss$F_MISS)),col = "#f47c7c",lty = 2,lwd = 2)
dev.off()
imiss_109 <- openxlsx::read.xlsx("E:\\Pinus SNP paper\\109sampleFinal_专利写作中只用到了109个个体的信息.xlsx",colNames = T,sheet = 3)
pdf("E:\\Pinus SNP paper\\109_sample_poly\\imiss.pdf",width = 12,height = 8)
par(mfrow=c(2,3))
for(p in sp){
  d <- imiss_109[which(imiss_109$species == p),5]
  d_new <- 1-d
  hist(d_new,main = NULL,xlab = p,ylab = "Number of samples",col = "#70a1d7",border = "black", cex.lab = 2,cex.axis = 1.5)
  abline(v = mean(d_new),col = "#f47c7c",lty = 2,lwd = 2)
}
dev.off()
lmiss_109 <- read.table("E:\\Pinus SNP paper\\109_sample_poly\\site.lmiss",header = T,stringsAsFactors = F)
pdf("E:\\Pinus SNP paper\\109_sample_poly\\lmiss.pdf",width = 4,height = 4)
hist((1-lmiss_109$F_MISS),main = NULL,xlab = p,ylab = "Number of samples",col = "#70a1d7",border = "black", cex.lab = 1.5,cex.axis = 1.5)
abline(v = mean((1-lmiss_109$F_MISS)),col = "#f47c7c",lty = 2,lwd = 2)
dev.off()
