for(i in sp[c(3,6)]){
  data <- gea_p[which(gea_p$sp == i),"V2"]
  i_new <- gsub(" ", "", i)
  write.table(data,file = paste0("E:\\Pinus SNP paper\\gea_sample\\",i_new,".txt"),col.names = F,row.names = F,quote = F,sep = "\t")
}
for(i in c(1:6)){
  sp[i] <- gsub("\t", "", sp[i])
}
data <- data.frame(name = c(rep("21Pk",4),rep("27Pm",4),rep("19Pd",4),rep("36Py",4),rep("3Pt",4),rep("3Pyv",4)),stat_name = rep(c("He","HardyHo","Ho","Pi"),6),mean = rep(NA,24),sd = rep(NA,24))
for(i in sp){
  He <- read.table(dir_list[grep(paste0(i,"_He.txt"),dir_list)],header = F,stringsAsFactors = F)
  He_value_m <- mean(as.numeric(He[,3]),na.rm = T)
  He_vale_d <- sd(as.numeric(He[,3]),na.rm = T)
  hardyHo <- read.table(dir_list[grep(paste0(i,"_hardyHo.txt"),dir_list)],header = F,stringsAsFactors = F)
  Ha_value_m <- mean(as.numeric(hardyHo[,3]),na.rm = T)
  Ha_vale_d <- sd(as.numeric(hardyHo[,3]),na.rm = T)
  Ho <- read.table(dir_list[grep(paste0(i,"_Ho.txt"),dir_list)],header = F,stringsAsFactors = F) 
  Ho_value_m <- mean(as.numeric(Ho[,2]),na.rm = T)
  Ho_vale_d <- sd(as.numeric(Ho[,2]),na.rm = T)
  pi <- read.table(dir_list[grep(paste0(i,"_pi.sites"),dir_list)],header = T,stringsAsFactors = F)
  pi_value_m <- mean(as.numeric(pi[,3]),na.rm = T)
  pi_vale_d <- sd(as.numeric(pi[,3]),na.rm = T)
  data[which(data[,1] == i),"mean"] <- c(He_value_m,Ha_value_m,Ho_value_m,pi_value_m)
  data[which(data[,1] == i),"sd"] <- c(He_vale_d,Ha_vale_d,Ho_vale_d,pi_vale_d)
}
