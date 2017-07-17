library(foreign)

toweka <- cbind(wm01_02l[[2]],bighlpcrps[[1]][[2]][[2]])
toweka <- rbind(toweka,cbind(bighlpopgr[[1]][[2]],bighlpcrps[[1]][[2]][[3]]))
toweka <- rbind(toweka,cbind(bighlpopgr[[2]][[2]],bighlpcrps[[1]][[2]][[5]]))
colnames(toweka) <- c(paste("cus",seq(1,20),sep=""),"crps","wv45")
rownames(toweka) <- c(paste("rnd",seq(1,160),sep=""),paste("opt",seq(1,2),"_",seq(1,16),sep=""))

write.arff(toweka,'biglp2weka001.arff',eol='\n',relation=('grp'))
