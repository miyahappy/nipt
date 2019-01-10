#####for NIPT data gc correction,LOESS regression method,by huliang,2018.3.29#####
#set#
library("ggplot2")
library("dplyr")
rm(list = ls())
file.remove(paste("after_gcc/",list.files("/after_gcc/"),sep=""))
file.remove(paste("stat/",list.files("/summary/"),sep=""))
file.remove(paste("plot/",list.files("/plot/"),sep=""))
path <- "/home/userroot/nipt/gcc"
setwd(path)
set_bed <- "gc_bin_set/gc_bin.bed"
file_names <- list.files("before_gcc")
plot_names <- paste(substring(file_names,1,18),".tiff",sep="")
csv_names <-  paste(substring(file_names,1,18),".csv",sep="")
gcc_counts <- paste(substring(file_names,1,18),".counts",sep="")
plot_dir <- paste("plot/",plot_names,sep="")
csv1_dir <- paste("summary/",csv_names,sep="")
csv2_dir <- paste("stat/",csv_names,sep="")
gcc_dir <- paste("after_gcc/","gcc_",gcc_counts,sep="")
dir <- paste("before_gcc/",file_names,sep="")
n <- length(dir)
#set LOESS regression parameters#
los_span <- 0.75
los_degree <- 1
#bed#
gc <- read.csv(set_bed,header = F,sep = "\t")[-1, ]
a <- as.character(gc$V1)
b <- as.character(gc$V2)
c <- as.character(gc$V3)
bin_name <- paste(a,b,c,sep = "-")
gc <- data.frame(bin_name,as.numeric(as.character(gc$V4)))
##
for (i in 1:n){
    counts <- read.csv(dir[i],header = F,sep = "\t")[ ,-(5:7)]
    a <- as.character(counts$V1)
    b <- as.character(counts$V2)
    cc <- as.character(counts$V3)
    bin_name <- paste(a,b,cc,sep = "-")
    counts <- data.frame(bin_name,as.numeric(as.character(counts$V4)))
    data <- merge(gc,counts,by = "bin_name")
    data1 <- cbind(data,z = "before_gcc")
    data1 <- filter(data1,(data1$as.numeric.as.character.counts.V4.. > 100) & (data1$as.numeric.as.character.counts.V4.. < 5000) ) 
    names(data1) <- c("bin_name","x","y","z")
    x <- data1$x
    y <- data1$y
    los <- loess(y~x,span = los_span,degree = los_degree)
    p <- predict(los,x)
    e <- mean(y)
    #y2 <- y-p+e
    y2 <- y*e/p
    data2 <- data.frame(data1$bin_name,x,y2,z="after_gcc")
    names(data2) <- c("bin_name","x","y","z")
    data3 <- rbind(data1,data2)
    ggplot(data3,aes(x=x,y=y,colour = z))+geom_point(size = 0.5,shape = 18)+xlim(0.25,0.65)+ylim(0,4000)+geom_smooth(method = lm)      
    ggsave(plot_dir[i], dpi=300)
    total_y <-sum(y)
    total_y2 <- sum(y2)
    gc_change <- (total_y2-total_y)/total_y
    gc_pc_before <- sum(y*35*x)/(sum(y)*35)
    gc_pc_after <- sum(y2*35*x)/(sum(y2)*35)
    sample_summary <- data.frame(total_y,total_y2,gc_change,gc_pc_before,gc_pc_after)
    names(sample_summary) <- c("before_reads","after_reads","reads_change_rate","before_GC_percent","after_GC_percent")
    write.csv(data2,csv1_dir[i],row.names = F)
    #write.table(data2[ ,1:3],gcc_dir[i])
    chr <- chartr (" ","-",c(formatC(1:22,flag="-",width=2),"X ","Y ","M "))
    chr <- paste("chr",chr,sep = "")
    chr1_pc <- c()
    chr2_pc <- c()
    chr1_gc_pc <- c()
    chr2_gc_pc <- c()
    for (j in 1:25){
        data1_chr <- filter(data1,substring(bin_name,1,5) == chr[j])
        chr1_pc[j] <- sum(data1_chr$y)/total_y
        chr1_gc_pc[j] <- sum(data1_chr$y*35*data1_chr$x)/(sum(data1_chr$y)*35)
        data2_chr <- filter(data2,substring(bin_name,1,5) == chr[j])
        chr2_pc[j] <- sum(data2_chr$y)/total_y
        chr2_gc_pc[j] <- sum(data2_chr$y*35*data2_chr$x)/(sum(data2_chr$y)*35)
    }
        gc_chr <- data.frame(chr,chr1_pc,chr2_pc,(as.numeric(chr2_pc)-as.numeric(chr1_pc))/as.numeric(chr1_pc),
                      chr1_gc_pc,chr2_gc_pc,(as.numeric(chr2_gc_pc)-as.numeric(chr1_gc_pc))/as.numeric(chr1_gc_pc))
        names(gc_chr) <- c("chromosome","chr%_before_gcc","chr%_after_gcc","chr%_change_rate","gc_percent_before_gcc","gc_percent_after_gcc","gc_percent_change_rate")
    #print(csv2_dir[i])    
    write.csv(gc_chr,csv2_dir[i],row.names = F)
}
