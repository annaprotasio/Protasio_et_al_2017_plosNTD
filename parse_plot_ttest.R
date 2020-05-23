# This script processes raw data from various qPCR files. 
# Collectivelly, it provides the statistical/mathematical steps that were taken to produced the qPCR figues in the paper.

# load libraries (make sure you have these installed in your R session)
library(gdata)
library(plotrix)
library(ggplot2)
library(metafolio)


# female worms time course ------------------------------------------------------

count_files = list.files("raw_data/",full.names=TRUE)
count_files_F = count_files[grep("F",count_files)] #keep only xls files
data_F = lapply(count_files_F, read.xls, sheet=2, skip=6, nrow=79, header=T)

data = data_F #this loop gets the mean of the delta-Ct for the biological replicates, and calculates the SEM for each mean

for (i in 1:length(count_files_F)) {
  
  data[[i]]=unique(data[[i]][,c(3,14)])
  data[[i]]=data[[i]][complete.cases(data[[i]]),]
  data[[i]]=cbind(
    dCt=aggregate(.~Target.Name,data=data[[i]],mean)[,2]
    , SD.dCt=aggregate(.~Target.Name,data=data[[i]],sd)[,2]
    , SE.dCt=aggregate(.~Target.Name,data=data[[i]],std.error)[,2]
  )
}
# cal = calibrator (endogenous control)
cal.mean=data[[1]][,1]
cal.SD=data[[1]][,2]
cal.SE=data[[1]][,3]

DF=as.data.frame(rbind(data[[1]],data[[2]],data[[3]],data[[4]]))
DF$ddCt=DF$dCt-cal.mean
DF$sd.ddCt=sqrt( (DF$SD.dCt^2) + (cal.SD^2))
DF$se.ddCt=sqrt( (DF$SE.dCt^2) + (cal.SE^2))

DF$RQ=2^(-DF$ddCt)
DF$RQmin=2^(-(DF$ddCt + (DF$SD.dCt * 2.92) )) #90% confidence level
DF$RQmax=2^(-(DF$ddCt - (DF$SD.dCt * 2.92) ))

DF$RQ.SDmin = 2^(- (DF$ddCt+DF$sd.ddCt)) #standard deviation
DF$RQ.SDmax = 2^(- (DF$ddCt-DF$sd.ddCt))

DF$RQ.SEmin = 2^(- (DF$ddCt+DF$se.ddCt)) #standard error of the mean
DF$RQ.SEmax = 2^(- (DF$ddCt-DF$se.ddCt))


#fix colour to miRNA to maintain consistency among figures
miRNAs=c("novel255","miRNA_2","miRNA_3","miRNA_4","miRNA_5","miRNA_6","miRNA_7")
col.palete=cbind(miRNAs,gg_color_hue(length(miRNAs)))

DF$Targets=rep(miRNAs,4)
DF$sample=c(rep("F_28",7),rep("F_35",7),rep("F_42",7),rep("F_49",7))

DF$colour=NA
for (i in 1:length(DF$Targets)) {
  DF$colour[i] = col.palete[which(col.palete[,1] == DF$Targets[i]),2]
}

DF=DF[order(DF$Targets),]
# keep only novel255
DF_f=DF[which(DF$Targets == "novel255"),]



# Female t-test -----------------------------------------------------------

# c(1:9) selects rows for Target 1 == novel255; c(2,3,14) selects Sample, Target and delta-Ct. F_28_1, F_28_2, F_28_3 are biological replicates

F_28_255 = data_F[[1]][c(1:9),c(2,3,14)] 
# remove duplication from techincal replicates and keep only the delta-Ct values
F_28_255 = unique(F_28_255[,3])  

F_49_255 = data_F[[4]][c(1:9),c(2,3,14)] 
# remove duplication from techincal replicates and keep only the delta-Ct values
F_49_255 = unique(F_49_255[,3])  

f.test = t.test(F_28_255, F_49_255, alternative = "g")


# male worm time course --------------------------------------------------------

count_files_M = count_files[grep("M",count_files)] #keep only xls files
data_M = lapply(count_files_M, read.xls, sheet=2, skip=6, nrow=79, header=T)
data = data_M #this loop gets the mean of the delta-Ct for the biological replicates, and calculates the SEM for each mean
for (i in 1:length(count_files_M)) {
  data[[i]]=unique(data[[i]][,c(3,14)])
  data[[i]]=data[[i]][complete.cases(data[[i]]),]
  data[[i]]=cbind(
    dCt=aggregate(.~Target.Name,data=data[[i]],mean)[,2]
    , SD.dCt=aggregate(.~Target.Name,data=data[[i]],sd)[,2]
    , SE.dCt=aggregate(.~Target.Name,data=data[[i]],std.error)[,2]
  )
}

cal.mean=data[[1]][,1]
cal.SD=data[[1]][,2]
cal.SE=data[[1]][,3]

DF=as.data.frame(rbind(data[[1]],data[[2]],data[[3]],data[[4]]))
DF$ddCt=DF$dCt-cal.mean
DF$sd.ddCt=sqrt( (DF$SD.dCt^2) + (cal.SD^2))
DF$se.ddCt=sqrt( (DF$SE.dCt^2) + (cal.SE^2))

DF$RQ=2^(-DF$ddCt)
DF$RQmin=2^(-(DF$ddCt + (DF$SD.dCt * 2.92) )) #90% confidence level
DF$RQmax=2^(-(DF$ddCt - (DF$SD.dCt * 2.92) ))

DF$RQ.SDmin = 2^(- (DF$ddCt+DF$sd.ddCt))
DF$RQ.SDmax = 2^(- (DF$ddCt-DF$sd.ddCt))

DF$RQ.SEmin = 2^(- (DF$ddCt+DF$se.ddCt))
DF$RQ.SEmax = 2^(- (DF$ddCt-DF$se.ddCt))

#fix colour to miRNA to maintain consistency among figures
miRNAs # was previoulsy declared

DF$Targets=rep(miRNAs,4)
DF$sample=c(rep("M_28",7),rep("M_35",7),rep("M_42",7),rep("M_49",7))

DF$colour=NA
for (i in 1:length(DF$Targets)) {
  DF$colour[i] = col.palete[which(col.palete[,1] == DF$Targets[i]),2]
}

DF_m=DF[which(DF$Targets == "novel255"),]



# Male t-test -----------------------------------------------------------

# c(1:9) selects rows for Target 1 == novel255; c(2,3,14) selects Sample, Target and delta-Ct. F_28_1, F_28_2, F_28_3 are biological replicates

M_28_255 = data_M[[1]][c(1:9),c(2,3,14)] 
# remove duplication from techincal replicates and keep only the delta-Ct values
M_28_255 = unique(M_28_255[,3])  

M_49_255 = data_M[[4]][c(1:9),c(2,3,14)] 
# remove duplication from techincal replicates and keep only the delta-Ct values
M_49_255 = unique(M_49_255[,3])  

m.test = t.test(M_28_255, M_49_255, alternative = "g")


# make barplot ------------------------------------------------------------



# prepare data for plot

df1 = rbind(DF_f, DF_m)
df1$sample = factor(gsub("\\w_(\\d\\d)", "\\1", df1$sample ))
df1$sex = factor(c(rep("Female",4),rep("Male",4)))

# plot

ggplot(df1, aes(x=as.factor(sample), y=RQ, fill=sex)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=RQ.SEmax, ymax=RQ.SEmin), width=.2,position=position_dodge(.9)) +
  xlab("Time of perfusion after infection (days)") +
  ylab("Fold change") +
  facet_wrap(~ sex)

#ggsave("qPCR_time_course_novel255.pdf")

m.test
f.test
