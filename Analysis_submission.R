#This code is used to generate the graphs used in the paper
#In general graphs are exported to an output folder in .pdf format, 
#from where they are further processed into the finalized figures using a vector graphics editing software such as InkScape

#R version used was version 3.6.1 (2019-07-05) -- "Action of the Toes" within RStudio 1.2.5001

#Required are an input folder at "/input" with the file containing the raw measurement data 
#"Raw_Data_for_analysis.csv" with annotations, the raw data for the parallelism and Quality Control performance 
#"Parallelism.csv.txt", "QC_IgG.csv.txt" and "QC_IgA.csv.txt"

#All export commands such as "write.table" or "pdf" have been disabled by addition of "#" at the start of the line.
#Plots will therefore appear in RStudio instead of exporting and tables will be in the environment
#Each subsection corresponds to one figure in the paper

#For column names SARS2 is the abbreviation used instead of SARS-CoV-2 and Euro stands for Euroimmun



#####
#####Library####
#Used for functions regarding colour and colourscales
library(fields)
library(RColorBrewer)
library(gplots)

#Used for data depiction
library(beeswarm)
library(eulerr)


#####
#####Read in of Data####
Data_full <- read.csv("input/Raw_data_for_analysis.csv.txt",h=T,sep=",")

#split IgG from IgA Data
Data_IgG <- Data_full[,1:28]
names(Data_IgG) <- gsub("IgG_","",names(Data_IgG))
Data_IgA <- Data_full[,c(1:10,29:46)]
names(Data_IgA) <- gsub("IgA_","",names(Data_IgA))



#Annotate assay classifiers for full data set based on COs 
#ROC analysis was performed and COs were set at the MFI values specified below 
#with focus on maximum specificity with acceptable sensitivity

#cl for IgG and IgA
cl_IgG <- Data_IgG$`SARS2_RBD`>450&Data_IgG$`SARS2_Spike_Trimer`>3000
cl_IgA <- Data_IgA$`SARS2_RBD`>250&Data_IgA$`SARS2_Spike_Trimer`>400
#Combined cl combining information from IgG and IgA with either IgG CO or S/CO>2 for IgA
cl <- (Data_IgA$`SARS2_RBD`>2*250&Data_IgA$`SARS2_Spike_Trimer`>2*400)|
  (Data_IgG$`SARS2_RBD`>450&Data_IgG$`SARS2_Spike_Trimer`>3000)


#####
#####Paper Figure 1 + Table 1 (+ Supplementary Fig.3) - Comparison to Commercial tests####


#The sample set for Figure 1 is a reduced set of only those samples 
#for which antibody tests from commercial vendors were also performed
F1_IgG <- cbind(Data_IgG,Data_full[,47:54])
F1_IgG <- F1_IgG[!is.na(F1_IgG$Signal_Roche),]#only retain samples for which data of commercial assays is available

F1_IgA <- cbind(Data_IgA,Data_full[,47:54])
F1_IgA <- F1_IgA[!is.na(F1_IgA$Signal_Roche),]#only retain samples for which data of commercial assays is available

#Should be 277 samples



#Fig1, Panel a - Boxplots for RBD & S

#which AGs to be displayed
AGs <- c("SARS2_Spike_Trimer","SARS2_RBD")

#Plotting commands
#pdf(paste("output/Fig1/Fig1a.pdf",sep=""),8,6)
par(mfrow=c(1,2),mar=c(5,4.25,3,4.25))
for (i in AGs){
  boxplot(log10(F1_IgG[F1_IgG$COVID19_infection=="+",i]),
          log10(F1_IgG[F1_IgG$COVID19_infection=="-",i]),
          col=c("salmon","deepskyblue2"),xaxt="n",yaxt="n",boxwex=0.7,xlim=c(0.5,5.5))
  #add axis on the left based on non-log values
  axis(2,at=log10(c(50,100,200,500,1000,2000,5000,10000,20000,50000)),labels=c(50,100,200,500,1000,2000,5000,10000,20000,50000))
  mtext("MFI(IgG)",2,2.2,cex=1.2)
  if (i=="SARS2_Spike_Trimer"){lines(c(0.5,2.5),log10(c(3000,3000)),lwd=2,col="#303030AA")#draw cutoff across boxes
  }else if(i=="SARS2_RBD"){lines(c(0.5,2.5),log10(c(450,450)),lwd=2,col="#303030AA")#draw cutoff across boxes
  }
  par(new=T)
  boxplot(log10(F1_IgA[F1_IgA$COVID19_infection=="+",i]),
          log10(F1_IgA[F1_IgA$COVID19_infection=="-",i]),
          col=c("salmon","deepskyblue2"),xaxt="n",yaxt="n",boxwex=0.7,xlim=c(0.5,5.5),at=c(4,5))
  #adds axis on the right for set non-log values
  axis(4,at=log10(c(50,100,200,500,1000,2000,5000,10000,20000,50000)),labels=c(50,100,200,500,1000,2000,5000,10000,20000,50000))
  mtext("MFI(IgA)",4,2.2,cex=1.2)
  if (i=="SARS2_Spike_Trimer"){lines(c(3.5,5.5),log10(c(400,400)),lwd=2,col="#303030AA")#draw cutoff across boxes
  }else if(i=="SARS2_RBD"){lines(c(3.5,5.5),log10(c(250,250)),lwd=2,col="#303030AA")#draw cutoff across boxes
  }
  axis(1,at=c(1,2,4,5),labels=F)
  mtext(c("+","-","+","-"),side=1,line=0.75,at=c(1,2,4,5),cex=1)
  mtext(paste("n=",c(nrow(F1_IgG[F1_IgG$COVID19_infection=="+",]),
                     nrow(F1_IgG[F1_IgG$COVID19_infection=="-",]),
                     nrow(F1_IgA[F1_IgA$COVID19_infection=="+",]),
                     nrow(F1_IgA[F1_IgA$COVID19_infection=="-",])),sep=""),
        side=1,line=-0.875,at=c(1,2,4,5),cex=0.8)
  mtext(c("COVID-19"),side=1,line=0.7,at=c(0.1),cex=0.7)
  mtext(c("IgG","IgA"),side=1,line=2,at=c(1.5,4.5),cex=1)
  mtext(i,side=3,line=1,cex=1)
}
#dev.off()
rm(i,AGs)



#Fig1, Panel b was created from the Raw_Data Table in excel and processed further with InkScape


#Table 1 is showing test performance, it was created in Excel, except for the Confidence intervals
#Calculation of Clopper-Pearson 95% Confidence intervals

#MultiCoV-AB IgG
#Sensitivity
binom.test(x=181,n=205)
#Specificity
binom.test(x=72,n=72)

#MultiCoV-AB IgA
#Sensitivity
binom.test(x=158,n=205)
#Specificity
binom.test(x=72,n=72)

#Roche
#Sensitivity
binom.test(x=173,n=205)
#Specificity
binom.test(x=72,n=72)

#Siemens
#Sensitivity
binom.test(x=170,n=205)
#Specificity
binom.test(x=72,n=72)

#Euroimmun IgG
#Sensitivity
binom.test(x=164,n=205)
#Specificity
binom.test(x=70,n=72)

#Euroimmun IgA
#Sensitivity
binom.test(x=157,n=205)
#Specificity
binom.test(x=62,n=72)



#Supplementary Fig. 3 (all 4 Panels)
#Plots comparing MultiCoV-Ab Spike Trimer Signals to commercials assays

#Spike Trimer vs Roche
#pdf(paste("output/Suppl/S_vs_Roche.pdf",sep=""),6,4.5)
par(mar=c(5,4,4,2),mfrow=c(1,1))
plot(F1_IgG[F1_IgG$COVID19_infection=="-","SARS2_Spike_Trimer"]/3000,
     F1_IgG[F1_IgG$COVID19_infection=="-","Signal_Roche"],
     log="xy",xlim=range(F1_IgG$`SARS2_Spike_Trimer`/3000),ylim=range(F1_IgG$Signal_Roche),
     main="Comparison S (IgG) vs Roche (total Ig)",
     xlab="S/CO SARS2_Spike_Trimer IgG",ylab="S/CO Roche total Ig")
points(F1_IgG[F1_IgG$COVID19_infection=="+","SARS2_Spike_Trimer"]/3000,
       F1_IgG[F1_IgG$COVID19_infection=="+","Signal_Roche"],col="red")
abline(h=1)
abline(v=1)
#Plot Legend
legend(10^par("usr")[1],10^par("usr")[4],
       c(paste("COVID-19 + (n=",nrow(F1_IgG[F1_IgG$COVID19_infection=="+",]),")",sep=""),
       paste("COVID-19 - (n=",nrow(F1_IgG[F1_IgG$COVID19_infection=="-",]),")",sep="")),
       pch=1,cex=0.8,col=c("red","black"))
#dev.off()

#Spike Trimer vs Siemens
#pdf(paste("output/Suppl/S_vs_Siemens.pdf",sep=""),6,4.5)
par(mar=c(5,4,4,2),mfrow=c(1,1))
plot(F1_IgG[F1_IgG$COVID19_infection=="-","SARS2_Spike_Trimer"]/3000,
     F1_IgG[F1_IgG$COVID19_infection=="-","Signal_Siemens"],
     log="xy",xlim=range(F1_IgG$`SARS2_Spike_Trimer`/3000),ylim=range(F1_IgG$Signal_Siemens),
     main="Comparison S (IgG) vs Siemens (total Ig)",
     xlab="S/CO SARS2_Spike_Trimer IgG",ylab="S/CO Siemens total Ig")
points(F1_IgG[F1_IgG$COVID19_infection=="+","SARS2_Spike_Trimer"]/3000,
       F1_IgG[F1_IgG$COVID19_infection=="+","Signal_Siemens"],col="red")
abline(h=1)
abline(v=1)
#Plot Legend
legend(10^par("usr")[1],10^par("usr")[4],
       c(paste("COVID-19 + (n=",nrow(F1_IgG[F1_IgG$COVID19_infection=="+",]),")",sep=""),
         paste("COVID-19 - (n=",nrow(F1_IgG[F1_IgG$COVID19_infection=="-",]),")",sep="")),
       pch=1,cex=0.8,col=c("red","black"))
#dev.off()

#Spike Trimer vs Euroimmun
#pdf(paste("output/Suppl/S_vs_EuroImmun_IgG.pdf",sep=""),6,4.5)
par(mar=c(5,4,4,2),mfrow=c(1,1))
plot(F1_IgG[F1_IgG$COVID19_infection=="-","SARS2_Spike_Trimer"]/3000,
     F1_IgG[F1_IgG$COVID19_infection=="-","Signal_Euro_IgG"],
     log="xy",xlim=range(F1_IgG$`SARS2_Spike_Trimer`/3000),ylim=range(F1_IgG$Signal_Euro_IgG),
     main="Comparison S (IgG) vs EuroImmun (IgG)",
     xlab="S/CO SARS2_Spike_Trimer IgG",ylab="S/CO EuroImmun IgG")
points(F1_IgG[F1_IgG$COVID19_infection=="+","SARS2_Spike_Trimer"]/3000,
       F1_IgG[F1_IgG$COVID19_infection=="+","Signal_Euro_IgG"],col="red")
abline(h=1)
abline(v=1)
#Plot Legend
legend(10^par("usr")[1],10^par("usr")[4],
       c(paste("COVID-19 + (n=",nrow(F1_IgG[F1_IgG$COVID19_infection=="+",]),")",sep=""),
         paste("COVID-19 - (n=",nrow(F1_IgG[F1_IgG$COVID19_infection=="-",]),")",sep="")),
       pch=1,cex=0.8,col=c("red","black"))
#dev.off()

#Spike Trimer IgA vs Euroimmun IgA, Cutoff at 400
#pdf(paste("output/Suppl/S_vs_EuroImmun_IgA.pdf",sep=""),6,4.5)
par(mar=c(5,4,4,2),mfrow=c(1,1))
plot(F1_IgA[F1_IgA$COVID19_infection=="-","SARS2_Spike_Trimer"]/400,
     F1_IgA[F1_IgA$COVID19_infection=="-","Signal_Euro_IgA"],
     log="xy",xlim=range(F1_IgA$`SARS2_Spike_Trimer`/400),ylim=range(F1_IgA$Signal_Euro_IgA),
     main="Comparison S (IgA) vs EuroImmun (IgA)",
     xlab="S/CO SARS2_Spike_Trimer IgA",ylab="S/CO EuroImmun IgA")
points(F1_IgA[F1_IgA$COVID19_infection=="+","SARS2_Spike_Trimer"]/400,
       F1_IgA[F1_IgA$COVID19_infection=="+","Signal_Euro_IgA"],col="red")
abline(h=1)
abline(v=1)
#Plot Legend
legend(10^par("usr")[1],10^par("usr")[4],
       c(paste("COVID-19 + (n=",nrow(F1_IgA[F1_IgA$COVID19_infection=="+",]),")",sep=""),
         paste("COVID-19 - (n=",nrow(F1_IgA[F1_IgA$COVID19_infection=="-",]),")",sep="")),
       pch=1,cex=0.8,col=c("red","black"))
#dev.off()

#
rm(F1_IgG,F1_IgA)

#####
#####Paper Figure 2 + Table 3 - Screening results, cutoffs, comparison among SARS-CoV-2 antigens####


#For Figure 2, use all Samples
F2_IgG <- Data_IgG
F2_IgA <- Data_IgA


#Table 2 (and Supplementary Table 2) is a simple table with sample numbers which was not generated using R

#Table 3
#Create Table for Sens/Spec of the respective cutoffs
a <- data.frame(matrix(NA,ncol=5,nrow=8)) 
names(a) <- c("Cutoff","Correct Positives","Correct Negatives","%Sens","%Spec")
#How many samples in total
a[1,] <- c("All Samples",
           nrow(F2_IgG[F2_IgG$COVID19_infection=="+",]),
           nrow(F2_IgG[F2_IgG$COVID19_infection=="-",]),
           100,
           100)
#Cutoff IgG SARS2_Spike_Trimer only
a[2,] <- c("IgG S",
           nrow(F2_IgG[F2_IgG$COVID19_infection=="+"&F2_IgG$`SARS2_Spike_Trimer`>3000,]),
           nrow(F2_IgG[F2_IgG$COVID19_infection=="-"&!F2_IgG$`SARS2_Spike_Trimer`>3000,]),
           round(nrow(F2_IgG[F2_IgG$COVID19_infection=="+"&F2_IgG$`SARS2_Spike_Trimer`>3000,])/
             nrow(F2_IgG[F2_IgG$COVID19_infection=="+",])*100,2),
           round(nrow(F2_IgG[F2_IgG$COVID19_infection=="-"&!F2_IgG$`SARS2_Spike_Trimer`>3000,])/
             nrow(F2_IgG[F2_IgG$COVID19_infection=="-",])*100,2))
#Cutoff IgG SARS2_RBD only
a[3,] <- c("IgG RBD",
           nrow(F2_IgG[F2_IgG$COVID19_infection=="+"&F2_IgG$`SARS2_RBD`>450,]),
           nrow(F2_IgG[F2_IgG$COVID19_infection=="-"&!F2_IgG$`SARS2_RBD`>450,]),
           round(nrow(F2_IgG[F2_IgG$COVID19_infection=="+"&F2_IgG$`SARS2_RBD`>450,])/
                   nrow(F2_IgG[F2_IgG$COVID19_infection=="+",])*100,2),
           round(nrow(F2_IgG[F2_IgG$COVID19_infection=="-"&!F2_IgG$`SARS2_RBD`>450,])/
                   nrow(F2_IgG[F2_IgG$COVID19_infection=="-",])*100,2))
#Cutoff IgG combined
a[4,] <- c("IgG combined",
           nrow(F2_IgG[F2_IgG$COVID19_infection=="+"&cl_IgG,]),
           nrow(F2_IgG[F2_IgG$COVID19_infection=="-"&!cl_IgG,]),
           round(nrow(F2_IgG[F2_IgG$COVID19_infection=="+"&cl_IgG,])/
                   nrow(F2_IgG[F2_IgG$COVID19_infection=="+",])*100,2),
           round(nrow(F2_IgG[F2_IgG$COVID19_infection=="-"&!cl_IgG,])/
                   nrow(F2_IgG[F2_IgG$COVID19_infection=="-",])*100,2))
#Cutoff IgA SARS2_Spike_Trimer only
a[5,] <- c("IgA S",
           nrow(F2_IgA[F2_IgA$COVID19_infection=="+"&F2_IgA$`SARS2_Spike_Trimer`>400,]),
           nrow(F2_IgA[F2_IgA$COVID19_infection=="-"&!F2_IgA$`SARS2_Spike_Trimer`>400,]),
           round(nrow(F2_IgA[F2_IgA$COVID19_infection=="+"&F2_IgA$`SARS2_Spike_Trimer`>400,])/
                   nrow(F2_IgA[F2_IgA$COVID19_infection=="+",])*100,2),
           round(nrow(F2_IgA[F2_IgA$COVID19_infection=="-"&!F2_IgA$`SARS2_Spike_Trimer`>400,])/
                   nrow(F2_IgA[F2_IgA$COVID19_infection=="-",])*100,2))
#Cutoff IgA SARS2_RBD only
a[6,] <- c("IgA RBD",
           nrow(F2_IgA[F2_IgA$COVID19_infection=="+"&F2_IgA$`SARS2_RBD`>250,]),
           nrow(F2_IgA[F2_IgA$COVID19_infection=="-"&!F2_IgA$`SARS2_RBD`>250,]),
           round(nrow(F2_IgA[F2_IgA$COVID19_infection=="+"&F2_IgA$`SARS2_RBD`>250,])/
                   nrow(F2_IgA[F2_IgA$COVID19_infection=="+",])*100,2),
           round(nrow(F2_IgA[F2_IgA$COVID19_infection=="-"&!F2_IgA$`SARS2_RBD`>250,])/
                   nrow(F2_IgA[F2_IgA$COVID19_infection=="-",])*100,2))
#Cutoff IgA combined
a[7,] <- c("IgA combined",
           nrow(F2_IgA[F2_IgA$COVID19_infection=="+"&cl_IgA,]),
           nrow(F2_IgA[F2_IgA$COVID19_infection=="-"&!cl_IgA,]),
           round(nrow(F2_IgA[F2_IgA$COVID19_infection=="+"&cl_IgA,])/
                   nrow(F2_IgA[F2_IgA$COVID19_infection=="+",])*100,2),
           round(nrow(F2_IgA[F2_IgA$COVID19_infection=="-"&!cl_IgA,])/
                   nrow(F2_IgA[F2_IgA$COVID19_infection=="-",])*100,2))
#Cutoff IgG & IgA combined
a[8,] <- c("combined IgA & IgG",
           nrow(F2_IgA[F2_IgA$COVID19_infection=="+"&cl,]),
           nrow(F2_IgA[F2_IgA$COVID19_infection=="-"&!cl,]),
           round(nrow(F2_IgA[F2_IgA$COVID19_infection=="+"&cl,])/
                   nrow(F2_IgA[F2_IgA$COVID19_infection=="+",])*100,2),
           round(nrow(F2_IgA[F2_IgA$COVID19_infection=="-"&!cl,])/
                   nrow(F2_IgA[F2_IgA$COVID19_infection=="-",])*100,2))

#Annotate Clopper-Pearson 95% Confidence intervals and PPV/NPV at 3% seroprevalence
Sens_CI <- rep("",8)
Spec_CI <- rep("",8)
PPV_3 <- rep("",8)
NPV_3 <- rep("",8) 
for (i in 1:8){
  Sens_CI[i] <- paste("(",round(binom.test(as.numeric(a[i,"Correct Positives"]),as.numeric(a[1,"Correct Positives"]))$conf[1]*100,1),
                   " % - ",
                   round(binom.test(as.numeric(a[i,"Correct Positives"]),as.numeric(a[1,"Correct Positives"]))$conf[2]*100,1),
                   " %)",sep="")
  Spec_CI[i] <- paste("(",round(binom.test(as.numeric(a[i,"Correct Negatives"]),as.numeric(a[1,"Correct Negatives"]))$conf[1]*100,1),
                   " % - ",
                   round(binom.test(as.numeric(a[i,"Correct Negatives"]),as.numeric(a[1,"Correct Negatives"]))$conf[2]*100,1),
                   " %)",sep="")
  PPV_3[i] <- round((3*as.numeric(a[i,"%Sens"]))/((3*as.numeric(a[i,"%Sens"]))+(97*(100-as.numeric(a[i,"%Spec"]))))*100,1)
  NPV_3[i] <- round((97*as.numeric(a[i,"%Spec"]))/((97*as.numeric(a[i,"%Spec"]))+(3*(100-as.numeric(a[i,"%Sens"]))))*100,1)
  }
a$Sens_CI <- Sens_CI
a$Spec_CI <- Spec_CI
a$PPV_3 <- PPV_3
a$NPV_3 <- NPV_3


#write.table(a,"output/Tables/Table3.csv",sep=",",row.names=F,col.names=T)
rm(a,Sens_CI,Spec_CI,PPV_3,NPV_3,i)







#Figure 2, Panel a - Displays the necessity of a combined cutoff - Spike + RBD
#pdf(paste("output/Fig2/Fig2a.pdf",sep=""),6,4.5)
par(mar=c(5,4,4,2),mfrow=c(1,1))
plot(F2_IgG[F2_IgG$COVID19_infection=="-","SARS2_Spike_Trimer"]/3000,
     F2_IgG[F2_IgG$COVID19_infection=="-","SARS2_RBD"]/450,
     log="xy",cex=1,xlim=range(F2_IgG$`SARS2_Spike_Trimer`/3000),ylim=range(F2_IgG$`SARS2_RBD`/450),
     main="Spike Trimer (IgG) vs Spike RBD (IgG)",
     xlab="S/CO SARS-CoV-2 Spike Trimer IgG",ylab="S/CO SARS-CoV-2 Spike RBD IgG")
#Cutoffs
abline(h=1)
abline(v=1)
#RED=All COVID +
points(F2_IgG[F2_IgG$COVID19_infection=="+","SARS2_Spike_Trimer"]/3000,
       F2_IgG[F2_IgG$COVID19_infection=="+","SARS2_RBD"]/450,col="red",cex=1)
#Plot Legend
legend(10^par("usr")[1],10^par("usr")[4],
       c(paste("SARS-CoV-2 + (n=",nrow(F2_IgG[F2_IgG$COVID19_infection=="+",]),")",sep=""),
         paste("SARS-CoV-2 - (n=",nrow(F2_IgG[F2_IgG$COVID19_infection=="-",]),")",sep="")),
       pch=1,col=c("red","black"),cex=0.8)
#dev.off()





#Panel b - IgA Response as supplement to IgG - SARS2_Spike_Trimer vs SARS2_RBD
#pdf(paste("output/Fig2/Fig2b.pdf",sep=""),6,4.5)
par(mar=c(5,4,4,2),mfrow=c(1,1))
plot(F2_IgA[F2_IgA$COVID19_infection=="-","SARS2_Spike_Trimer"]/400,
     F2_IgA[F2_IgA$COVID19_infection=="-","SARS2_RBD"]/250,
     log="xy",cex=1,xlim=range(F2_IgA$`SARS2_Spike_Trimer`/400),ylim=range(F2_IgA$`SARS2_RBD`/250),
     main="Spike Trimer (IgA) vs Spike RBD (IgA)",
     xlab="S/CO SARS-CoV-2 Spike Trimer IgA",ylab="S/CO SARS-CoV-2 Spike RBD IgA")
#Cutoffs
#Regular IgA cutoff as dashed line
abline(h=1,lty=2)
abline(v=1,lty=2)
#S/CO>2 cutoff for combined classifier = Additional classification of strongly positive IgA samples only
abline(h=2)
abline(v=2)
#CYAN= points classified by IgG = redundant information through IgA
points(F2_IgA[cl_IgG,"SARS2_Spike_Trimer"]/400,
       F2_IgA[cl_IgG,"SARS2_RBD"]/250,col="cyan",cex=1)
#RED= points not classified though IgG
points(F2_IgA[F2_IgA$COVID19_infection=="+"&!cl_IgG,"SARS2_Spike_Trimer"]/400,
       F2_IgA[F2_IgA$COVID19_infection=="+"&!cl_IgG,"SARS2_RBD"]/250,col="red",cex=1)
#Plot Legend
legend(10^par("usr")[1],10^par("usr")[4],
       c(paste("SARS-CoV-2 + IgG + (n=",nrow(F2_IgA[F2_IgA$COVID19_infection=="+"&cl_IgG,]),")",sep=""),
         paste("SARS-CoV-2 + IgG - (n=",nrow(F2_IgA[F2_IgA$COVID19_infection=="+"&!cl_IgG,]),")",sep=""),
         paste("SARS-CoV-2 - IgG - (n=",nrow(F2_IgA[F2_IgA$COVID19_infection=="-",]),")",sep="")),
       pch=1,col=c("cyan","red","black"),cex=0.8,bg="white")
#dev.off()







#Panel c - Comparison of Spike subdomains (only IgG)
#Display S2 vs S1:

#pdf(paste("output/Fig2/Fig2c.pdf",sep=""),6,4.5)
#IgG
par(mar=c(5,4,4,2),mfrow=c(1,1))
plot(F2_IgG[F2_IgG$COVID19_infection=="-","SARS2_S2"],
     F2_IgG[F2_IgG$COVID19_infection=="-","SARS2_S1"],
     log="xy",cex=1,xlim=range(F2_IgG$`SARS2_S2`),ylim=range(F2_IgG$`SARS2_S1`),
     main="Spike S2 (IgG) vs Spike S1 (IgG)",
     xlab="MFI SARS-CoV-2 Spike S2 IgG",ylab="MFI SARS-CoV-2 Spike S1 IgG")
#CYAN=All SARS-CoV-2 infected & Classified positive by IgG/IgA combined (denoted as Spike +)
points(F2_IgG[F2_IgG$COVID19_infection=="+"&cl,"SARS2_S2"],
       F2_IgG[F2_IgG$COVID19_infection=="+"&cl,"SARS2_S1"],col="cyan",cex=1)
#RED=All COVID + & Classified -
points(F2_IgG[F2_IgG$COVID19_infection=="+"&!cl,"SARS2_S2"],
       F2_IgG[F2_IgG$COVID19_infection=="+"&!cl,"SARS2_S1"],col="red",cex=1)
#Plot Legend
legend(10^par("usr")[1],10^par("usr")[4],
       c(paste("SARS-CoV-2 + S + (n=",nrow(F2_IgG[F2_IgG$COVID19_infection=="+"&cl,]),")",sep=""),
         paste("SARS-CoV-2 + S - (n=",nrow(F2_IgG[F2_IgG$COVID19_infection=="+"&!cl,]),")",sep=""),
         paste("SARS-CoV-2 - S - (n=",nrow(F2_IgG[F2_IgG$COVID19_infection=="-",]),")",sep="")),
       pch=1,col=c("cyan","red","black"),cex=0.8)
#dev.off()





#Panel d - Nucleocapsid antigens (IgG)

#Display N vs NTD and highlight Spike positive samples
#pdf(paste("output/Fig2/Fig2d.pdf",sep=""),6,4.5)
par(mar=c(5,4,4,2),mfrow=c(1,1))
plot(F2_IgG[F2_IgG$COVID19_infection=="-","SARS2_N"],
     F2_IgG[F2_IgG$COVID19_infection=="-","SARS2_N_NTD"],
     log="xy",cex=1,xlim=range(F2_IgG$`SARS2_N`),ylim=range(F2_IgG$`SARS2_N_NTD`),
     main="Nucleocapsid (IgG) vs Nucleocapsid NTD (IgG)",
     xlab="MFI SARS-CoV-2 Nucleocapsid IgG",ylab="MFI SARS-CoV-2 Nucleocapsid NTD IgG")
#CYAN=All COVID + & Classified +
points(F2_IgG[F2_IgG$COVID19_infection=="+"&cl,"SARS2_N"],
       F2_IgG[F2_IgG$COVID19_infection=="+"&cl,"SARS2_N_NTD"],col="cyan",cex=1)
#RED=All COVID + & Classified -
points(F2_IgG[F2_IgG$COVID19_infection=="+"&!cl,"SARS2_N"],
       F2_IgG[F2_IgG$COVID19_infection=="+"&!cl,"SARS2_N_NTD"],col="red",cex=1)
#Plot Legend
legend(10^par("usr")[1],10^par("usr")[4],
       c(paste("SARS-CoV-2 + S + (n=",nrow(F2_IgG[F2_IgG$COVID19_infection=="+"&cl,]),")",sep=""),
         paste("SARS-CoV-2 + S - (n=",nrow(F2_IgG[F2_IgG$COVID19_infection=="+"&!cl,]),")",sep=""),
         paste("SARS-CoV-2 - S - (n=",nrow(F2_IgG[F2_IgG$COVID19_infection=="-",]),")",sep="")),
       pch=1,col=c("cyan","red","black"),cex=0.8)
#dev.off()


#
rm(F2_IgG,F2_IgA)


#####
#####Paper Figure 3 - SARS-CoV-2 Samples Trends#######


#For Figure 3, use all Samples
F3_IgG <- Data_IgG
F3_IgA <- Data_IgA



#Fig3, Panel a - Ab response of patients with multiple samples

#Retain only samples from a time series for this panel using patient IDatients 
IgA <- F3_IgA[!is.na(F3_IgA$Time_series_PatientID),]
IgG <- F3_IgG[!is.na(F3_IgG$Time_series_PatientID),]

#Display time course of response for all 6 SARS2 antigens
AGs <- c("SARS2_Spike_Trimer","SARS2_RBD","SARS2_S1","SARS2_S2","SARS2_N","SARS2_N_NTD")
col <- colorRampPalette(brewer.pal(11,"Spectral")[c(8:11,1:3)])(5)
#pdf(paste("output/Fig3/Fig3a.pdf",sep=""),3,6.8)
par(mfrow=c(3,1),mar=c(5,5,4,4))
for (j in AGs){
  plot(1,1,cex=0,xlim=range(as.numeric(IgG$COVID19_dT)),ylim=range(IgG[,j]),
       xlab="dT",ylab=paste(j),log="y",main=paste("IgG",j))
  c <- 0
  for (i in unique(IgG$Time_series_PatientID)){
    c <- c+1#counter for colors
    points(IgG[IgG$Time_series_PatientID==i,"COVID19_dT"],IgG[IgG$Time_series_PatientID==i,j],pch=1,cex=1.5,col=col[c])
    lines(IgG[IgG$Time_series_PatientID==i,"COVID19_dT"],IgG[IgG$Time_series_PatientID==i,j],lty=1,col=col[c])
  }
  plot(1,1,cex=0,xlim=range(as.numeric(IgA$COVID19_dT)),ylim=range(IgA[,j]),
       xlab="dT",ylab=paste(j),log="y",main=paste("IgA",j))
  c <- 0
  for (i in unique(IgG$Time_series_PatientID)){
    c <- c+1#counter for colors
    points(IgA[IgA$Time_series_PatientID==i,"COVID19_dT"],IgA[IgA$Time_series_PatientID==i,j],pch=1,cex=1.5,col=col[c])
    lines(IgA[IgA$Time_series_PatientID==i,"COVID19_dT"],IgA[IgA$Time_series_PatientID==i,j],lty=1,col=col[c])
  }
}
#dev.off()
rm(AGs,c,i,j,IgG,IgA,col)
#re-sort plots and append legend manually in InkScape








#Fig3, Panel b - Effect of Hospitalisation

AGs <- c("SARS2_Spike_Trimer","SARS2_RBD","SARS2_N")
#pdf(paste("output/Fig3/Fig3b.pdf",sep=""),8,6)
par(mfrow=c(1,2),mar=c(5,4.25,3,4.25))
for (i in AGs){
  boxplot(log10(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="+",i]),
          log10(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="-",i]),
          col=c("salmon","deepskyblue2"),xaxt="n",yaxt="n",boxwex=0.7,xlim=c(0.5,5.5))
  #add axis on the left based on non-log values
  axis(2,at=log10(c(50,100,200,500,1000,2000,5000,10000,20000,50000)),labels=c(50,100,200,500,1000,2000,5000,10000,20000,50000))
  mtext("MFI(IgG)",2,2.2,cex=1.2)
  if (i=="SARS2_Spike_Trimer"){lines(c(0.5,2.5),log10(c(3000,3000)),lwd=2,col="#303030AA")#draw cutoff across boxes
  }else if(i=="SARS2_RBD"){lines(c(0.5,2.5),log10(c(450,450)),lwd=2,col="#303030AA")#draw cutoff across boxes
  }
  par(new=T)
  boxplot(log10(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="+",i]),
          log10(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="-",i]),
          col=c("salmon","deepskyblue2"),xaxt="n",yaxt="n",boxwex=0.7,xlim=c(0.5,5.5),at=c(4,5))
  #adds axis on the right for set non-log values
  axis(4,at=log10(c(50,100,200,500,1000,2000,5000,10000,20000,50000)),labels=c(50,100,200,500,1000,2000,5000,10000,20000,50000))
  mtext("MFI(IgA)",4,2.2,cex=1.2)
  if (i=="SARS2_Spike_Trimer"){lines(c(3.5,5.5),log10(c(400,400)),lwd=2,col="#303030AA")#draw cutoff across boxes
  }else if(i=="SARS2_RBD"){lines(c(3.5,5.5),log10(c(250,250)),lwd=2,col="#303030AA")#draw cutoff across boxes
  }
  axis(1,at=c(1,2,4,5),labels=F)
  mtext(c("+","-","+","-"),side=1,line=0.75,at=c(1,2,4,5),cex=1)
  mtext(paste("n=",c(nrow(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="+",]),
                     nrow(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="-",]),
                     nrow(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="+",]),
                     nrow(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="-",])),sep=""),
        side=1,line=-0.875,at=c(1,2,4,5),cex=0.8)
  mtext(c("Hospitalisation"),side=1,line=0.7,at=c(0.1),cex=0.7)
  mtext(c("IgG","IgA"),side=1,line=2,at=c(1.5,4.5),cex=1)
  mtext(i,side=3,line=1,cex=1)
}
#dev.off()
rm(i,AGs)








#Fig3, Panel c - Effect of Age

AGs <- c("SARS2_Spike_Trimer","SARS2_RBD","SARS2_N")
#pdf(paste("output/Fig3/Fig3c.pdf",sep=""),10,6)
par(mfrow=c(1,2),mar=c(5,4.25,3,4.25))
for (i in AGs){
  boxplot(log10(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-",i]),
          log10(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59",i]),
          log10(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+",i]),
          col=c("gray90","gray60","gray30"),xaxt="n",yaxt="n",boxwex=0.7,xlim=c(0.5,7.5))
  axis(2,at=log10(c(50,100,200,500,1000,2000,5000,10000,20000,50000)),labels=c(50,100,200,500,1000,2000,5000,10000,20000,50000))
  mtext("MFI(IgG)",2,2.2,cex=1.2)
  if (i=="SARS2_Spike_Trimer"){lines(c(0.5,3.5),log10(c(3000,3000)),lwd=2,col="#303030AA")#draw cutoff across boxes
  }else if(i=="SARS2_RBD"){lines(c(0.5,3.5),log10(c(450,450)),lwd=2,col="#303030AA")#draw cutoff across boxes
  }
  par(new=T)
  boxplot(log10(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-",i]),
          log10(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59",i]),
          log10(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+",i]),
          col=c("gray90","gray60","gray30"),xaxt="n",yaxt="n",boxwex=0.7,xlim=c(0.5,7.5),at=c(5,6,7))
  axis(4,at=log10(c(50,100,200,500,1000,2000,5000,10000,20000,50000)),labels=c(50,100,200,500,1000,2000,5000,10000,20000,50000))
  mtext("MFI(IgA)",4,2.2,cex=1.2)
  if (i=="SARS2_Spike_Trimer"){lines(c(4.5,7.5),log10(c(400,400)),lwd=2,col="#303030AA")#draw cutoff across boxes
  }else if(i=="SARS2_RBD"){lines(c(4.5,7.5),log10(c(250,250)),lwd=2,col="#303030AA")#draw cutoff across boxes
  }
  axis(1,at=c(1,2,4,5),labels=F)
  mtext(c("<=40","40-59",">=60","<=40","40-59",">=60"),side=1,line=0.75,at=c(1,2,3,5,6,7),cex=0.8)
  mtext(paste("n=",c(nrow(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-",]),
                     nrow(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59",]),
                     nrow(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+",]),
                     nrow(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-",]),
                     nrow(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59",]),
                     nrow(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+",])),sep=""),
        side=1,line=-0.875,at=c(1,2,3,5,6,7),cex=0.8)
  mtext(c("Age"),side=1,line=0.7,at=c(0.1),cex=0.7)
  mtext(c("IgG","IgA"),side=1,line=2,at=c(2,6),cex=1)
  mtext(i,side=3,line=1,cex=1)
}
#dev.off()
rm(i,AGs)





#Statistics for Panel b & c - Mann-Whitney U Tests comparing two boxes each
#For export into table and manual annotation into Figure 3 plots afterwards

#Create Table and fill in each cell, then export
MWU <- data.frame(matrix(0,ncol=7,nrow=4))
names(MWU) <- c("Param","Spike IgG","Spike IgA","RBD IgG","RBD IgA","N IgG","N IgA")
MWU$Param <- c("Hosp. vs Non-Hosp.","Young vs Middle","Middle vs Old","Young vs Old")

#Hospitalisation has an effect on Signals
#Spike Trimer IgG
MWU[1,"Spike IgG"] <- wilcox.test(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="+","SARS2_Spike_Trimer"],
                                  F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="-","SARS2_Spike_Trimer"],
                                  alternative="t")$p.value
#Spike Trimer IgA
MWU[1,"Spike IgA"] <- wilcox.test(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="+","SARS2_Spike_Trimer"],
                                  F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="-","SARS2_Spike_Trimer"],
                                  alternative="t")$p.value
#RBD IgG
MWU[1,"RBD IgG"] <- wilcox.test(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="+","SARS2_RBD"],
                                F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="-","SARS2_RBD"],
                                alternative="t")$p.value
#RBD IgA
MWU[1,"RBD IgA"] <- wilcox.test(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="+","SARS2_RBD"],
                                F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="-","SARS2_RBD"],
                                alternative="t")$p.value
#N IgG
MWU[1,"N IgG"] <- wilcox.test(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="+","SARS2_N"],
                              F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="-","SARS2_N"],
                              alternative="t")$p.value
#N IgA
MWU[1,"N IgA"] <- wilcox.test(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="+","SARS2_N"],
                              F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$COVID19_hospitalized=="-","SARS2_N"],
                              alternative="t")$p.value
#Age
#Spike Trimer IgG
#Young vs Middle
MWU[2,"Spike IgG"] <- wilcox.test(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-","SARS2_Spike_Trimer"],
                                  F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59","SARS2_Spike_Trimer"],
                                  alternative="t")$p.value
#Middle vs Old
MWU[3,"Spike IgG"] <- wilcox.test(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59","SARS2_Spike_Trimer"],
                                  F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+","SARS2_Spike_Trimer"],
                                  alternative="t")$p.value
#Young vs Old
MWU[4,"Spike IgG"] <- wilcox.test(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-","SARS2_Spike_Trimer"],
                                  F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+","SARS2_Spike_Trimer"],
                                  alternative="t")$p.value
#Spike Trimer IgA
#Young vs Middle
MWU[2,"Spike IgA"] <- wilcox.test(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-","SARS2_Spike_Trimer"],
                                  F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59","SARS2_Spike_Trimer"],
                                  alternative="t")$p.value
#Middle vs Old
MWU[3,"Spike IgA"] <- wilcox.test(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59","SARS2_Spike_Trimer"],
                                  F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+","SARS2_Spike_Trimer"],
                                  alternative="t")$p.value
#Young vs Old
MWU[4,"Spike IgA"] <- wilcox.test(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-","SARS2_Spike_Trimer"],
                                  F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+","SARS2_Spike_Trimer"],
                                  alternative="t")$p.value
#RBD IgG
#Young vs Middle
MWU[2,"RBD IgG"] <- wilcox.test(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-","SARS2_RBD"],
                                  F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59","SARS2_RBD"],
                                  alternative="t")$p.value
#Middle vs Old
MWU[3,"RBD IgG"] <- wilcox.test(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59","SARS2_RBD"],
                                  F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+","SARS2_RBD"],
                                  alternative="t")$p.value
#Young vs Old
MWU[4,"RBD IgG"] <- wilcox.test(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-","SARS2_RBD"],
                                  F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+","SARS2_RBD"],
                                  alternative="t")$p.value
#RBD IgA
#Young vs Middle
MWU[2,"RBD IgA"] <- wilcox.test(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-","SARS2_RBD"],
                                  F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59","SARS2_RBD"],
                                  alternative="t")$p.value
#Middle vs Old
MWU[3,"RBD IgA"] <- wilcox.test(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59","SARS2_RBD"],
                                  F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+","SARS2_RBD"],
                                  alternative="t")$p.value
#Young vs Old
MWU[4,"RBD IgA"] <- wilcox.test(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-","SARS2_RBD"],
                                  F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+","SARS2_RBD"],
                                  alternative="t")$p.value
#N IgG
#Young vs Middle
MWU[2,"N IgG"] <- wilcox.test(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-","SARS2_N"],
                                F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59","SARS2_N"],
                                alternative="t")$p.value
#Middle vs Old
MWU[3,"N IgG"] <- wilcox.test(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59","SARS2_N"],
                                F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+","SARS2_N"],
                                alternative="t")$p.value
#Young vs Old
MWU[4,"N IgG"] <- wilcox.test(F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-","SARS2_N"],
                                F3_IgG[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+","SARS2_N"],
                                alternative="t")$p.value
#N IgA
#Young vs Middle
MWU[2,"N IgA"] <- wilcox.test(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-","SARS2_N"],
                                F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59","SARS2_N"],
                                alternative="t")$p.value
#Middle vs Old
MWU[3,"N IgA"] <- wilcox.test(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="40-59","SARS2_N"],
                                F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+","SARS2_N"],
                                alternative="t")$p.value
#Young vs Old
MWU[4,"N IgA"] <- wilcox.test(F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="39-","SARS2_N"],
                                F3_IgA[F3_IgG$COVID19_infection=="+"&F3_IgG$Age=="60+","SARS2_N"],
                                alternative="t")$p.value
#write.table(MWU,"output/Fig3/Mann_Whitney_U_Hosp_Age.csv",col.names=T,row.names=F,sep=",")
rm(MWU)

#
rm(F3_IgA,F3_IgG)

#####
#####Paper Figure 4 and 5 - hCoV Analysis for cross-reactivity#####

#Data to work with for Fig4

#For Figure 4 (and 5), use all Samples
F4_IgG <- Data_IgG
F4_IgA <- Data_IgA


#Get a scaled version of the Data for clustering and to compare signals against hCoVs amongst one another
Sc_IgG <- F4_IgG
Sc_IgG[,11:28] <- scale(log(Sc_IgG[,11:28]))
Sc_IgA <- F4_IgA
Sc_IgA[,11:28] <- scale(log(Sc_IgA[,11:28]))



#Fig4, Panel a - Correlation Heatmap for all Antigens using IgG Data

x <- F4_IgG[,c("SARS2_Spike_Trimer","SARS2_RBD","SARS2_S1","SARS2_S2","SARS2_N","SARS2_N_NTD",
                "229E_N","229E_N_NTD","229E_S1",
                "NL63_N","NL63_N_NTD","NL63_S1",
                "OC43_N","OC43_N_NTD","OC43_S1",
                "HKU1_N","HKU1_N_NTD","HKU1_S1")]
x <- cor(x,method="spearman")
hmcol <- colorRampPalette(c(col2hex("red"),col2hex("white"),col2hex("blue")))(50)
#pdf("output/Fig4/Fig4a.pdf",8,8)
heatmap.2(x,main="Correlation of Antigens",scale='none',
          Colv=T,Rowv=T,dendrogram='row',
          labRow="",labCol=colnames(x),
          col=hmcol,trace='none',margins=c(10,10),breaks=51,symbreaks=T,cexCol=1,cexRow=1,
          density.info="none",key.xlab="R(Spearman)")
#dev.off()
rm(x,hmcol)
#Further editing done in inkscape


#Supplementary Fig. 6 - Panel a - Correlation Heatmap for all Antigens using IgA Data

x <- F4_IgA[,c("SARS2_Spike_Trimer","SARS2_RBD","SARS2_S1","SARS2_S2","SARS2_N","SARS2_N_NTD",
                "229E_N","229E_N_NTD","229E_S1",
                "NL63_N","NL63_N_NTD","NL63_S1",
                "OC43_N","OC43_N_NTD","OC43_S1",
                "HKU1_N","HKU1_N_NTD","HKU1_S1")]
x <- cor(x,method="spearman")
hmcol <- colorRampPalette(c(col2hex("red"),col2hex("white"),col2hex("blue")))(50)
#pdf("output/Suppl/SupplFig6a.pdf",8,8)
heatmap.2(x,main="Correlation of Antigens",scale='none',
          Colv=T,Rowv=T,dendrogram='row',
          labRow="",labCol=colnames(x),
          col=hmcol,trace='none',margins=c(10,10),breaks=51,symbreaks=T,cexCol=1,cexRow=1,
          density.info="none",key.xlab="R(Spearman)")
#dev.off()
rm(x,hmcol)




#Fig4, Panel b & c + Supplementary Fig. 6b - Display hCoV AGs signals for COVID+ vs COVID- 
#     as Boxplots (N NTDs are shown in Supplementary as Figure 6b)

for (i in c("S1","N_NTD","N")){
  #pdf(paste("output/Fig4/Boxplots_hCOV_",i,".pdf",sep=""),8,6)
  par(mfrow=c(1,2),mar=c(5,4.25,3,4.25))
  for (j in c("NL63","229E","OC43","HKU1")){
    boxplot(log10(F4_IgG[F4_IgG$COVID19_infection=="+",paste(j,i,sep="_")]),
            log10(F4_IgG[F4_IgG$COVID19_infection=="-",paste(j,i,sep="_")]),
            col=c("salmon","deepskyblue2"),xaxt="n",yaxt="n",boxwex=0.7,xlim=c(0.5,5.5))
    #add axis on the left based on non-log values
    axis(2,at=log10(c(50,100,200,500,1000,2000,5000,10000,20000,50000)),labels=c(50,100,200,500,1000,2000,5000,10000,20000,50000))
    mtext("MFI(IgG)",2,2.2,cex=1.2)
    par(new=T)
    boxplot(log10(F4_IgA[F4_IgA$COVID19_infection=="+",paste(j,i,sep="_")]),
            log10(F4_IgA[F4_IgA$COVID19_infection=="-",paste(j,i,sep="_")]),
            col=c("salmon","deepskyblue2"),xaxt="n",yaxt="n",boxwex=0.7,xlim=c(0.5,5.5),at=c(4,5))
    #adds axis on the right for set non-log values
    axis(4,at=log10(c(50,100,200,500,1000,2000,5000,10000,20000,50000)),labels=c(50,100,200,500,1000,2000,5000,10000,20000,50000))
    mtext("MFI(IgA)",4,2.2,cex=1.2)
    axis(1,at=c(1,2,4,5),labels=F)
    mtext(c("+","-","+","-"),side=1,line=0.75,at=c(1,2,4,5),cex=1)
    mtext(paste("n=",c(nrow(F4_IgG[F4_IgG$COVID19_infection=="+",]),
                       nrow(F4_IgG[F4_IgG$COVID19_infection=="-",]),
                       nrow(F4_IgA[F4_IgA$COVID19_infection=="+",]),
                       nrow(F4_IgA[F4_IgA$COVID19_infection=="-",])),sep=""),
          side=1,line=-0.875,at=c(1,2,4,5),cex=0.8)
    mtext(c("SARS-CoV-2"),side=1,line=0.75,at=-0.75,cex=1)
    mtext(c("IgG","IgA"),side=1,line=2,at=c(1.5,4.5),cex=1)
    mtext(paste(j,i),side=3,line=1,cex=1)
  }
  #dev.off()
}
rm(i,j)
#




#Panel d & e -  Analysis of hCoV response to SARS-CoV-2 IgG Falsepositives and Falsenegatives

#Work with scaled Data

#Get target groups in data frame
#Spike False Positives
FalsePosS <- Sc_IgG[F4_IgG$COVID19_infection=="-"&F4_IgG$`SARS2_Spike_Trimer`>3000,]
#Infected but not classified
FalseNeg <- Sc_IgG[F4_IgG$COVID19_infection=="+"&!cl,]

#Beeswarm on top of boxplot
for (i in c("S1","N")){
  #pdf(paste("output/Fig4/Boxplots_FalsePosNeg_",i,".pdf",sep=""),8,5)
  par(mfrow=c(1,1),mar=c(5,4.25,3,4.25))
  beeswarm(list(FalsePosS[,paste("NL63",i,sep="_")],FalseNeg[,paste("NL63",i,sep="_")],
                FalsePosS[,paste("229E",i,sep="_")],FalseNeg[,paste("229E",i,sep="_")],
                FalsePosS[,paste("OC43",i,sep="_")],FalseNeg[,paste("OC43",i,sep="_")],
                FalsePosS[,paste("HKU1",i,sep="_")],FalseNeg[,paste("HKU1",i,sep="_")]),
           col=c("deepskyblue3","salmon"),pch=16,ylim=c(-3.5,3.5),main=paste(i),xaxt="n",at=c(1,2,3.5,4.5,6,7,8.5,9.5))
  boxplot(FalsePosS[,paste("NL63",i,sep="_")],FalseNeg[,paste("NL63",i,sep="_")],
          FalsePosS[,paste("229E",i,sep="_")],FalseNeg[,paste("229E",i,sep="_")],
          FalsePosS[,paste("OC43",i,sep="_")],FalseNeg[,paste("OC43",i,sep="_")],
          FalsePosS[,paste("HKU1",i,sep="_")],FalseNeg[,paste("HKU1",i,sep="_")],
          outline=F,add=T,lwd=0.75,col="#00000000",xaxt="n",at=c(1,2,3.5,4.5,6,7,8.5,9.5))
  abline(h=0,lty=2)
  legend(0.1,3.9,paste(c("Spike Trimer False-Positive","False-Negative")," n=",c(nrow(FalsePosS),nrow(FalseNeg)),sep=""),
         cex=1,col=c("deepskyblue3","salmon"),pch=16,bty="n")
  axis(side=1,at=c(1.5,4,6.5,9),labels=c("NL63","229E","OC43","HKU1"))
  #dev.off()
}
rm(i)


# Supplementary Fig 6c - Same analysis as above for N NTD
for (i in c("N_NTD")){
  #pdf(paste("output/Suppl/Boxplots_FalsePosNeg_",i,".pdf",sep=""),8,5)
  par(mfrow=c(1,1),mar=c(5,4.25,3,4.25))
  beeswarm(list(FalsePosS[,paste("NL63",i,sep="_")],FalseNeg[,paste("NL63",i,sep="_")],
                FalsePosS[,paste("229E",i,sep="_")],FalseNeg[,paste("229E",i,sep="_")],
                FalsePosS[,paste("OC43",i,sep="_")],FalseNeg[,paste("OC43",i,sep="_")],
                FalsePosS[,paste("HKU1",i,sep="_")],FalseNeg[,paste("HKU1",i,sep="_")]),
           col=c("deepskyblue3","salmon"),pch=16,ylim=c(-4,4),main=paste(i),xaxt="n",at=c(1,2,3.5,4.5,6,7,8.5,9.5))
  boxplot(FalsePosS[,paste("NL63",i,sep="_")],FalseNeg[,paste("NL63",i,sep="_")],
          FalsePosS[,paste("229E",i,sep="_")],FalseNeg[,paste("229E",i,sep="_")],
          FalsePosS[,paste("OC43",i,sep="_")],FalseNeg[,paste("OC43",i,sep="_")],
          FalsePosS[,paste("HKU1",i,sep="_")],FalseNeg[,paste("HKU1",i,sep="_")],
          outline=F,add=T,lwd=0.75,col="#00000000",xaxt="n",at=c(1,2,3.5,4.5,6,7,8.5,9.5))
  abline(h=0,lty=2)
  legend(0.1,3.9,paste(c("Spike Trimer False-Positive","False-Negative")," n=",c(nrow(FalsePosS),nrow(FalseNeg)),sep=""),
         cex=1,col=c("deepskyblue3","salmon"),pch=16,bty="n")
  axis(side=1,at=c(1.5,4,6.5,9),labels=c("NL63","229E","OC43","HKU1"))
  #dev.off()
}
rm(i)

#
rm(FalseNeg,FalsePosS)



#Fig5 - Grouping into hCoV high and low responders (IgG Response)


#Create High and low responder subgroups, scaled signals have to be >0 for all 4 alpha/beta-CoV antigens
#Also create the overlapping group with classification positives as well for the following analysis (ending in "P")
AlphaUp <- Sc_IgG[Sc_IgG$`NL63_S1`>0&Sc_IgG$`229E_S1`>0&
                    Sc_IgG$`NL63_N`>0&Sc_IgG$`229E_N`>0,]
AlphaUpP <- Sc_IgG[Sc_IgG$`NL63_S1`>0&Sc_IgG$`229E_S1`>0&
                     Sc_IgG$`NL63_N`>0&Sc_IgG$`229E_N`>0&cl,]
BetaUp <- Sc_IgG[Sc_IgG$`OC43_S1`>0&Sc_IgG$`HKU1_S1`>0&
                   Sc_IgG$`OC43_N`>0&Sc_IgG$`HKU1_N`>0,]
BetaUpP <- Sc_IgG[Sc_IgG$`OC43_S1`>0&Sc_IgG$`HKU1_S1`>0&
                    Sc_IgG$`OC43_N`>0&Sc_IgG$`HKU1_N`>0&cl,]
AlphaDown <- Sc_IgG[Sc_IgG$`NL63_S1`<0&Sc_IgG$`229E_S1`<0&
                      Sc_IgG$`NL63_N`<0&Sc_IgG$`229E_N`<0,]
AlphaDownP <- Sc_IgG[Sc_IgG$`NL63_S1`<0&Sc_IgG$`229E_S1`<0&
                       Sc_IgG$`NL63_N`<0&Sc_IgG$`229E_N`<0&cl,]
BetaDown <- Sc_IgG[Sc_IgG$`OC43_S1`<0&Sc_IgG$`HKU1_S1`<0&
                     Sc_IgG$`OC43_N`<0&Sc_IgG$`HKU1_N`<0,]
BetaDownP <- Sc_IgG[Sc_IgG$`OC43_S1`<0&Sc_IgG$`HKU1_S1`<0&
                      Sc_IgG$`OC43_N`<0&Sc_IgG$`HKU1_N`<0&cl,]
#MultiCoV-Ab positives
P <- Sc_IgG[cl,]



#Display group phenotypes as boxplot
#pdf(paste("output/Fig5/Fig5a.pdf",sep=""),10,8)
par(mar=c(4,3,4,2),mfrow=c(2,2))
boxplot(AlphaUp[,"NL63_S1"],AlphaUp[,"229E_S1"],
        AlphaUp[,"OC43_S1"],AlphaUp[,"HKU1_S1"],
        AlphaUp[,"NL63_N"],AlphaUp[,"229E_N"],
        AlphaUp[,"OC43_N"],AlphaUp[,"HKU1_N"],
        ylim=c(-5.5,5.5),
        col=brewer.pal(9,"Set1")[1],names=c("NL63","229E","OC43","HKU1","NL63","229E","OC43","HKU1"))
mtext(paste("a-CoV High-responders n =",nrow(AlphaUp[,])),3,1)
mtext(c("S1","N"),1,2.5,at=c(2.5,6.5))
abline(h=0,lty=2)
boxplot(AlphaDown[,"NL63_S1"],AlphaDown[,"229E_S1"],
        AlphaDown[,"OC43_S1"],AlphaDown[,"HKU1_S1"],
        AlphaDown[,"NL63_N"],AlphaDown[,"229E_N"],
        AlphaDown[,"OC43_N"],AlphaDown[,"HKU1_N"],
        ylim=c(-5.5,5.5),
        col=brewer.pal(9,"Set1")[2],names=c("NL63","229E","OC43","HKU1","NL63","229E","OC43","HKU1"))
mtext(paste("a-CoV Low-responders n =",nrow(AlphaDown[,])),3,1)
mtext(c("S1","N"),1,2.5,at=c(2.5,6.5))
abline(h=0,lty=2)
boxplot(BetaUp[,"NL63_S1"],BetaUp[,"229E_S1"],
        BetaUp[,"OC43_S1"],BetaUp[,"HKU1_S1"],
        BetaUp[,"NL63_N"],BetaUp[,"229E_N"],
        BetaUp[,"OC43_N"],BetaUp[,"HKU1_N"],
        ylim=c(-5.5,5.5),
        col=brewer.pal(9,"Set1")[3],names=c("NL63","229E","OC43","HKU1","NL63","229E","OC43","HKU1"))
mtext(paste("b-CoV High-responders n =",nrow(BetaUp[,])),3,1)
mtext(c("S1","N"),1,2.5,at=c(2.5,6.5))
abline(h=0,lty=2)
boxplot(BetaDown[,"NL63_S1"],BetaDown[,"229E_S1"],
        BetaDown[,"OC43_S1"],BetaDown[,"HKU1_S1"],
        BetaDown[,"NL63_N"],BetaDown[,"229E_N"],
        BetaDown[,"OC43_N"],BetaDown[,"HKU1_N"],
        ylim=c(-5.5,5.5),
        col=brewer.pal(9,"Set1")[4],names=c("NL63","229E","OC43","HKU1","NL63","229E","OC43","HKU1"))
mtext(paste("b-CoV Low-responders n =",nrow(BetaDown[,])),3,1)
mtext(c("S1","N"),1,2.5,at=c(2.5,6.5))
abline(h=0,lty=2)
#dev.off()


#Panel b - venn Diagrams of overlap

#Get a mixture of the group colors and "ivory" for the venn diagrams
col1 <- data.frame(as.character(brewer.pal(6, "Set1")),NA)
for (i in 1:6){col1[i,2] <- designer.colors(n=6, col=c(as.character(col1[i,1]),"ivory"))[5]}
rm(i)


#Venn Diagrams
#pdf(paste("output/Fig5/Fig5b.pdf",sep=""),8,5)
plot(euler(c("a-CoV+"=nrow(AlphaUp)-nrow(AlphaUpP),"SARS2+"=nrow(P)-nrow(AlphaUpP),"a-CoV+&SARS2+"=nrow(AlphaUpP))),
     input="disjoint",fill=c(paste(col1[1,1],"70",sep=""),paste(col1[1,2],"70",sep="")),lwd=2)
plot(euler(c("a-CoV-"=nrow(AlphaDown)-nrow(AlphaDownP),"SARS2+"=nrow(P)-nrow(AlphaDownP),"a-CoV-&SARS2+"=nrow(AlphaDownP))),
     input="disjoint",fill=c(paste(col1[2,1],"70",sep=""),paste(col1[2,2],"70",sep="")),lwd=2)
plot(euler(c("b-CoV+"=nrow(BetaUp)-nrow(BetaUpP),"SARS2+"=nrow(P)-nrow(BetaUpP),"b-CoV+&SARS2+"=nrow(BetaUpP))),
     input="disjoint",fill=c(paste(col1[3,1],"70",sep=""),paste(col1[3,2],"70",sep="")),lwd=2)
plot(euler(c("b-CoV-"=nrow(BetaDown)-nrow(BetaDownP),"SARS2+"=nrow(P)-nrow(BetaDownP),"b-CoV-&SARS2+"=nrow(BetaDownP))),
     input="disjoint",fill=c(paste(col1[4,1],"70",sep=""),paste(col1[4,2],"70",sep="")),lwd=2)
#dev.off()


#Fishers Exact Test to see if enrichment occurs, the 4 p-values are manually annotated to Fig.5 panel b in inkscape

#use formatC in order to get scientific notation
#Positives in above average alpha-CoV responders
data = matrix(data=c(nrow(AlphaUpP),nrow(AlphaUp)-nrow(AlphaUpP),nrow(P),nrow(Sc_IgG)-nrow(P)),nrow=2)
formatC(fisher.test(data,alternative="two.sided")$p.value,format="e")
#Positives in below average alpha-CoV responders
data = matrix(data=c(nrow(AlphaDownP),nrow(AlphaDown)-nrow(AlphaDownP),nrow(P),nrow(Sc_IgG)-nrow(P)),nrow=2)
formatC(fisher.test(data,alternative="two.sided")$p.value,format="e")
#Positives in above average beta-CoV responders
data = matrix(data=c(nrow(BetaUpP),nrow(BetaUp)-nrow(BetaUpP),nrow(P),nrow(Sc_IgG)-nrow(P)),nrow=2)
formatC(fisher.test(data,alternative="two.sided")$p.value,format="e")
#Positives in below average beta-CoV responders
data = matrix(data=c(nrow(BetaDownP),nrow(BetaDown)-nrow(BetaDownP),nrow(P),nrow(Sc_IgG)-nrow(P)),nrow=2)
formatC(fisher.test(data,alternative="two.sided")$p.value,format="e")

rm(AlphaDown,AlphaDownP,AlphaUp,AlphaUpP,P,BetaDown,BetaDownP,BetaUp,BetaUpP,col1,data)
#
rm(F4_IgG,F4_IgA,Sc_IgG,Sc_IgA)

#####
#####Supplementary Figure 5 - delta T vs MFI####


#pdf("output/Suppl/SupplFig5_dT_effects.pdf",9,6)
par(mfrow=c(2,3),mar=c(4,4,4,4))
for (i in c("SARS2_Spike_Trimer","SARS2_RBD","SARS2_S1","SARS2_S2","SARS2_N","SARS2_N_NTD")){
  plot(Data_IgG$COVID19_dT,Data_IgG[,i],log="y",
       main=paste("IgG -",i),xlab="Sample dT (days)",ylab=paste("Signal IgG",i,"(MFI)"),cex=0)
  points(Data_IgG[Data_IgG$COVID19_dT_type=="PCR Specimen","COVID19_dT"],
         Data_IgG[Data_IgG$COVID19_dT_type=="PCR Specimen",i],col="#FF000088")
  points(Data_IgG[Data_IgG$COVID19_dT_type=="PCR Date","COVID19_dT"],
         Data_IgG[Data_IgG$COVID19_dT_type=="PCR Date",i],col="#00FF0088")
  points(Data_IgG[Data_IgG$COVID19_dT_type=="Symptom Onset","COVID19_dT"],
         Data_IgG[Data_IgG$COVID19_dT_type=="Symptom Onset",i],col="#0000FF88")
}
for (i in c("SARS2_Spike_Trimer","SARS2_RBD","SARS2_S1","SARS2_S2","SARS2_N","SARS2_N_NTD")){
  plot(Data_IgA$COVID19_dT,Data_IgA[,i],log="y",
       main=paste("IgA -",i),xlab="Sample dT (days)",ylab=paste("Signal IgA",i,"(MFI)"),cex=0)
  points(Data_IgA[Data_IgA$COVID19_dT_type=="PCR Specimen","COVID19_dT"],
         Data_IgA[Data_IgA$COVID19_dT_type=="PCR Specimen",i],col="#FF000088")
  points(Data_IgA[Data_IgA$COVID19_dT_type=="PCR Date","COVID19_dT"],
         Data_IgA[Data_IgA$COVID19_dT_type=="PCR Date",i],col="#00FF0088")
  points(Data_IgA[Data_IgA$COVID19_dT_type=="Symptom Onset","COVID19_dT"],
         Data_IgA[Data_IgA$COVID19_dT_type=="Symptom Onset",i],col="#0000FF88")
}
#dev.off()
rm(i)
#


#####
#####Supplementary Figure 4 - ROC Analysis####


#ROC-Analysis for SARS CoV 2 Antigens

#R package Desctools AUC function, calculates AUC with a given FPR and TPR
AUC <- function(x, y, from=min(x, na.rm=TRUE), to = max(x, na.rm=TRUE), 
                method=c("trapezoid", "step", "spline", "linear"), 
                absolutearea = FALSE, subdivisions = 100, na.rm = FALSE, ...) {
  
  # calculates Area unter the curve
  # example:
  #   AUC( x=c(1,2,3,5), y=c(0,1,1,2))
  #   AUC( x=c(2,3,4,5), y=c(0,1,1,2))
  
  if(na.rm) {
    idx <- na.omit(cbind(x,y))
    x <- x[idx]
    y <- y[idx]
  }
  
  if (length(x) != length(y))
    stop("length x must equal length y")
  
  idx <- order(x)
  x <- x[idx]
  y <- y[idx]
  
  switch( match.arg( arg=method, choices=c("trapezoid","step","spline","linear") )
          , "trapezoid" = { a <- sum((apply( cbind(y[-length(y)], y[-1]), 1, mean))*(x[-1] - x[-length(x)])) }
          , "step" = { a <- sum( y[-length(y)] * (x[-1] - x[-length(x)])) }
          , "linear" = {
            a <- MESS_auc(x, y, from = from , to = to, type="linear", 
                          absolutearea=absolutearea, subdivisions=subdivisions, ...)
          }
          , "spline" = { 
            a <- MESS_auc(x, y, from = from , to = to, type="spline", 
                          absolutearea=absolutearea, subdivisions=subdivisions, ...)
            # a <- integrate(splinefun(x, y, method="natural"), lower=min(x), upper=max(x))$value 
          }
  )
  return(a)
}
MESS_auc <- function(x, y, from = min(x, na.rm=TRUE), to = max(x, na.rm=TRUE), type=c("linear", "spline"), 
                     absolutearea=FALSE, subdivisions =100, ...) {
  
  type <- match.arg(type)
  
  # Sanity checks
  stopifnot(length(x) == length(y))
  stopifnot(!is.na(from))
  
  if (length(unique(x)) < 2)
    return(NA)
  
  if (type=="linear") {
    ## Default option
    if (absolutearea==FALSE) {
      values <- approx(x, y, xout = sort(unique(c(from, to, x[x > from & x < to]))), ...)
      res <- 0.5 * sum(diff(values$x) * (values$y[-1] + values$y[-length(values$y)]))
    } else { ## Absolute areas
      ## This is done by adding artificial dummy points on the x axis
      o <- order(x)
      ox <- x[o]
      oy <- y[o]
      
      idx <- which(diff(oy >= 0)!=0)
      newx <- c(x, x[idx] - oy[idx]*(x[idx+1]-x[idx]) / (y[idx+1]-y[idx]))
      newy <- c(y, rep(0, length(idx)))
      values <- approx(newx, newy, xout = sort(unique(c(from, to, newx[newx > from & newx < to]))), ...)
      res <- 0.5 * sum(diff(values$x) * (abs(values$y[-1]) + abs(values$y[-length(values$y)])))
    }
    
  } else { ## If it is not a linear approximation
    if (absolutearea)
      myfunction <- function(z) { abs(splinefun(x, y, method="natural")(z)) }
    
    else
      myfunction <- splinefun(x, y, method="natural")
    
    res <- integrate(myfunction, lower=from, upper=to, subdivisions=subdivisions)$value
    
  }
  
  res
  
}




#ROC for SARS-CoV-2 Antigens IgG detection
AGs <- c("SARS2_Spike_Trimer","SARS2_RBD","SARS2_S1","SARS2_S2","SARS2_N","SARS2_N_NTD")
#List for export of performance
Export <- list()
#pdf("output/Suppl/SupplFig4_ROC_Analysis.pdf",9,6)
par(mfrow=c(2,3),mar=c(4,4,4,4))
for (k in c("IgG","IgA")){
  for (j in AGs){
    if (k=="IgG"){data <- Data_IgG[,c("COVID19_infection",j)]
    }else if (k=="IgA"){data <- Data_IgA[,c("COVID19_infection",j)]
    }
    data <- data[order(data[,j],decreasing=T),]
    TPR <- rep(0,nrow(data)+1)
    FPR <- rep(0,nrow(data)+1)
    for (i in 1:nrow(data)){#Loop to determine TPR and FPR in case the cutoff was above sample i
      tmp <- data[1:i,]
      TPR[i+1] <- nrow(tmp[tmp$COVID19_infection=="+",])/nrow(data[data$COVID19_infection=="+",])
      FPR[i+1] <- nrow(tmp[tmp$COVID19_infection=="-",])/nrow(data[data$COVID19_infection=="-",])
    }
    #Data frame with TPR and FPR for cutoffs set ABOVE the respective sample
    x <- cbind("Duplicated"=duplicated(data[,j]),data,"TPR"=TPR[1:(nrow(data))],"FPR"=FPR[1:(nrow(data))])
    #If two or more samples are tied in score, no cutoff can be set inbetween them, 
    #we must therefore remove these cutoff values from FPR and TPR and leave the cutoffs above and below the tied group
    #We exclude all "duplicated" MFI values, which will be the second/third/etc. in order from the top as given by the duplicated 
    #function and since the TPR and FPR are set to reflect cutoff ABOVE the respective sample, 
    #the lower ones of the group should also be the ones to get removed, sicn ethey reflect an impossible cutoff
    x <- x[x$Duplicated==F,2:5]
    Performance <- x[,2:4] 
    names(Performance) <- c("CO (> X MFI)","Sensitivity (%)","Specificity (%)") 
    Performance$`Specificity (%)` <- 1-Performance$`Specificity (%)` #FPR into Specificity
    Performance[,2:3] <- Performance[,2:3]*100 # In Percent
    #add a column with overall correctly classified samples, 
    #can be calculated via sens and spec and numbers of true positive and neagtive samples
    Performance$`Correct (%)` <- (Performance$`Sensitivity (%)`/100*nrow(data[data$COVID19_infection=="+",])+
                                    Performance$`Specificity (%)`/100*nrow(data[data$COVID19_infection=="-",]))/nrow(data)*100
    #plot ROC curve
    plot(1,1,cex=0,xlim=c(1,0),ylim=c(0,1),main=paste("ROC Analyis\n",j,k),
         xlab="1 - False Positive Rate",ylab="True positive Rate")
    lines(1-x$FPR,x$TPR)
    abline(1,-1)
    mtext(paste("AUC =",round(AUC(x$FPR,x$TPR),4)),side=1,at=0,adj=1,line=-2,cex=0.7)
    rm(x,i,FPR,TPR,data,tmp)
    #save performance in list for export
    Export[[paste(j,k)]] <- Performance
    rm(Performance)
  }
}
#dev.off()
#Export of List of possible cutoffs with respective peformance
for(i in AGs){
  #write.table(Export[[paste(i,"IgG")]],paste("output/Suppl/ROC_",i,"_IgG.csv",sep=""),sep=",",row.names=F,col.names=T)
  #write.table(Export[[paste(i,"IgA")]],paste("output/Suppl/ROC_",i,"_IgA.csv",sep=""),sep=",",row.names=F,col.names=T)
}






rm(AGs,i,j,k,Export,AUC,MESS_auc)

#####
#####Supplementary Figure 7 - Time series for hCoV antigens####


#Retain only patients with multiple samples for this panel
IgA <- Data_IgA[!is.na(Data_IgA$Time_series_PatientID),]
IgG <- Data_IgG[!is.na(Data_IgG$Time_series_PatientID),]

#Display time course of response for all 6 SARS2 antigens
AGs <- c("NL63_S1","NL63_N","NL63_N_NTD","229E_S1","229E_N","229E_N_NTD",
         "OC43_S1","OC43_N","OC43_N_NTD","HKU1_S1","HKU1_N","HKU1_N_NTD")
col <- colorRampPalette(brewer.pal(11,"Spectral")[c(8:11,1:3)])(5)
#pdf(paste("output/Suppl/SupplFig7_Time_Series_hCoVs.pdf",sep=""),6,6.8)
par(mfrow=c(3,2),mar=c(5,5,4,4))
for (j in AGs){
  plot(1,1,cex=0,xlim=range(as.numeric(IgG$COVID19_dT)),ylim=range(IgG[,j]),
       xlab="dT",ylab=paste(j),log="y",main=paste("IgG",j))
  c <- 0
  for (i in unique(IgG$Time_series_PatientID)){
    c <- c+1#counter for colors
    points(IgG[IgG$Time_series_PatientID==i,"COVID19_dT"],IgG[IgG$Time_series_PatientID==i,j],pch=1,cex=1.5,col=col[c])
    lines(IgG[IgG$Time_series_PatientID==i,"COVID19_dT"],IgG[IgG$Time_series_PatientID==i,j],lty=1,col=col[c])
  }
  plot(1,1,cex=0,xlim=range(as.numeric(IgA$COVID19_dT)),ylim=range(IgA[,j]),
       xlab="dT",ylab=paste(j),log="y",main=paste("IgA",j))
  c <- 0
  for (i in unique(IgG$Time_series_PatientID)){
    c <- c+1#counter for colors
    points(IgA[IgA$Time_series_PatientID==i,"COVID19_dT"],IgA[IgA$Time_series_PatientID==i,j],pch=1,cex=1.5,col=col[c])
    lines(IgA[IgA$Time_series_PatientID==i,"COVID19_dT"],IgA[IgA$Time_series_PatientID==i,j],lty=1,col=col[c])
  }
}
#dev.off()
rm(AGs,c,i,j,IgG,IgA,col)
#Reorder graphs and append legend manually in InkScape









#####
#####Extended Data Figure 2 - Assay QCs and Parallelism####


#Extended Data Figure 2, Panel b - Parallelism

#An experiment where multiple samples were diluted over multiple steps, including paired serum/plasma samples
#Parallelism is shown graphically only

#Read in of raw data
Para <- read.csv("input/Parallelism.csv.txt",h=T,sep=",",stringsAsFactors=F)
names(Para)[1] <- "Analyte"
#Matrix type and Donor ID are denoted in column title and are translated into colors and shapes for the plot
col <- c("#4DAF4A","#FF7F00","#FFFF33","#984EA3","#984EA3","#E41A1C","#E41A1C","#377EB8","#377EB8","#377EB8")
pchlist <- c(1,1,1,1,2,1,4,1,4,2)


#Plots MFI vs Dilution for all, Note that only Spike Trimer and RBD are shown in the paper
#pdf(paste("output/Suppl/SupplFig2_Parallelism.pdf",sep=""),12,6)
par(mfrow=c(1,2))
for (i in unique(Para$Analyte)){
  par(mfrow=c(1,2))
  plot(1,1,cex=0,ylim=range(c(Para[Para$Analyte==i&Para$det=="IgG",4:13])),xlim=c(100,12800),log="xy",
       ylab="MFI",xlab="Sample Dilution",xaxt="n",main=paste("IgG -",i))
  axis(1,at=c(100,200,400,800,1600,3200,6400,12800))
  c <- 0
  for(j in names(Para)[4:13]){
    c <- c+1
    points(Para[Para$Analyte==i&Para$det=="IgG","DF"],Para[Para$Analyte==i&Para$det=="IgG",j],col=col[c],pch=pchlist[c])
    lines(Para[Para$Analyte==i&Para$det=="IgG","DF"],Para[Para$Analyte==i&Para$det=="IgG",j],col=paste(col[c],"70",sep=""),lwd=2)
  }
  legend("topright",c("Serum","Plasma (EDTA)","Plasma (Heparin)"),pch=c(1,2,4),cex=0.8,bg="#00000000")
  
  plot(1,1,cex=0,ylim=range(c(Para[Para$Analyte==i&Para$det=="IgA",4:13])),xlim=c(100,12800),log="xy",
       ylab="MFI",xlab="Sample Dilution",xaxt="n",main=paste("IgA -",i))
  axis(1,at=c(100,200,400,800,1600,3200,6400,12800))
  c <- 0
  for(j in names(Para)[4:13]){
    c <- c+1
    points(Para[Para$Analyte==i&Para$det=="IgA","DF"],Para[Para$Analyte==i&Para$det=="IgA",j],col=col[c],pch=pchlist[c])
    lines(Para[Para$Analyte==i&Para$det=="IgA","DF"],Para[Para$Analyte==i&Para$det=="IgA",j],col=paste(col[c],"70",sep=""),lwd=2)
  }
  legend("topright",c("Serum","Plasma (EDTA)","Plasma (Heparin)"),pch=c(1,2,4),cex=0.8,bg="#00000000")
}
#dev.off()
rm(c,i,j,col,pchlist,Para)


#Extended Data Table 1 & Extended Data Figure 2, Panel a - QC performance

#QC sample performance over 17 plates (which includes measurement of all samples in this paper)
#Plate 14 Blank 1 was an outlier and was contaminated and will thus be excluded

#Read-in of QC Data
QC_IgA <- read.csv("input/QC_IgA.csv.txt",h=T,sep=",")
names(QC_IgA) <- gsub("X229","229",names(QC_IgA))#remove appended x from 229E

QC_IgG <- read.csv("input/QC_IgG.csv.txt",h=T,sep=",")
names(QC_IgG) <- gsub("X229","229",names(QC_IgG))#remove appended x from 229E


#Dataframe for export of QC Performance (Extended Data Table 1)
rep <- data.frame(c("IgG Mean(QC1)","IgG %CV(QC1)","IgG Mean(QC2)","IgG %CV(QC2)",
                    "IgG Mean(QC3)","IgG %CV(QC3)","IgG Mean(Blank)","IgG %CV(Blank)",
                    "IgA Mean(QC1)","IgA %CV(QC1)","IgA Mean(QC2)","IgA %CV(QC2)",
                    "IgA Mean(QC3)","IgA %CV(QC3)","IgA Mean(Blank)","IgA %CV(Blank)"))
QC_IgG <- QC_IgG[c(1:110,112:136),]#Row 111 is the outlier
QC_IgA <- QC_IgA[c(1:110,112:136),]#Row 111 is the outlier

#Display for the QCs of each individual antigen - RBD and Spike Trimer are shown in Ext. Data Fig 2
#Loop at the same time complies a summary table
nplates <- 17
#pdf(paste("output/Suppl/SupplFig2_QC_performance.pdf",sep=""),7.84,5.66)
for (i in names(QC_IgG)[3:ncol(QC_IgG)]){
  plot(1,1,log="y",xaxt="none",ylab="MFI",xlab="Plate Number",ylim=c(10,50000),xlim=c(1-0.2*(nplates-1),nplates),
       main=paste("QC performance -",i))  
  axis(1,1:nplates)
  #IgG
  #QC1
  tmp <- rep(0,nplates)
  for (j in 1:nplates){tmp[j]<-mean(QC_IgG[QC_IgG$Assay_Plate==j&QC_IgG$SampleID_Assay=="QC1",i])}
  lines(1:nplates,tmp,col="dodgerblue")
  points(QC_IgG[QC_IgG$SampleID_Assay=="QC1","Assay_Plate"],QC_IgG[QC_IgG$SampleID_Assay=="QC1",i],col="dodgerblue")
  m1 <- mean(QC_IgG[QC_IgG$SampleID_Assay=="QC1",i])
  cv1 <- round(sd(QC_IgG[QC_IgG$SampleID_Assay=="QC1",i])/m1*100,1)
  m1 <- round(m1,0)
  #QC2
  tmp <- rep(0,nplates)
  for (j in 1:nplates){tmp[j]<-mean(QC_IgG[QC_IgG$Assay_Plate==j&QC_IgG$SampleID_Assay=="QC2",i])}
  lines(1:nplates,tmp,col="purple")
  points(QC_IgG[QC_IgG$SampleID_Assay=="QC2","Assay_Plate"],QC_IgG[QC_IgG$SampleID_Assay=="QC2",i],col="purple")
  m2 <- mean(QC_IgG[QC_IgG$SampleID_Assay=="QC2",i])
  cv2 <- round(sd(QC_IgG[QC_IgG$SampleID_Assay=="QC2",i])/m2*100,1)
  m2 <- round(m2,0)
  #QC3
  tmp <- rep(0,nplates)
  for (j in 1:nplates){tmp[j]<-mean(QC_IgG[QC_IgG$Assay_Plate==j&QC_IgG$SampleID_Assay=="QC3",i])}
  lines(1:nplates,tmp,col="hotpink3")
  points(QC_IgG[QC_IgG$SampleID_Assay=="QC3","Assay_Plate"],QC_IgG[QC_IgG$SampleID_Assay=="QC3",i],col="hotpink3")
  m3 <- mean(QC_IgG[QC_IgG$SampleID_Assay=="QC3",i])
  cv3 <- round(sd(QC_IgG[QC_IgG$SampleID_Assay=="QC3",i])/m3*100,1)
  m3 <- round(m3,0)
  #Blank
  tmp <- rep(0,nplates)
  for (j in 1:nplates){tmp[j]<-mean(QC_IgG[QC_IgG$Assay_Plate==j&QC_IgG$SampleID_Assay=="Blank",i])}
  lines(1:nplates,tmp,col="gray20")
  points(QC_IgG[QC_IgG$SampleID_Assay=="Blank","Assay_Plate"],QC_IgG[QC_IgG$SampleID_Assay=="Blank",i],col="gray20")
  m4 <- mean(QC_IgG[QC_IgG$SampleID_Assay=="Blank",i])
  cv4 <- round(sd(QC_IgG[QC_IgG$SampleID_Assay=="Blank",i])/m4*100,1)
  m4 <- round(m4,0)
  #IgA
  #QC1
  tmp <- rep(0,nplates)
  for (j in 1:nplates){tmp[j]<-mean(QC_IgA[QC_IgA$Assay_Plate==j&QC_IgA$SampleID_Assay=="QC1",i])}
  lines(1:nplates,tmp,col="dodgerblue")
  points(QC_IgA[QC_IgA$SampleID_Assay=="QC1","Assay_Plate"],QC_IgA[QC_IgA$SampleID_Assay=="QC1",i],col="dodgerblue",pch=4)
  m5 <- mean(QC_IgA[QC_IgA$SampleID_Assay=="QC1",i])
  cv5 <- round(sd(QC_IgA[QC_IgA$SampleID_Assay=="QC1",i])/m5*100,1)
  m5 <- round(m5,0)
  #QC2
  tmp <- rep(0,nplates)
  for (j in 1:nplates){tmp[j]<-mean(QC_IgA[QC_IgA$Assay_Plate==j&QC_IgA$SampleID_Assay=="QC2",i])}
  lines(1:nplates,tmp,col="purple")
  points(QC_IgA[QC_IgA$SampleID_Assay=="QC2","Assay_Plate"],QC_IgA[QC_IgA$SampleID_Assay=="QC2",i],col="purple",pch=4)
  m6 <- mean(QC_IgA[QC_IgA$SampleID_Assay=="QC2",i])
  cv6 <- round(sd(QC_IgA[QC_IgA$SampleID_Assay=="QC2",i])/m6*100,1)
  m6 <- round(m6,0)
  #QC3
  tmp <- rep(0,nplates)
  for (j in 1:nplates){tmp[j]<-mean(QC_IgA[QC_IgA$Assay_Plate==j&QC_IgA$SampleID_Assay=="QC3",i])}
  lines(1:nplates,tmp,col="hotpink3")
  points(QC_IgA[QC_IgA$SampleID_Assay=="QC3","Assay_Plate"],QC_IgA[QC_IgA$SampleID_Assay=="QC3",i],col="hotpink3",pch=4)
  m7 <- mean(QC_IgA[QC_IgA$SampleID_Assay=="QC3",i])
  cv7 <- round(sd(QC_IgA[QC_IgA$SampleID_Assay=="QC3",i])/m7*100,1)
  m7 <- round(m7,0)
  #Blank
  tmp <- rep(0,nplates)
  for (j in 1:nplates){tmp[j]<-mean(QC_IgA[QC_IgA$Assay_Plate==j&QC_IgA$SampleID_Assay=="Blank",i])}
  lines(1:nplates,tmp,col="gray20")
  points(QC_IgA[QC_IgA$SampleID_Assay=="Blank","Assay_Plate"],QC_IgA[QC_IgA$SampleID_Assay=="Blank",i],col="gray20",pch=4)
  m8 <- mean(QC_IgA[QC_IgA$SampleID_Assay=="Blank",i])
  cv8 <- round(sd(QC_IgA[QC_IgA$SampleID_Assay=="Blank",i])/m8*100,1)
  m8 <- round(m8,0)
  #Legend
  legend((1-0.23*(nplates-1)),15000,"IgG",bty="n")
  legend((1-0.255*(nplates-1)),5500,legend=c("","","",""),
         col=c("dodgerblue","purple","hotpink3","gray20"),bty="n",cex=0.7,pch=1)
  legend((1-0.265*(nplates-1)),5500,c("QC1","QC2","QC3","Blank"),bty="n",cex=0.7)
  legend((1-0.15*(nplates-1)),8000,"Mean",bty="n",cex=0.7,adj=1)
  legend((1-0.15*(nplates-1)),5500,round(c(m1,m2,m3,m4),0),bty="n",cex=0.7,adj = 1)
  legend((1-0.09*(nplates-1)),8000,"%CV",bty="n",cex=0.7,adj=1)
  legend((1-0.09*(nplates-1)),5500,round(c(cv1,cv2,cv3,cv4),1),bty="n",cex=0.7,adj = 1)
  legend((1-0.23*(nplates-1)),750,"IgA",bty="n")
  legend((1-0.255*(nplates-1)),275,legend=c("","","",""),
         col=c("dodgerblue","purple","hotpink3","gray20"),bty="n",cex=0.7,pch=4)
  legend((1-0.265*(nplates-1)),275,c("QC1","QC2","QC3","Blank"),bty="n",cex=0.7)
  legend((1-0.15*(nplates-1)),400,"Mean",bty="n",cex=0.7,adj=1)
  legend((1-0.15*(nplates-1)),275,round(c(m5,m6,m7,m8),0),bty="n",cex=0.7,adj = 1)
  legend((1-0.09*(nplates-1)),400,"%CV",bty="n",cex=0.7,adj=1)
  legend((1-0.09*(nplates-1)),275,round(c(cv5,cv6,cv7,cv8),1),bty="n",cex=0.7,adj = 1)
  tmp <- c(m1,cv1,m2,cv2,m3,cv3,m4,cv4,m5,cv5,m6,cv6,m7,cv7,m8,cv8)
  rep <- cbind(rep,tmp)
}
#dev.off()

#Export stats as .csv file, these values are the measure for inter assay variance and are also displayed in Supplementary Table 1
rownames(rep) <- rep[,1]
rep <- rep[,2:ncol(rep)]
names(rep) <- names(QC_IgG)[3:ncol(QC_IgG)]
rm(i,j,tmp,cv1,cv2,cv3,cv4,m1,m2,m3,m4,cv5,cv6,cv7,cv8,m5,m6,m7,m8,nplates)
#write.table(rep,"output/Suppl/SupplFig2_QC_Performance.csv",row.names=T,col.names=NA,sep=",")
rm(QC_IgA,QC_IgG,rep)



#####
