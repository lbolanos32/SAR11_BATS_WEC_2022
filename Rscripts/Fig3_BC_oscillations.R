#FIGURE BRAY-Curtis distances#

In the oscillations plot, we normalized all the ASVs using Deseq so we can take into account ASVs that may be minimal and present in deep coverage samples. Evaluate whether these ASVs are appearing and dissapeing from early years in the survey compared to the latest samples.

###################################################
#####BRAY CURTIS dissimilarities vs TIME ##########
###################################################

############
####BATS####
############
library("phyloseq")
library("ggplot2")
library("dplyr")
library("tidyr")
library("reshape")
library("see")
library("cowplot")
library("grid")
library("gridExtra")
library("DESeq2")
library("ggpubr")

#####Generate the phyloseq object BATScf #####

#Create the phyloseq object

count_tabS11phyTS5m <- read.table("BIOSotu5.otu", header=T, row.names=1, check.names=F)
sample_info_tabS11phyTS5m <- read.table("BIOSenvts5.env", header=T, row.names=1, check.names=F, sep ="\t")
tax_tabS11phyTS5m <- as.matrix(read.table("BIOStax5.tax", header=T, row.names=1, check.names=F, sep="\t",na.strings = "#NA"))

#Use as date
sample_info_tabS11phyTS5m$Date <- as.Date(with(sample_info_tabS11phyTS5m, paste(year, month, day,sep="-")), "%Y-%m-%d")

S11OTUphyTS5m = otu_table(count_tabS11phyTS5m, taxa_are_rows = TRUE)
S11TAXphyTS5m = tax_table(tax_tabS11phyTS5m)
S11SAMphyTS5m= sample_data(sample_info_tabS11phyTS5m)

phyTS5m<-phyloseq(S11OTUphyTS5m,S11TAXphyTS5m,S11SAMphyTS5m)

phyTS5mV1 <- prune_samples(sample_sums(phyTS5m)>=1000, phyTS5m) #prune samples less than 1000 reads (eliminate 17 --> 146/163)

phyTS5mpruneV2<-filter_taxa(phyTS5mV1, function(x) sum(x >= 3) > (0.02 *length(x)), TRUE) #V1=916/2875 Again from phyTS5mprune, V2 (146 samples, .01) ~500 ASVs

####################SAR11 ECOTYPES###################

SAR11_BATSphy = subset_taxa(phyTS5mpruneV2, Order=="SAR11_clade")

####MISEQ###

SAR11_BATSphyMiSeq<-subset_samples(SAR11_BATSphy, tech=="MiSeq")

###Rarefy 7200 ######
set.seed(717)
phyeuphminrarWSAR11_BATSphyMiSeq = rarefy_even_depth(SAR11_BATSphyMiSeq, sample.size = 7200)

#Generate the Bray Curtis distances of the samples
bray_not_na_BATSSAR11MiSeq <- phyloseq::distance(physeq = phyeuphminrarWSAR11_BATSphyMiSeq, method = "bray")

library(reshape2)
df_BC <- melt(as.matrix(bray_not_na_BATSSAR11MiSeq), varnames = c("row_BC", "col_BC"))
names(df_BC)[3] <- "Bray_Curtis"
rownames(df_BC)<-paste(df_BC$row, df_BC$col, sep="_")

#Generate the days distances between samples to match bray Curtis
DFminrarWSAR11_BATSphyMiSeq<-data.frame(sample_data(phyeuphminrarWSAR11_BATSphyMiSeq))
DFminrarWSAR11_BATSphyMiSeq$rnames<-rownames(DFminrarWSAR11_BATSphyMiSeq)

#Create a Matrix of days between samples
days_s2010 <- difftime(DFminrarWSAR11_BATSphyMiSeq$Date,as.Date("2016-03-01","%Y-%m-%d"))
dist_days <- as.matrix(dist(days_s2010,diag=TRUE,upper=TRUE))
rownames(dist_days) <- DFminrarWSAR11_BATSphyMiSeq$rnames; colnames(dist_days) <- DFminrarWSAR11_BATSphyMiSeq$rnames

dist_days[upper.tri(dist_days)] <- NA
dist_days <- melt(as.matrix(dist_days), varnames = c("row_day", "col_day"))
names(dist_days)[3] <- "days"
rownames(dist_days)<-paste(dist_days$row, dist_days$col, sep="_")

#merge the 2 generated dfs
days_BC<-merge(dist_days, df_BC, by=0, all=TRUE)
days_BC_lowt<-na.omit(days_BC) #Remove NAs
days_BC_lowt_final<-days_BC_lowt[,c(2:7)]
days_BC_lowt_final_no0<-days_BC_lowt_final[days_BC_lowt_final$days != 0, ]
#add a column with consecutive numbers so we can label 
days_BC_lowt_final_no0$comparison <- 1:nrow(days_BC_lowt_final_no0)


####FLEX###
SAR11_BATSphy454<-subset_samples(SAR11_BATSphy, tech=="FLX454")

###Rarefy 465 ######
set.seed(717)
phyeuphminrarWSAR11_BATSphy454 = rarefy_even_depth(SAR11_BATSphy454, sample.size = 465) 

bray_not_na_BATSSAR11_454 <- phyloseq::distance(physeq = phyeuphminrarWSAR11_BATSphy454, method = "bray")
df_BC454 <- melt(as.matrix(bray_not_na_BATSSAR11_454), varnames = c("row_BC", "col_BC"))
names(df_BC454)[3] <- "Bray_Curtis"
rownames(df_BC454)<-paste(df_BC454$row, df_BC454$col, sep="_")

#Generate the days distances between samples to match bray Curtis
DFminrarWSAR11_BATSphy454<-data.frame(sample_data(phyeuphminrarWSAR11_BATSphy454))
DFminrarWSAR11_BATSphy454$rnames<-rownames(DFminrarWSAR11_BATSphy454)

days_s2010_454 <- difftime(DFminrarWSAR11_BATSphy454$Date,as.Date("1990-03-01","%Y-%m-%d"))
dist_days454 <- as.matrix(dist(days_s2010_454,diag=TRUE,upper=TRUE))
rownames(dist_days454) <- DFminrarWSAR11_BATSphy454$rnames; colnames(dist_days454) <- DFminrarWSAR11_BATSphy454$rnames

dist_days454[upper.tri(dist_days454)] <- NA
days_s2010_454 <- melt(as.matrix(dist_days454), varnames = c("row_day", "col_day"))
names(days_s2010_454)[3] <- "days"
rownames(days_s2010_454)<-paste(days_s2010_454$row, days_s2010_454$col, sep="_")

#merge the 2 generated dfs
days_BC454<-merge(days_s2010_454, df_BC454, by=0, all=TRUE)
days_BC454_lowt<-na.omit(days_BC454) #Remove NAs
days_BC454_lowt_final<-days_BC454_lowt[,c(2:7)]
days_BC454_lowt_final_no0<-days_BC454_lowt_final[days_BC454_lowt_final$days != 0, ]
#add a column with consecutive numbers so we can label 
days_BC454_lowt_final_no0$comparison <- 1:nrow(days_BC454_lowt_final_no0)

##PLOT de MiSEQ rarefy 465 #Max 4551.042 days
BC_454_Bats<-ggscatter(days_BC454_lowt_final_no0, x = "days", y = "Bray_Curtis",color = "black", shape = 21, size = 2)+ scale_x_continuous(breaks=seq(0, 4600, by = 365.25),limits=(c(0,4600)),labels=c("0","1","2","3","4","5","6","7","8","9","10","11","12"))+xlab("years")+ylab("Bray-Curtis dissimilarity")+ylim(c(0,1))

#############
#####WEC#####
#############


countwec <- read.table("WEC_mar.otu", header=T, row.names=1, check.names=F)
sampleinfowec <- read.table("WEC_mar.env", header=T, row.names=1, check.names=F, sep ="\t")
taxwec <- as.matrix(read.table("WEC_marV2.tax", header=T, row.names=1, check.names=F, na.strings="", sep="\t")) 

OTUmar = otu_table(countwec, taxa_are_rows = TRUE)
TAXmar = tax_table(taxwec)
SAMar= sample_data(sampleinfowec)

WEC<-phyloseq(OTUmar,TAXmar,SAMar)
WECcf= filter_taxa(WEC, function(x) sum(x > 1) > (0.0015*length(x)), TRUE)

#####WEC_SAR11####

SAR11_WECphy = subset_taxa(WECcf, Order=="SAR11_clade") #min 15 and max 32147

####Using rarefied dataset. 

set.seed(717)
phyeuphminrarWECcfS11 = rarefy_even_depth(SAR11_WECphy, sample.size = 1000)

bray_WECS11 <- phyloseq::distance(physeq = phyeuphminrarWECcfS11, method = "bray")

library(reshape2)
df_WECfS11 <- melt(as.matrix(bray_WECS11), varnames = c("row_Wf", "col_Wf"))
names(df_WECfS11)[3] <- "Bray_Curtis"
rownames(df_WECfS11)<-paste(df_WECfS11$row_Wf, df_WECfS11$col_Wf, sep="_") #df_WECfS11 has the BRAY-CURTIS distances 

#####NOW GENERATE THE Days distance
DFminrarWECS11<-data.frame(sample_data(phyeuphminrarWECcfS11))
DFminrarWECS11$rnames<-rownames(DFminrarWECS11)

DFminrarWECS11$Date<-as.Date(DFminrarWECS11$Date,"%m/%d/%Y") 

#Create a Matrix of days between samples
days_s2010 <- difftime(DFminrarWECS11$Date,as.Date("2010-01-01","%Y-%m-%d"))
dist_daysWECS11 <- as.matrix(dist(days_s2010,diag=TRUE,upper=TRUE))

rownames(dist_daysWECS11) <- DFminrarWECS11$rnames; colnames(dist_daysWECS11) <- DFminrarWECS11$rnames

dist_daysWECS11[upper.tri(dist_daysWECS11)] <- NA
dist_daysWECS11 <- melt(as.matrix(dist_daysWECS11), varnames = c("row_day", "col_day"))
names(dist_daysWECS11)[3] <- "days"
rownames(dist_daysWECS11)<-paste(dist_daysWECS11$row, dist_daysWECS11$col, sep="_")

###MERGE BC and days
days_BCWECS11<-merge(dist_daysWECS11, df_WECfS11, by=0, all=TRUE)
days_BCWECS11_lowt<-na.omit(days_BCWECS11) #Remove NAs
days_BCWECS11_lowt_final<-days_BCWECS11_lowt[,c(2:7)]
days_BCWECS11_lowt_final_no0<-days_BCWECS11_lowt_final[days_BCWECS11_lowt_final$days != 0, ] #DIM! 34191,6
#add a column with consecutive numbers so we can label 
days_BCWECS11_lowt_final_no0$comparison <- 1:nrow(days_BCWECS11_lowt_final_no0)

#PLOT WEC 2506 days 
WEC_BC<-ggscatter(days_BCWECS11_lowt_final_no0, x = "days", y = "Bray_Curtis",color = "black", size = .05)+ scale_x_continuous(breaks=seq(0, 4600, by = 365.25),limits=(c(0,4600)),labels=c("0","1","2","3","4","5","6","7","8","9","10","11","12"))+xlab("time between samples (years)")+ylab("Bray-Curtis dissimilarity")+ylim(c(0,1))+theme(panel.background = element_blank(),legend.position = "none",strip.background = element_blank(),axis.text.y=element_text(colour = 'black',size=22), axis.text.x=element_text(colour = 'black',size=26),axis.line = element_line(colour = "black"),axis.title.x=element_text(colour = 'black',size=26),axis.title.y=element_text(colour = 'black',size=20) )

##PLOT de MiSEQ rarefy 7200 #MAX in days = 1012
BC_MiSeq_Bats<-ggscatter(days_BC_lowt_final_no0, x = "days", y = "Bray_Curtis",color = "black", shape = 21, size = .8)+ scale_x_continuous(breaks=seq(0, 4600, by = 365.25),limits=(c(0,4600)),labels=c("0","1","2","3","4","5","6","7","8","9","10","11","12"))+xlab("time between samples (years)")+ylab("Bray-Curtis dissimilarity")+ylim(c(0,1))+theme(panel.background = element_blank(),legend.position = "none",strip.background = element_blank(),axis.text.y=element_text(colour = 'black',size=22), axis.text.x=element_text(colour = 'black',size=26),axis.line = element_line(colour = "black"),axis.title.x=element_text(colour = 'black',size=26),axis.title.y=element_text(colour = 'black',size=20) )

##PLOT de MiSEQ rarefy 465 #Max 4551.042 days
BC_454_Bats<-ggscatter(days_BC454_lowt_final_no0, x = "days", y = "Bray_Curtis",color = "black", shape = 21, size = .8)+ scale_x_continuous(breaks=seq(0, 4600, by = 365.25),limits=(c(0,4600)),labels=c("0","1","2","3","4","5","6","7","8","9","10","11","12"))+xlab("time between samples (years)")+ylab("Bray-Curtis dissimilarity")+ylim(c(0,1))+theme(panel.background = element_blank(),legend.position = "none",strip.background = element_blank(),axis.text.y=element_text(colour = 'black',size=22), axis.text.x=element_text(colour = 'black',size=26),axis.line = element_line(colour = "black"),axis.title.x=element_text(colour = 'black',size=26),axis.title.y=element_text(colour = 'black',size=20) )


plot_grid(BC_454_Bats,BC_MiSeq_Bats,WEC_BC, ncol=1, align="hv", labels = c('a', 'b', 'c'),label_size = 24)

