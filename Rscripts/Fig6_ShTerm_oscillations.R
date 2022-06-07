#FIGURE 6 TO ANALYZE the Influence of short-time environmental variability on Bray-Curtis SAR11 dissimilarities in the WCO.

#In the below analyses, environmental data were previously standardized to values between 0 and 1, based on the minimum and maximum values of each variable, using the formula: X'=x-min(x)/max(x)-min(x) -> [normalized = (x-min(x))/(max(x)-min(x))]

###Strategy###

#A) Generate the Bray Curtis Matrix -> two columns 
	#Grab the dates 
	#Generate a list of comparisons (n+1). Extract this list of comporisons from the two column. 

#############
#####WEC#####
#############

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
library(zoo)
library("ggdendro")

#The environmental file was changed for wavelet analysis so we can have exact 7 day distances  A method to reveal the periodic properties of a time series, such as wavelet transformation, will only yield valid results if observations are made at equidistant points in time. What can be accepted as “equidistant” depends on the particular application. For example, daily observations of a stock index are usually considered as being made at equidistant epochs, even though no observation is available at weekends when the stock exchange is closed.

#To overcome the lack of periodicity between the sampling dates, we converted our irregularly distributed observations into a fixed interval dataset (following Carey et al., 2016) by reassigning sampling days to the closest regularly spaced day on a 3-days interval throughout the time-series, and linearly interpolating the missing data for the rare occasions on which sampling did not occur.

setwd("/Users/luisbolanos/Documents/MyDrafts/InProgress/SAR11TS/Analysis/WEC")
countwec <- read.table("WEC_mar.otu", header=T, row.names=1, check.names=F)
sampleinfowec <- read.table("WEC_mar.env", header=T, row.names=1, check.names=F, sep ="\t")
taxwec <- as.matrix(read.table("WEC_marV2.tax", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

OTUmar = otu_table(countwec, taxa_are_rows = TRUE)
TAXmar = tax_table(taxwec)
SAMar= sample_data(sampleinfowec)

WEC<-phyloseq(OTUmar,TAXmar,SAMar)
WECcf= filter_taxa(WEC, function(x) sum(x > 1) > (0.0015*length(x)), TRUE)

#####WEC_SAR11####

SAR11_WECphy = subset_taxa(WECcf, Order=="SAR11_clade")

####Using rarefied dataset. 

set.seed(717)
phyeuphminrarWECcfS11 = rarefy_even_depth(SAR11_WECphy, sample.size = 1000)

bray_WECS11 <- phyloseq::distance(physeq = phyeuphminrarWECcfS11, method = "bray")

library(reshape2)
df_WECfS11 <- melt(as.matrix(bray_WECS11), varnames = c("row_Wf", "col_Wf"))
names(df_WECfS11)[3] <- "Bray_Curtis"
rownames(df_WECfS11)<-paste(df_WECfS11$row_Wf, df_WECfS11$col_Wf, sep="_") #df_WECfS11 has the BRAY-CURTIS distance
df_WECfS11$distname<-rownames(df_WECfS11)

#####NOW GENERATE THE distance between the days

DFminrarWECS11<-data.frame(sample_data(phyeuphminrarWECcfS11))
DFminrarWECS11$rnames<-rownames(DFminrarWECS11)

DFminrarWECS11$Date<-as.Date(DFminrarWECS11$Date,"%m/%d/%Y") 

#Sort by Date
DFminrarWECS11_sorted<-(DFminrarWECS11[order(DFminrarWECS11$Date),])

#Create a two column of the consecutive rnames
for(i in 1:nrow(DFminrarWECS11_sorted)) {
  DFminrarWECS11_sorted$rnames2[i]<- DFminrarWECS11_sorted$rnames[i+1]}

DFminrarWECS11_sorted$distname<-paste(DFminrarWECS11_sorted$rnames, DFminrarWECS11_sorted$rnames2, sep="_")

#Create a Matrix of days between samples
days_s2010 <- difftime(DFminrarWECS11$Date,as.Date("2010-01-01","%Y-%m-%d"))
dist_daysWECS11 <- as.matrix(dist(days_s2010,diag=TRUE,upper=TRUE))

rownames(dist_daysWECS11) <- DFminrarWECS11$rnames; colnames(dist_daysWECS11) <- DFminrarWECS11$rnames

dist_daysWECS11 <- melt(as.matrix(dist_daysWECS11), varnames = c("row_day", "col_day"))
names(dist_daysWECS11)[3] <- "days"
rownames(dist_daysWECS11)<-paste(dist_daysWECS11$row, dist_daysWECS11$col, sep="_")

###MERGE BC and days
days_BCWECS11<-merge(dist_daysWECS11, df_WECfS11, by=0, all=TRUE)
days_BCWECS11_lowt<-na.omit(days_BCWECS11) #Remove NAs
days_BCWECS11_lowt_final<-days_BCWECS11_lowt[,c(2:8)]
days_BCWECS11_lowt_final_no0<-days_BCWECS11_lowt_final[days_BCWECS11_lowt_final$days != 0, ] #DIM! 34191,6

##I have in the original data frame DFminrarWECS11_sorted the consecutive distances I want.

#NOW I need to extract the consecutive bray Curtis based on DFminrarWECS11_sorted$distname

BC_consecutive<-days_BCWECS11_lowt_final_no0 %>% filter(distname %in% DFminrarWECS11_sorted$distname)

#BC_consecutive needs to be pasted in the new tale WEC_mar.env

#####Get the only same exact periodicity and extrapolating the values 

#(a) GET ALL Mondays from 2/6/12 to 12/17/18, (b) set a DF with the equidistant weekly Time Frame, (c) put the samples names to the correspondent date or otherwise just N1:Nn to those whithout sample names, (c) do a subsequent distant consecutive numbers and (d) linear interpolation of all this data frame. 

#Function Mondays per year 
getAllMondays <- function(year) {
    days <- as.POSIXlt(paste(year, 1:366, sep="-"), format="%Y-%j")
    Ms <- days[days$wday==1]
    Ms[!is.na(Ms)]  # Needed to remove NA from day 366 in non-leap years
}

getAllMondays(2012)

Mday12<-data.frame(date_wv = getAllMondays(2012))
Mday13<-data.frame(date_wv = getAllMondays(2013))
Mday14<-data.frame(date_wv = getAllMondays(2014))
Mday15<-data.frame(date_wv = getAllMondays(2015))
Mday16<-data.frame(date_wv = getAllMondays(2016))
Mday17<-data.frame(date_wv = getAllMondays(2017))
Mday18<-data.frame(date_wv = getAllMondays(2018))

rbind(Mday12,Mday13,Mday14,Mday15,Mday16,Mday17,Mday18)


#Do the same as pasting the B-C distances to the 7-day table

#Once I have the distances I evaluate a linear interpolation of missing values for all of them. 

#At the end just cut till KT16S452 9/3/18
#Sort by Days
BC_consecutive[order(BC_consecutive$days),]

#Before merging remove the consecutive distances BC longer than 10 days. Although these are consecutive in the sampling, they are not strictly weekly
BC_consecutive_7day<- subset(BC_consecutive, BC_consecutive[ ,"days"] < 12)  #dim 186 distances from 350... 

####READ THE NEW 7 day table 

wl_7day_df <- read.table("/Users/luisbolanos/Documents/MyDrafts/InProgress/SAR11TS/Analysis/WEC/WEC_mar_wavelet.env.tab", header=T, row.names=1, check.names=F, sep ="\t") #This is sorted from the file I did

wl_7day_df_short<-wl_7day_df[ , c("ID","Station","Type","Date", "Date_WaveLet","season","DateEnv","NITRITE_0","NITRATE.NIT_0","AMMONIA_0","SILICATE_0","PHOSPHATE_0","Density_kg.m3_2.5","Par_uE.m2.s_2.5","salinity_PSU_2.5","oxygen_uM_2.5","temp_C_2.5","Fluorescence_volts_2.5","Transmission_._2.5","Chl_0m","maximum.wind.speed","average.weekly.pressure","total.rainfall.for.week","average.weekly.temperature","average.weekly.wind.speed","average.weekly.wave.height","average.weekly.river.flow")]

#####We need to inspect
wl_7day_df_short_log<-cbind(wl_7day_df_short[,c(1:7)],!is.na(wl_7day_df_short[,c(8:27)]))
wl_7day_df_short_log$rnames<-rownames(wl_7day_df_short_log)

#Get the samples after rarefying from DFminrarWECS11. This should be added as amplicon data available in the table. If rnames of DFminrarWECS11 are present in wl_7day_df add a column named "amplicons" with a value TRUE if not value should be FALSE

wl_7day_df_short_log$amplicon<-wl_7day_df_short_log$rnames %in% sort(DFminrarWECS11$rnames) #compare the two columns of rnames from different DFs
logical_7day<-wl_7day_df_short_log[c(5,8:27,29)]

logical_melt<-melt(logical_7day,id.vars = c("Date_WaveLet")) #Now ready for heat map

logical_melt$Date_WaveLet <- as.Date(logical_melt$Date_WaveLet,"%m/%d/%Y")
   
md_gral<-ggplot(data = logical_melt, aes(x = Date_WaveLet, y=variable,fill=value)) + geom_tile(colour="white",size=0.25)+scale_x_date( breaks = "1 year",minor_breaks = "1 week", date_labels = "%Y")+scale_fill_manual(values=c("TRUE"="black","FALSE"="white"),labels = c("NA", "Av"))+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=12,colour = "black"), axis.text.x=element_text(size=14,colour = "black"),legend.text=element_text(size=14))+geom_vline(xintercept=as.Date(c("2015-04-19","2017-04-25")),linetype=1, colour="firebrick3",size=2)

### Limit to the analysis we want to do and removing variables "Fluorescence_volts_2.5" and "Transmission_._2.5"
logical_melt_clean<-logical_melt[!(logical_melt$variable == "Fluorescence_volts_2.5"| logical_melt$variable == "Transmission_._2.5"),]

md_constr<-ggplot(data = logical_melt_clean, aes(x = Date_WaveLet, y=variable,fill=value)) + geom_tile(colour="white",size=0.25)+
scale_x_date( breaks = "1 month",minor_breaks = "1 week", date_labels = "%m/%y",limits = as.Date(c("2015-04-19","2017-04-15")))+scale_fill_manual(values=c("TRUE"="black","FALSE"="white"),labels = c("NA", "Av"))+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=12,colour = "black"), axis.text.x=element_text(size=7,colour = "black"),legend.text=element_text(size=12))

#Plot the two panels

plot_grid(md_gral,md_constr,ncol=1,labels = c('a', 'b'),label_size = 20)


############################
#Generate the distance names
############################
#Start from wl_7day_df_short
 #for interpolation

wl_7day_df_short$rnames<-rownames(wl_7day_df_short)

###Generate the normalize measurements

#Function to normalize all the environmental columns
A <- function(x) (x - min(x,na.rm = TRUE))/(max(x,na.rm = TRUE) - min(x,na.rm = TRUE))

wl_7day_df_short_norm<-data.frame(wl_7day_df_short[1:7], lapply(wl_7day_df_short[8:27], A), wl_7day_df_short[28])

n<-dim(wl_7day_df_short_norm)[1]
wl_7day_df_short_norm<-wl_7day_df_short_norm[1:(n-2),] #remove the last two rows that belong to N96 and N97. This makes problems when linearly interpolating

wl_7day_df_short_norm_interp<- data.frame(wl_7day_df_short_norm[1:7], na.approx(wl_7day_df_short_norm[8:27]), wl_7day_df_short_norm[28]) #Inerpolate linearly the missing values

#Ammonia has a blank in the first row and that sucks (NA), add the average of the whole vector in that position
addNH5_firstpos<-mean(wl_7day_df_short_norm$AMMONIA_0, na.rm=T)
wl_7day_df_short_norm_interp[1,10]<-addNH5_firstpos


#####EUCLIDEAN_ALL####
wl_7day_df_short_norm_interp$EucDist.All<-sqrt(c(0,diff(wl_7day_df_short_norm_interp$NITRITE_0)^2) + c(0,diff(wl_7day_df_short_norm_interp$NITRATE.NIT_0)^2)+ c(0,diff(wl_7day_df_short_norm_interp$AMMONIA_0)^2)+ c(0,diff(wl_7day_df_short_norm_interp$SILICATE_0)^2)+ c(0,diff(wl_7day_df_short_norm_interp$PHOSPHATE_0)^2)+ c(0,diff(wl_7day_df_short_norm_interp$Par_uE.m2.s_2.5)^2)+ c(0,diff(wl_7day_df_short_norm_interp$salinity_PSU_2.5)^2)+ c(0,diff(wl_7day_df_short_norm_interp$oxygen_uM_2.5)^2)+ c(0,diff(wl_7day_df_short_norm_interp$temp_C_2.5)^2)+ c(0,diff(wl_7day_df_short_norm_interp$Chl_0m)^2)+ c(0,diff(wl_7day_df_short_norm_interp$average.weekly.pressure)^2)+ c(0,diff(wl_7day_df_short_norm_interp$total.rainfall.for.week)^2)+ c(0,diff(wl_7day_df_short_norm_interp$average.weekly.wind.speed)^2)+ c(0,diff(wl_7day_df_short_norm_interp$average.weekly.wave.height)^2)+ c(0,diff(wl_7day_df_short_norm_interp$average.weekly.river.flow)^2))


#Function to get the distances between consecutive days by each independent variable
B<-function(x) abs(c(0,diff(x))) #If I do it in this way, distances will start in position 2 and finish in position 'N' 

wl_7day_df_short_norm_difs<-data.frame(wl_7day_df_short_norm_interp[1:7], lapply(wl_7day_df_short_norm_interp[8:27], B), wl_7day_df_short_norm_interp[28:29])


#names to match with Bray_Curtis distances
wl_7day_df_short_norm_difs$lag.value <- c(NA, wl_7day_df_short_norm_difs$rnames[-nrow(wl_7day_df_short_norm_difs)])

wl_7day_df_short_norm_difs$distnameR<-paste(wl_7day_df_short_norm_difs$rnames, wl_7day_df_short_norm_difs$lag.value, sep="_")
wl_7day_df_short_norm_difs$distnamebckwd<-paste(wl_7day_df_short_norm_difs$lag.value, wl_7day_df_short_norm_difs$rnames, sep="_")

#Merge distances and Euclidean to Bray Curtis                            
zz <- merge(wl_7day_df_short_norm_difs, BC_consecutive_7day, by.x='distnamebckwd', by.y='distname', all = TRUE)    #359 y #186


#write.table(zz,"/Users/luisbolanos/Documents/MisDrafts/InProgress/SAR11TS/Analysis/WEC/WEC_mar_wavelet1.env.tab", quote=FALSE, sep = "\t")
#After sorting the data by date. We have a period of time from 4/20/15 to 4/24/17 that we would only need to interpolate 15/105 weeks. So we have 90 real data points     

#merge also the ecotype relative abundances per sample. These are already normalized
zz$Date_WaveLet<-as.Date(zz$Date_WaveLet,"%m/%d/%Y") 
wl_7day_BC_sorted<-zz[order(zz$Date_WaveLet),]

wl_7day_BC_sorted$Bray_Curtis.apx<-c("NA","NA","NA",na.approx(wl_7day_BC_sorted$Bray_Curtis))

###wl_7day_BC_sorted#### THIS FILE NEEDS to be merged to the distances between relative abundances of ecotypes. 

####Relative abundance by ecotype

head(otu_table(phyeuphminrarWECcfS11))

glomV1filt<-tax_glom(phyeuphminrarWECcfS11, taxrank="Genus")
glomV1filtrel<-transform_sample_counts(glomV1filt, function(x){x / sum(x)})

SAR11tax<-as.data.frame(tax_table(glomV1filtrel)[,6])

SAR11_ASV<-as.data.frame(otu_table(glomV1filtrel))SAR11_ASV[ "Taxa" ] <- SAR11tax[,1]#dim(SAR11_ASV)#[1]  15 243ASV_frame2 <- SAR11_ASV[,-243]rownames(ASV_frame2) <- SAR11_ASV[,243]ASV_frw<-t(ASV_frame2) #Only select: IIa.B,I;Ia,Ia.1,I;Ia,Ia.3 

Ecotypes<-data.frame(ASV_frw[,c(1:9)])

md_to_add1a<-as.data.frame(sample_data(phyeuphminrarWECcfS11))[,4]

Ecotypes_date<-cbind(Ecotypes,md_to_add1a)

#names to match distances
Ecotypes_date$rnames<-rownames(Ecotypes_date)

#names to match distances
Ecotypes_date$lag.value <- c(NA, Ecotypes_date$rnames[-nrow(Ecotypes)])

Ecotypes_date$distnameR<-paste(Ecotypes_date$rnames, Ecotypes_date$lag.value, sep="_")
Ecotypes_date$distnameEco<-paste(Ecotypes_date$lag.value, Ecotypes_date$rnames, sep="_")

changes_relContr_ecotypes<-data.frame(Ecotypes_date[10:14], lapply(Ecotypes_date[1:9], B))
ecot_changes<-changes_relContr_ecotypes[,c(1,4:13)]

ecot_changes$Date<-as.Date(ecot_changes$Date,"%m/%d/%Y") 
ecot_changes_sorted<-ecot_changes[order(ecot_changes$Date),]

#Merge distances and Euclidean to Bray Curtis                            
Yy <- merge(wl_7day_BC_sorted, ecot_changes_sorted, by.x='distnamebckwd', by.y='distnameEco', all.x = TRUE)

Yy2<-subset(Yy, select=-c(rnames,lag.value,distnameR.x,row_day,col_day,days,row_Wf,col_Wf,Date.y,distnameR.y))

trial1<-subset(Yy2, Date_WaveLet>as.Date("2015-04-19") & Date_WaveLet<as.Date("2017-04-25"))
trial1_sorted<-trial1[order(trial1$Date_WaveLet),]

#linear interpolation for ecotypes
trial1_sorted$I.Ia.Ia.3.apx<-c("NA",na.approx(trial1_sorted$I.Ia.Ia.3)) 
trial1_sorted$I.Ia.Ia.1.apx<-c("NA",na.approx(trial1_sorted$I.Ia.Ia.1)) 
trial1_sorted$II.IIa.IIa.B.apx<-c("NA",na.approx(trial1_sorted$II.IIa.IIa.B)) 
trial1_sorted$IV.apx<-c("NA",na.approx(trial1_sorted$IV)) 
trial1_sorted$I.Ib.apx<-c("NA",na.approx(trial1_sorted$I.Ib)) 
trial1_sorted$I.apx<-c("NA",na.approx(trial1_sorted$I))
trial1_sorted$III.IIIa.apx<-c("NA",na.approx(trial1_sorted$III.IIIa)) 
trial1_sorted$II.apx<-c("NA",na.approx(trial1_sorted$II))  


#Heatmap of the normalized variables (use wl_7day_df_short_norm_difs)
#subset
wl_hm<-wl_7day_df_short_norm_interp[,c(5,8:17,20:23,25:27)]
wl_hm_melt<-melt(wl_hm, id.vars=c("Date_WaveLet"))
wl_hm_melt$Date_WaveLet <- as.Date(wl_hm_melt$Date_WaveLet,"%m/%d/%Y")

wl_hm_plot<-ggplot(data = wl_hm_melt, aes(x = Date_WaveLet, y=variable,fill=value)) + geom_tile(colour="white",size=0.25)+
scale_x_date( breaks = "1 month",minor_breaks = "1 week", date_labels = "%m/%y",limits = as.Date(c("2015-04-19","2017-04-15")))+scale_fill_gradient2(low = "grey", mid = "white", high = "brown", midpoint = .5)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=12,colour = "black"), axis.text.x=element_text(size=7,colour = "black"),legend.text=element_text(size=12))

#Divide wl_hm_melt in (a) strictly seasonal and (b) variable 

wl_hm_melt$variable<-factor(wl_hm_melt$variable,levels=rev(c("total.rainfall.for.week","average.weekly.river.flow","average.weekly.wave.height","average.weekly.wind.speed","maximum.wind.speed","average.weekly.pressure","salinity_PSU_2.5","oxygen_uM_2.5","temp_C_2.5","Density_kg.m3_2.5","PHOSPHATE_0","SILICATE_0","NITRATE.NIT_0","NITRITE_0","AMMONIA_0", "Par_uE.m2.s_2.5","Chl_0m")))        

wl_hm_plot<-ggplot(data = wl_hm_melt, aes(x = Date_WaveLet, y=variable,fill=value)) + geom_tile(colour="white",size=0.25)+
scale_x_date( breaks = "1 month",minor_breaks = "1 week", date_labels = "%m/%y",limits = as.Date(c("2015-04-19","2017-04-15")))+scale_fill_gradient2(low = "white", mid = "khaki2", high = "brown", midpoint = .5)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=16,colour = "black"), axis.text.x=element_text(size=8,colour = "black"),legend.text=element_text(size=12))+geom_vline(xintercept=as.Date(c("2016-01-01","2017-01-01")),linetype=1, colour="black",size=.5)

#shade by season and draw a line where years are

wl_hm_plot


#################
#####WAVELET#####
#################
#****trial1_sorted****# This file has all the distances we want to try. 

library(WaveletComp)

#############################Bray_Curtis#############################

my.date <- seq(as.POSIXct("2015-04-20 00:00:00", format = "%F %T"),by = "week",length.out = 106)

WeeklyBray<-as.numeric(trial1_sorted$Bray_Curtis.apx)

#Wavelet is kind of weird it only analyze taking into account date if that is the only extra column 

date_BC_wvl<-data.frame(date = my.date, Bray_Curtis.apx = WeeklyBray)

my.w <- analyze.wavelet(date_BC_wvl, "Bray_Curtis.apx",loess.span = 0,dt = 1, dj = 1/250, lowerPeriod = 1, upperPeriod = 64, make.pval = TRUE, n.sim = 1000) #1000 simulations

ticks <- seq(as.POSIXct("2015-04-20 00:00:00", format = "%F %T"), as.POSIXct("2017-04-25 23:00:00", format = "%F %T"), by = "month")labels <- seq(as.Date("2015-04-20"), as.Date("2017-04-25"), by = "month")
wt.image(my.w, periodlab = "periods (weeks)", legend.params = list(lab = "wavelet power levels"), label.time.axis = TRUE, show.date = TRUE, date.format = "%F %T",spec.time.axis = list(at = ticks, labels = labels, las = 2))


reconstruct(my.w, plot.waves = FALSE, lwd = c(1,2), legend.coords = "topleft",show.date = TRUE, date.format = "%F %T",spec.time.axis = list(at = ticks, labels = labels, las = 2))


#############################Euclidean Distances#############################

WeeklyEucl<-as.numeric(trial1_sorted$EucDist.All)
date_Euc_wvl<-data.frame(date = my.date, EucDist.All = WeeklyEucl)

my.y <- analyze.wavelet(date_Euc_wvl, "EucDist.All",loess.span = 0,dt = 1, dj = 1/250, lowerPeriod = 1, upperPeriod = 64, make.pval = TRUE, n.sim = 1000) #1000 simulations

wt.image(my.y, periodlab = "periods (weeks)", legend.params = list(lab = "wavelet power levels"), label.time.axis = TRUE, show.date = TRUE, date.format = "%F %T",spec.time.axis = list(at = ticks, labels = labels, las = 2))

reconstruct(my.y, plot.waves = FALSE, lwd = c(1,2), legend.coords = "topleft",show.date = TRUE, date.format = "%F %T",spec.time.axis = list(at = ticks, labels = labels, las = 2))

#I MERGED THE BRAY AND EUC DISTANCES IN INKSCAPE

###3 Analysis of a bivariate time series: Coherence analysis Bray-Eucl

# 1st merge date_BC_wvl and date_Euc_wvl

date_BC_Euc_wvl<-merge(date_BC_wvl, date_Euc_wvl, by="date", all=TRUE)

#my.wc <- analyze.coherency(date_BC_Euc_wvl, my.pair = c("Bray_Curtis.apx","EucDist.All"),loess.span = 0,dt = 1, dj = 1/250, lowerPeriod = 1, upperPeriod = 64, make.pval = TRUE, n.sim = 1000) #it takes about ~30 min

#saveRDS(my.wc, "/Users/luisbolanos/Documents/MisDrafts/InProgress/SAR11TS/Analysis/WEC/Wavelet/wc_BC_Euc.rds") 

####We can just start from here instead of estimating the coherence every time we run this script

my.wc<-readRDS("/Users/luisbolanos/Documents/MyDrafts/InProgress/SAR11TS/Analysis/WEC/Wavelet/wc_BC_Euc.rds")

#Cross-Wavelet power spectrum of the series with interval color key and restricted arrow area

wc.image(my.wc, n.levels = 250, color.key = "interval", siglvl.contour = 0.1, siglvl.arrow = 0.05, which.arrow.sig = "wt", legend.params = list(lab = "cross-wavelet power levels"),show.date = TRUE, date.format = "%F %T",spec.time.axis = list(at = ticks, labels = labels, las = 2))



#Phases and phase differences at period 16, default phase axis phase axis labels expressed in time units
#Phase distribution of the two time-series at the 16 week period

wc.sel.phases(my.wc, sel.period = 16,only.sig = TRUE,siglvl = 0.05,phaselim = c(-pi,+pi),legend.coords = "topleft", legend.horiz = FALSE,main = "", sub = "", show.date = TRUE, date.format = "%F %T",spec.time.axis = list(at = ticks, labels = labels, las = 2))


#Cross-wavelet average power
wc.avg(my.wc, siglvl = 0.01, sigcol = "red", sigpch = 20,periodlab = "period (weeks)")
