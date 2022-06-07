#This R script reproduce the Figure X: Linear models of the different SAR11 ASVS from the Bermuda Atlantic Time Series and The Western English Channel.

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
library("stringr")
library("directlabels")


#####Create the PhyloSeq object####

countwec <- read.table("WEC_mar.otu", header=T, row.names=1, check.names=F)
sampleinfowec <- read.table("WEC_mar.env", header=T, row.names=1, check.names=F, sep ="\t")
taxwec <- as.matrix(read.table("WEC_marV2.tax", header=T, row.names=1, check.names=F, na.strings="", sep="\t")) 

OTUmar = otu_table(countwec, taxa_are_rows = TRUE)
TAXmar = tax_table(taxwec)
SAMar= sample_data(sampleinfowec)

WEC<-phyloseq(OTUmar,TAXmar,SAMar)
WECcf= filter_taxa(WEC, function(x) sum(x > 1) > (0.0015*length(x)), TRUE)

#####WEC_SAR11####
SAR11_WECphy = subset_taxa(WECcf, Order=="SAR11_clade") #264 samples 

#####Remove samples with less than 1000 SAR11 reads and convert to %#####

SAR11_WECphy1k<-prune_samples(sample_sums(SAR11_WECphy)>1000, SAR11_WECphy) #242 samples
SAR11_WECphy1krel<-transform_sample_counts(SAR11_WECphy1k, function(x){x / sum(x)})

###Get Sample-ASVs-relAbund-Tax####
#TAXA
SAR11tax<-as.data.frame(tax_table(SAR11_WECphy1krel)[,6])
SAR11tax$ASV<-rownames(SAR11tax)

#OTU
SAR11_ASV<-as.data.frame(otu_table(SAR11_WECphy1krel))
SAR11_ASVt<-t(SAR11_ASV)

#####Create a data.frame cols=ASV, rows= samples -> add date to this and estimate the sin/cos thing 
SAR11tax<-as.data.frame(tax_table(SAR11_WECphy1krel)[,6]) ###Keep this one to add to the melted object

SAR11_ASV<-as.data.frame(otu_table(SAR11_WECphy1krel))ASV_frw<-as.data.frame(x = t(SAR11_ASV), stringsAsFactors = FALSE)md_to_add<-as.data.frame(sample_data(SAR11_WECphy1krel))[,c(4)] #ID,Date, DateEnvfinal_2a<-cbind(ASV_frw,md_to_add) #I'll use Date to estimate stuff
final_2a$Date<-as.Date(final_2a$Date,"%m/%d/%Y")

final_2a_sort<-final_2a[order(final_2a$Date),] ##sort from least recent to most recent
final_2a_sort$Date<-as.Date(final_2a_sort$Date,"%Y-%m-%d")

##Add the Julian day 
final_2a_sort$Jul<-format(final_2a_sort$Date, "%j")

###Serial Day
#Starting from
WECTS_start<- as.Date('2012-01-01')

final_2a_sort$Serial<- as.numeric(-difftime(WECTS_start,final_2a_sort$Date)+1)

##Days since winter solstice (21 Dic)

#Add the winter solstice (Basically you could arrange it for any date you may want to explore
final_2a_sort$wintersols[final_2a_sort$Date>as.Date("2012-01-01") & final_2a_sort$Date<as.Date("2012-12-21")] <- "2011-12-21"
final_2a_sort$wintersols[final_2a_sort$Date>as.Date("2013-01-01") & final_2a_sort$Date<as.Date("2013-12-21")] <- "2012-12-21"
final_2a_sort$wintersols[final_2a_sort$Date>as.Date("2014-01-01") & final_2a_sort$Date<as.Date("2014-12-21")] <- "2013-12-21"
final_2a_sort$wintersols[final_2a_sort$Date>as.Date("2015-01-01") & final_2a_sort$Date<as.Date("2015-12-21")] <- "2014-12-21"
final_2a_sort$wintersols[final_2a_sort$Date>as.Date("2016-01-01") & final_2a_sort$Date<as.Date("2016-12-21")] <- "2015-12-21"
final_2a_sort$wintersols[final_2a_sort$Date>as.Date("2017-01-01") & final_2a_sort$Date<as.Date("2017-12-21")] <- "2016-12-21"
final_2a_sort$wintersols[final_2a_sort$Date>as.Date("2018-01-01") & final_2a_sort$Date<as.Date("2018-12-21")] <- "2017-12-21" #The most recent date is 2018-12-17

final_2a_sort$DaysSincewintersols<-as.numeric(-difftime(final_2a_sort$wintersols,final_2a_sort$Date)+1) #Add the numeric column

#jultonumeric
final_2a_sort$Jul<-as.numeric(final_2a_sort$Jul)

#models based on the sin and cos
final_2a_sort$cosDX1<-as.numeric(cos(2*pi*((final_2a_sort$DaysSincewintersols)/365))) #Peaking in midwinter
final_2a_sort$sinDX1<-as.numeric(sin(2*pi*((final_2a_sort$Jul)/365))) #Peaking other time

#Get the linear model for all the columns
#FIRST: loop which regresses each ASV column to cosDX and save these regressions on a list (1:144) 

#[145] "Date" [146] "Jul"[147]"Serial" [148]"wintersols" [149]"DaysSincewintersols" [150] "cosDX1" [151] "senDX1"

storage <- list()
for(i in names(final_2a_sort)[1:144]){
  storage[[i]] <- lm(get(i) ~ sinDX1+cosDX1, final_2a_sort)
}

#Create a data set with the regressions
fitted.values(storage[[1]])

#Get the significant regressions (THIS MUST BE THE SAME ASVS as I did with g.fisher
pvalues_lm<-sapply(storage, function(x){
  ff <- summary(x)$fstatistic
  pf(ff[1], df1 = ff[2], df2 = ff[3], lower.tail = FALSE)
})

pvalues_lm05<-pvalues_lm[pvalues_lm < 0.05] #The names here should be ~ to g.fisher #103/144 = 71%

#All the significant 
sign_names<-str_remove(names(pvalues_lm05), ".value") #extract this from both datasets

#Number of the column that matches a list of names get the position of the names <0.05
signll<-match(sign_names,names(final_2a_sort))

#subset list of lists 
storage05<-storage[signll]

#I need the constant, sin and cos of the models, the cosDX and SinDX originals and
#Get the constant, sin and cos of the models:
coeffs05<-data.frame(coeffs=sapply(storage05, FUN=function(item){item$coefficients}))
colnames(coeffs05)<-str_remove(colnames(coeffs05), "coeffs.")

#Subset original day, sin, cos
newdf_model <- final_2a_sort[, c("Serial","cosDX1","sinDX1")]

#automatizise this :P 
newdf_model$WEC1<-coeffs05$WEC1[1]+(newdf_model$sinDX1*coeffs05$WEC1[2])+(newdf_model$cosDX1*coeffs05$WEC1[3])
newdf_model$WEC3<-coeffs05$WEC3[1]+(newdf_model$sinDX1*coeffs05$WEC3[2])+(newdf_model$cosDX1*coeffs05$WEC3[3])
newdf_model$WEC6<-coeffs05$WEC6[1]+(newdf_model$sinDX1*coeffs05$WEC6[2])+(newdf_model$cosDX1*coeffs05$WEC6[3])
newdf_model$WEC7<-coeffs05$WEC7[1]+(newdf_model$sinDX1*coeffs05$WEC7[2])+(newdf_model$cosDX1*coeffs05$WEC7[3])
newdf_model$WEC44<-coeffs05$WEC44[1]+(newdf_model$sinDX1*coeffs05$WEC44[2])+(newdf_model$cosDX1*coeffs05$WEC44[3])
newdf_model$WEC79<-coeffs05$WEC79[1]+(newdf_model$sinDX1*coeffs05$WEC79[2])+(newdf_model$cosDX1*coeffs05$WEC79[3])
newdf_model$WEC109<-coeffs05$WEC109[1]+(newdf_model$sinDX1*coeffs05$WEC109[2])+(newdf_model$cosDX1*coeffs05$WEC109[3])
newdf_model$WEC131<-coeffs05$WEC131[1]+(newdf_model$sinDX1*coeffs05$WEC131[2])+(newdf_model$cosDX1*coeffs05$WEC131[3])
newdf_model$WEC138<-coeffs05$WEC138[1]+(newdf_model$sinDX1*coeffs05$WEC138[2])+(newdf_model$cosDX1*coeffs05$WEC138[3])
newdf_model$WEC140<-coeffs05$WEC140[1]+(newdf_model$sinDX1*coeffs05$WEC140[2])+(newdf_model$cosDX1*coeffs05$WEC140[3])
newdf_model$WEC142<-coeffs05$WEC142[1]+(newdf_model$sinDX1*coeffs05$WEC142[2])+(newdf_model$cosDX1*coeffs05$WEC142[3])
newdf_model$WEC183<-coeffs05$WEC183[1]+(newdf_model$sinDX1*coeffs05$WEC183[2])+(newdf_model$cosDX1*coeffs05$WEC183[3])
newdf_model$WEC191<-coeffs05$WEC191[1]+(newdf_model$sinDX1*coeffs05$WEC191[2])+(newdf_model$cosDX1*coeffs05$WEC191[3])
newdf_model$WEC197<-coeffs05$WEC197[1]+(newdf_model$sinDX1*coeffs05$WEC197[2])+(newdf_model$cosDX1*coeffs05$WEC197[3])
newdf_model$WEC201<-coeffs05$WEC201[1]+(newdf_model$sinDX1*coeffs05$WEC201[2])+(newdf_model$cosDX1*coeffs05$WEC201[3])
newdf_model$WEC206<-coeffs05$WEC206[1]+(newdf_model$sinDX1*coeffs05$WEC206[2])+(newdf_model$cosDX1*coeffs05$WEC206[3])
newdf_model$WEC286<-coeffs05$WEC286[1]+(newdf_model$sinDX1*coeffs05$WEC286[2])+(newdf_model$cosDX1*coeffs05$WEC286[3])
newdf_model$WEC302<-coeffs05$WEC302[1]+(newdf_model$sinDX1*coeffs05$WEC302[2])+(newdf_model$cosDX1*coeffs05$WEC302[3])
newdf_model$WEC319<-coeffs05$WEC319[1]+(newdf_model$sinDX1*coeffs05$WEC319[2])+(newdf_model$cosDX1*coeffs05$WEC319[3])
newdf_model$WEC321<-coeffs05$WEC321[1]+(newdf_model$sinDX1*coeffs05$WEC321[2])+(newdf_model$cosDX1*coeffs05$WEC321[3])
newdf_model$WEC332<-coeffs05$WEC332[1]+(newdf_model$sinDX1*coeffs05$WEC332[2])+(newdf_model$cosDX1*coeffs05$WEC332[3])
newdf_model$WEC343<-coeffs05$WEC343[1]+(newdf_model$sinDX1*coeffs05$WEC343[2])+(newdf_model$cosDX1*coeffs05$WEC343[3])
newdf_model$WEC360<-coeffs05$WEC360[1]+(newdf_model$sinDX1*coeffs05$WEC360[2])+(newdf_model$cosDX1*coeffs05$WEC360[3])
newdf_model$WEC371<-coeffs05$WEC371[1]+(newdf_model$sinDX1*coeffs05$WEC371[2])+(newdf_model$cosDX1*coeffs05$WEC371[3])
newdf_model$WEC384<-coeffs05$WEC384[1]+(newdf_model$sinDX1*coeffs05$WEC384[2])+(newdf_model$cosDX1*coeffs05$WEC384[3])
newdf_model$WEC386<-coeffs05$WEC386[1]+(newdf_model$sinDX1*coeffs05$WEC386[2])+(newdf_model$cosDX1*coeffs05$WEC386[3])
newdf_model$WEC396<-coeffs05$WEC396[1]+(newdf_model$sinDX1*coeffs05$WEC396[2])+(newdf_model$cosDX1*coeffs05$WEC396[3])
newdf_model$WEC402<-coeffs05$WEC402[1]+(newdf_model$sinDX1*coeffs05$WEC402[2])+(newdf_model$cosDX1*coeffs05$WEC402[3])
newdf_model$WEC420<-coeffs05$WEC420[1]+(newdf_model$sinDX1*coeffs05$WEC420[2])+(newdf_model$cosDX1*coeffs05$WEC420[3])
newdf_model$WEC435<-coeffs05$WEC435[1]+(newdf_model$sinDX1*coeffs05$WEC435[2])+(newdf_model$cosDX1*coeffs05$WEC435[3])
newdf_model$WEC440<-coeffs05$WEC440[1]+(newdf_model$sinDX1*coeffs05$WEC440[2])+(newdf_model$cosDX1*coeffs05$WEC440[3])
newdf_model$WEC523<-coeffs05$WEC523[1]+(newdf_model$sinDX1*coeffs05$WEC523[2])+(newdf_model$cosDX1*coeffs05$WEC523[3])
newdf_model$WEC524<-coeffs05$WEC524[1]+(newdf_model$sinDX1*coeffs05$WEC524[2])+(newdf_model$cosDX1*coeffs05$WEC524[3])
newdf_model$WEC571<-coeffs05$WEC571[1]+(newdf_model$sinDX1*coeffs05$WEC571[2])+(newdf_model$cosDX1*coeffs05$WEC571[3])
newdf_model$WEC573<-coeffs05$WEC573[1]+(newdf_model$sinDX1*coeffs05$WEC573[2])+(newdf_model$cosDX1*coeffs05$WEC573[3])
newdf_model$WEC576<-coeffs05$WEC576[1]+(newdf_model$sinDX1*coeffs05$WEC576[2])+(newdf_model$cosDX1*coeffs05$WEC576[3])
newdf_model$WEC581<-coeffs05$WEC581[1]+(newdf_model$sinDX1*coeffs05$WEC581[2])+(newdf_model$cosDX1*coeffs05$WEC581[3])
newdf_model$WEC584<-coeffs05$WEC584[1]+(newdf_model$sinDX1*coeffs05$WEC584[2])+(newdf_model$cosDX1*coeffs05$WEC584[3])
newdf_model$WEC637<-coeffs05$WEC637[1]+(newdf_model$sinDX1*coeffs05$WEC637[2])+(newdf_model$cosDX1*coeffs05$WEC637[3])
newdf_model$WEC650<-coeffs05$WEC650[1]+(newdf_model$sinDX1*coeffs05$WEC650[2])+(newdf_model$cosDX1*coeffs05$WEC650[3])
newdf_model$WEC686<-coeffs05$WEC686[1]+(newdf_model$sinDX1*coeffs05$WEC686[2])+(newdf_model$cosDX1*coeffs05$WEC686[3])
newdf_model$WEC703<-coeffs05$WEC703[1]+(newdf_model$sinDX1*coeffs05$WEC703[2])+(newdf_model$cosDX1*coeffs05$WEC703[3])
newdf_model$WEC720<-coeffs05$WEC720[1]+(newdf_model$sinDX1*coeffs05$WEC720[2])+(newdf_model$cosDX1*coeffs05$WEC720[3])
newdf_model$WEC784<-coeffs05$WEC784[1]+(newdf_model$sinDX1*coeffs05$WEC784[2])+(newdf_model$cosDX1*coeffs05$WEC784[3])
newdf_model$WEC796<-coeffs05$WEC796[1]+(newdf_model$sinDX1*coeffs05$WEC796[2])+(newdf_model$cosDX1*coeffs05$WEC796[3])
newdf_model$WEC828<-coeffs05$WEC828[1]+(newdf_model$sinDX1*coeffs05$WEC828[2])+(newdf_model$cosDX1*coeffs05$WEC828[3])
newdf_model$WEC841<-coeffs05$WEC841[1]+(newdf_model$sinDX1*coeffs05$WEC841[2])+(newdf_model$cosDX1*coeffs05$WEC841[3])
newdf_model$WEC846<-coeffs05$WEC846[1]+(newdf_model$sinDX1*coeffs05$WEC846[2])+(newdf_model$cosDX1*coeffs05$WEC846[3])
newdf_model$WEC850<-coeffs05$WEC850[1]+(newdf_model$sinDX1*coeffs05$WEC850[2])+(newdf_model$cosDX1*coeffs05$WEC850[3])
newdf_model$WEC857<-coeffs05$WEC857[1]+(newdf_model$sinDX1*coeffs05$WEC857[2])+(newdf_model$cosDX1*coeffs05$WEC857[3])
newdf_model$WEC918<-coeffs05$WEC918[1]+(newdf_model$sinDX1*coeffs05$WEC918[2])+(newdf_model$cosDX1*coeffs05$WEC918[3])
newdf_model$WEC932<-coeffs05$WEC932[1]+(newdf_model$sinDX1*coeffs05$WEC932[2])+(newdf_model$cosDX1*coeffs05$WEC932[3])
newdf_model$WEC948<-coeffs05$WEC948[1]+(newdf_model$sinDX1*coeffs05$WEC948[2])+(newdf_model$cosDX1*coeffs05$WEC948[3])
newdf_model$WEC1010<-coeffs05$WEC1010[1]+(newdf_model$sinDX1*coeffs05$WEC1010[2])+(newdf_model$cosDX1*coeffs05$WEC1010[3])
newdf_model$WEC1034<-coeffs05$WEC1034[1]+(newdf_model$sinDX1*coeffs05$WEC1034[2])+(newdf_model$cosDX1*coeffs05$WEC1034[3])
newdf_model$WEC1070<-coeffs05$WEC1070[1]+(newdf_model$sinDX1*coeffs05$WEC1070[2])+(newdf_model$cosDX1*coeffs05$WEC1070[3])
newdf_model$WEC1095<-coeffs05$WEC1095[1]+(newdf_model$sinDX1*coeffs05$WEC1095[2])+(newdf_model$cosDX1*coeffs05$WEC1095[3])
newdf_model$WEC1125<-coeffs05$WEC1125[1]+(newdf_model$sinDX1*coeffs05$WEC1125[2])+(newdf_model$cosDX1*coeffs05$WEC1125[3])
newdf_model$WEC1146<-coeffs05$WEC1146[1]+(newdf_model$sinDX1*coeffs05$WEC1146[2])+(newdf_model$cosDX1*coeffs05$WEC1146[3])
newdf_model$WEC1197<-coeffs05$WEC1197[1]+(newdf_model$sinDX1*coeffs05$WEC1197[2])+(newdf_model$cosDX1*coeffs05$WEC1197[3])
newdf_model$WEC1302<-coeffs05$WEC1302[1]+(newdf_model$sinDX1*coeffs05$WEC1302[2])+(newdf_model$cosDX1*coeffs05$WEC1302[3])
newdf_model$WEC1320<-coeffs05$WEC1320[1]+(newdf_model$sinDX1*coeffs05$WEC1320[2])+(newdf_model$cosDX1*coeffs05$WEC1320[3])
newdf_model$WEC1361<-coeffs05$WEC1361[1]+(newdf_model$sinDX1*coeffs05$WEC1361[2])+(newdf_model$cosDX1*coeffs05$WEC1361[3])
newdf_model$WEC1363<-coeffs05$WEC1363[1]+(newdf_model$sinDX1*coeffs05$WEC1363[2])+(newdf_model$cosDX1*coeffs05$WEC1363[3])
newdf_model$WEC1375<-coeffs05$WEC1375[1]+(newdf_model$sinDX1*coeffs05$WEC1375[2])+(newdf_model$cosDX1*coeffs05$WEC1375[3])
newdf_model$WEC1415<-coeffs05$WEC1415[1]+(newdf_model$sinDX1*coeffs05$WEC1415[2])+(newdf_model$cosDX1*coeffs05$WEC1415[3])
newdf_model$WEC1441<-coeffs05$WEC1441[1]+(newdf_model$sinDX1*coeffs05$WEC1441[2])+(newdf_model$cosDX1*coeffs05$WEC1441[3])
newdf_model$WEC1447<-coeffs05$WEC1447[1]+(newdf_model$sinDX1*coeffs05$WEC1447[2])+(newdf_model$cosDX1*coeffs05$WEC1447[3])
newdf_model$WEC1462<-coeffs05$WEC1462[1]+(newdf_model$sinDX1*coeffs05$WEC1462[2])+(newdf_model$cosDX1*coeffs05$WEC1462[3])
newdf_model$WEC1530<-coeffs05$WEC1530[1]+(newdf_model$sinDX1*coeffs05$WEC1530[2])+(newdf_model$cosDX1*coeffs05$WEC1530[3])
newdf_model$WEC1611<-coeffs05$WEC1611[1]+(newdf_model$sinDX1*coeffs05$WEC1611[2])+(newdf_model$cosDX1*coeffs05$WEC1611[3])
newdf_model$WEC1665<-coeffs05$WEC1665[1]+(newdf_model$sinDX1*coeffs05$WEC1665[2])+(newdf_model$cosDX1*coeffs05$WEC1665[3])
newdf_model$WEC1738<-coeffs05$WEC1738[1]+(newdf_model$sinDX1*coeffs05$WEC1738[2])+(newdf_model$cosDX1*coeffs05$WEC1738[3])
newdf_model$WEC1801<-coeffs05$WEC1801[1]+(newdf_model$sinDX1*coeffs05$WEC1801[2])+(newdf_model$cosDX1*coeffs05$WEC1801[3])
newdf_model$WEC1856<-coeffs05$WEC1856[1]+(newdf_model$sinDX1*coeffs05$WEC1856[2])+(newdf_model$cosDX1*coeffs05$WEC1856[3])
newdf_model$WEC1909<-coeffs05$WEC1909[1]+(newdf_model$sinDX1*coeffs05$WEC1909[2])+(newdf_model$cosDX1*coeffs05$WEC1909[3])
newdf_model$WEC2033<-coeffs05$WEC2033[1]+(newdf_model$sinDX1*coeffs05$WEC2033[2])+(newdf_model$cosDX1*coeffs05$WEC2033[3])
newdf_model$WEC2057<-coeffs05$WEC2057[1]+(newdf_model$sinDX1*coeffs05$WEC2057[2])+(newdf_model$cosDX1*coeffs05$WEC2057[3])
newdf_model$WEC2163<-coeffs05$WEC2163[1]+(newdf_model$sinDX1*coeffs05$WEC2163[2])+(newdf_model$cosDX1*coeffs05$WEC2163[3])
newdf_model$WEC2205<-coeffs05$WEC2205[1]+(newdf_model$sinDX1*coeffs05$WEC2205[2])+(newdf_model$cosDX1*coeffs05$WEC2205[3])
newdf_model$WEC2275<-coeffs05$WEC2275[1]+(newdf_model$sinDX1*coeffs05$WEC2275[2])+(newdf_model$cosDX1*coeffs05$WEC2275[3])
newdf_model$WEC2313<-coeffs05$WEC2313[1]+(newdf_model$sinDX1*coeffs05$WEC2313[2])+(newdf_model$cosDX1*coeffs05$WEC2313[3])
newdf_model$WEC2517<-coeffs05$WEC2517[1]+(newdf_model$sinDX1*coeffs05$WEC2517[2])+(newdf_model$cosDX1*coeffs05$WEC2517[3])
newdf_model$WEC2688<-coeffs05$WEC2688[1]+(newdf_model$sinDX1*coeffs05$WEC2688[2])+(newdf_model$cosDX1*coeffs05$WEC2688[3])
newdf_model$WEC2692<-coeffs05$WEC2692[1]+(newdf_model$sinDX1*coeffs05$WEC2692[2])+(newdf_model$cosDX1*coeffs05$WEC2692[3])
newdf_model$WEC2874<-coeffs05$WEC2874[1]+(newdf_model$sinDX1*coeffs05$WEC2874[2])+(newdf_model$cosDX1*coeffs05$WEC2874[3])
newdf_model$WEC2960<-coeffs05$WEC2960[1]+(newdf_model$sinDX1*coeffs05$WEC2960[2])+(newdf_model$cosDX1*coeffs05$WEC2960[3])
newdf_model$WEC3242<-coeffs05$WEC3242[1]+(newdf_model$sinDX1*coeffs05$WEC3242[2])+(newdf_model$cosDX1*coeffs05$WEC3242[3])
newdf_model$WEC3325<-coeffs05$WEC3325[1]+(newdf_model$sinDX1*coeffs05$WEC3325[2])+(newdf_model$cosDX1*coeffs05$WEC3325[3])
newdf_model$WEC3360<-coeffs05$WEC3360[1]+(newdf_model$sinDX1*coeffs05$WEC3360[2])+(newdf_model$cosDX1*coeffs05$WEC3360[3])
newdf_model$WEC3609<-coeffs05$WEC3609[1]+(newdf_model$sinDX1*coeffs05$WEC3609[2])+(newdf_model$cosDX1*coeffs05$WEC3609[3])
newdf_model$WEC3749<-coeffs05$WEC3749[1]+(newdf_model$sinDX1*coeffs05$WEC3749[2])+(newdf_model$cosDX1*coeffs05$WEC3749[3])
newdf_model$WEC3877<-coeffs05$WEC3877[1]+(newdf_model$sinDX1*coeffs05$WEC3877[2])+(newdf_model$cosDX1*coeffs05$WEC3877[3])
newdf_model$WEC4009<-coeffs05$WEC4009[1]+(newdf_model$sinDX1*coeffs05$WEC4009[2])+(newdf_model$cosDX1*coeffs05$WEC4009[3])
newdf_model$WEC4308<-coeffs05$WEC4308[1]+(newdf_model$sinDX1*coeffs05$WEC4308[2])+(newdf_model$cosDX1*coeffs05$WEC4308[3])
newdf_model$WEC4909<-coeffs05$WEC4909[1]+(newdf_model$sinDX1*coeffs05$WEC4909[2])+(newdf_model$cosDX1*coeffs05$WEC4909[3])
newdf_model$WEC5154<-coeffs05$WEC5154[1]+(newdf_model$sinDX1*coeffs05$WEC5154[2])+(newdf_model$cosDX1*coeffs05$WEC5154[3])
newdf_model$WEC5652<-coeffs05$WEC5652[1]+(newdf_model$sinDX1*coeffs05$WEC5652[2])+(newdf_model$cosDX1*coeffs05$WEC5652[3])
newdf_model$WEC6115<-coeffs05$WEC6115[1]+(newdf_model$sinDX1*coeffs05$WEC6115[2])+(newdf_model$cosDX1*coeffs05$WEC6115[3])
newdf_model$WEC8824<-coeffs05$WEC8824[1]+(newdf_model$sinDX1*coeffs05$WEC8824[2])+(newdf_model$cosDX1*coeffs05$WEC8824[3])
newdf_model$WEC25800<-coeffs05$WEC25800[1]+(newdf_model$sinDX1*coeffs05$WEC25800[2])+(newdf_model$cosDX1*coeffs05$WEC25800[3])
newdf_model$WEC26018<-coeffs05$WEC26018[1]+(newdf_model$sinDX1*coeffs05$WEC26018[2])+(newdf_model$cosDX1*coeffs05$WEC26018[3])
newdf_model$WEC26223<-coeffs05$WEC26223[1]+(newdf_model$sinDX1*coeffs05$WEC26223[2])+(newdf_model$cosDX1*coeffs05$WEC26223[3])



#newdf_model for melt 

Tomelt<-(newdf_model[c(1,4:106)])

#Add taxa to this melted 
newdf_model_melted<-melt(Tomelt,id.vars="Serial")

#Agregar taxa extracting sign_names variable from the Phyloseq object, thereafter extract the tax table only with genus 

##### Subset from physeq sar11.

phy_sign05 = prune_taxa(sign_names, SAR11_WECphy)

GenusSAR11sign1<-data.frame(tax_table(phy_sign05))[,6,drop=FALSE]

GenusSAR11sign1$ASV<-rownames(GenusSAR11sign1) #GenusSAR11sign1 contains the neccesary to add taxa to melted

#change variable to ASV in newdf_model_melted 
names(newdf_model_melted)[2]<-"ASV"

#merge by ASV and add the column Genus 
newdf_model_melted_taxa<-merge(newdf_model_melted, GenusSAR11sign1[, c("Genus", "ASV")], by="ASV")

#Subset 1 to 365 
newdf_model_melted_taxa365<-newdf_model_melted_taxa[newdf_model_melted_taxa$Serial >1095 & newdf_model_melted_taxa$Serial <1460,]
newdf_model_melted_taxa365$Serial365<-(newdf_model_melted_taxa365$Serial-1094)

#plot
#unique(newdf_model_melted_taxa$Genus) = "I;Ia,Ia3","II;IIb","cladeN2","II;IIa;IIa.B" ,"II;IIa;IIa_NA","II","IV","I","I;Ib","I;Ia,Ia.1","I;Ic;Ic1","I;Ia,Ia.4","III;IIIa","II;IIa;IIa.A","I;Ic;Ic.2"

#Colors for SAR11 subclasses
coloresbarplot = c("I;Ia,Ia.1"="darkcyan","I;Ia,Ia.3"="blue","I"="black","I;Ia,Ia.4"="lightskyblue","II;IIa;IIa.A"="bisque4","II;IIa;IIa.B"="darkgoldenrod3","II;IIa;IIa_NA"="gold2","III;IIIa"="cornsilk3","II"="chocolate4","II;IIb"="olivedrab","IV"="darkgreen","I;Ic;Ic.2"="lightcoral","I;Ic;Ic.1"="hotpink","clade N2"="lavenderblush2", "I;Ib"="#d45d41")

#I will split the panels in top and bottom as the heat map. The only difference is that I´ll use 10 as value to divide. The 4 most abundant are in a larger scale

d<-colnames(SAR11_ASVt)[apply(SAR11_ASVt, 2, function(u) colSums(SAR11_ASVt)>10)] #Get the names of most abundant
d <- d[!is.na(d)] 
a<-colnames(SAR11_ASVt)[apply(SAR11_ASVt, 2, function(u) colSums(SAR11_ASVt)<=10)& colSums(SAR11_ASVt)>0.5] #Get the names of least abundant
a <- a[!is.na(a)]
t<-colnames(SAR11_ASVt)[apply(SAR11_ASVt, 2, function(u) colSums(SAR11_ASVt)<=.05)] #Get the names of least abundant
t<- t[!is.na(t)]

###Subset from the melted object

highSAR11<- newdf_model_melted_taxa365[newdf_model_melted_taxa365$ASV %in% d, ]
middleSAR11<- newdf_model_melted_taxa365[newdf_model_melted_taxa365$ASV %in% a, ]
lowSAR11<-newdf_model_melted_taxa365[newdf_model_melted_taxa365$ASV %in% t, ]

#Plot
myplothigh<-ggplot(highSAR11, aes(x = Serial365, y =  value, colour=Genus, group = ASV, label=ASV)) +geom_line(size=.8) + 
  directlabels::geom_dl(aes(label = ASV),  method = list(cex = .9, "smart.grid"))+ scale_color_manual(values = coloresbarplot)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=14,colour = "black"), axis.text.x=element_text(size=14,colour = "black"),legend.text=element_text(size=12))+labs(x = "day", y = "relative contribution")+ guides(color = guide_legend(override.aes = list(size = 3.5)))

myplotmid<-ggplot(middleSAR11, aes(x = Serial365, y =  value, colour=Genus, group = ASV, label=ASV)) + geom_line(size=.8) + scale_color_manual(values = coloresbarplot)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=14,colour = "black"), axis.text.x=element_text(size=14,colour = "black"),legend.text=element_text(size=12))+labs(x = "day", y = "relative contribution")+ guides(color = guide_legend(override.aes = list(size = 3.5)))

myplotlow<-ggplot(lowSAR11, aes(x = Serial365, y =  value, colour=Genus, group = ASV, label=ASV)) + geom_line(size=.8) + scale_color_manual(values = coloresbarplot)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_text(size=14,colour = "black"),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=14,colour = "black"), axis.text.x=element_text(size=14,colour = "black"),legend.text=element_text(size=12))+labs(x = "day", y = "relative contribution")+ guides(color = guide_legend(override.aes = list(size = 3.5)))

meds_WEC<-plot_grid(myplothigh+theme(legend.position="none"),myplotmid+theme(legend.position="none"),myplotlow+theme(legend.position="none"), ncol=1, align="hv", axis = "bl")

#Get legend
legend <- get_legend(
  # create some space to the left of the legend
  myplotlow + theme(legend.box.margin = margin(0, 0, 0, 12))
)

WEC_model_plot<-plot_grid(meds_WEC, legend, rel_widths = c(2.2, .5))

############################
############BATS############
############################

#The below code replicates the linear models use above for WEC but now for BATS 2016-2018

count_tabS11 <- read.table("/Users/luisbolanos/Documents/MyDrafts/InProgress/SAR11TS/Github_submission/BIOSotu5.otu", header=T, row.names=1, check.names=F)
sample_info_tabS11 <- read.table("/Users/luisbolanos/Documents/MyDrafts/InProgress/SAR11TS/Github_submission/BIOSenvts5.env", header=T, row.names=1, check.names=F, sep ="\t")
tax_tabS11 <- as.matrix(read.table("/Users/luisbolanos/Documents/MyDrafts/InProgress/SAR11TS/Github_submission/BIOStax5.tax", header=T, row.names=1, check.names=F, sep="\t",na.strings = "#NA")) #No dejar espacios en blanco

sample_info_tabS11$Date <- as.Date(with(sample_info_tabS11, paste(year, month, day,sep="-")), "%Y-%m-%d")

#Add tech label
sample_info_tabS11$tech<-ifelse(sample_info_tabS11$year>2015, "MiSeq", "FLX454") 

S11OTU = otu_table(count_tabS11, taxa_are_rows = TRUE)
S11TAX = tax_table(tax_tabS11)
S11SAM= sample_data(sample_info_tabS11)

S11phy<-phyloseq(S11OTU,S11TAX,S11SAM)

ts = get_variable(S11phy, "Type") %in% "BATS"sample_data(S11phy)$ts <- factor(ts)phyTS<-subset_samples(S11phy, ts %in% TRUE)

phyTS_prun<-prune_taxa(taxa_sums(phyTS) > 0, phyTS) # De 3193 a 2875 ASVs (318)

##Get all the 5m samples####
ts5 = get_variable(phyTS_prun, "Depth_C") %in% "0"sample_data(phyTS_prun)$ts5 <- factor(ts5)phyTS5m<-subset_samples(phyTS_prun, ts5 %in% TRUE) ###163 samples

phyTS5mV1 <- prune_samples(sample_sums(phyTS5m)>=1000, phyTS5m) #prune samples less than 1000 reads (eliminate 17 --> 146/163)

phyTS5mpruneV2<-filter_taxa(phyTS5mV1, function(x) sum(x >= 3) > (0.02 *length(x)), TRUE) #V1=916/2875 Again from phyTS5mprune, V2 (146 samples, .01) =500 ASVs

SAR11B5 = subset_taxa(phyTS5mpruneV2, Order=="SAR11_clade")

SAR11_BATSphyMiSeq<-subset_samples(SAR11B5, tech=="MiSeq") # Myseq 2016-2018

SAR11_BATSphyMiSeqrel<-transform_sample_counts(SAR11_BATSphyMiSeq, function(x){x / sum(x)}) #Transform to % only based on the contribution to total sar11

#TAXA
SAR11_BATStax<-as.data.frame(tax_table(SAR11_BATSphyMiSeqrel)[,8])
SAR11_BATStax$ASV<-rownames(SAR11_BATStax)

#OTU
SAR11_BATS_ASV<-as.data.frame(otu_table(SAR11_BATSphyMiSeqrel))
SAR11_BATS_ASVt<-t(SAR11_BATS_ASV)

#####Create a data.frame cols=ASV, rows= samples -> add date to this and estimate the sin/cos thing 
SAR11_BATStax<-as.data.frame(tax_table(SAR11_BATSphyMiSeqrel)[,8]) ###Keep this one to add to the melted object

ASV_frw_BATS<-as.data.frame(x = t(SAR11_BATS_ASV), stringsAsFactors = FALSE)

ASV_frw_BATS<-as.data.frame(x = t(SAR11_BATS_ASV), stringsAsFactors = FALSE)md_to_add_BATS<-as.data.frame(sample_data(SAR11_BATSphyMiSeqrel))[,c(31)] #ID,Date, DateEnvfinal_2a_BATS<-cbind(ASV_frw_BATS,md_to_add_BATS) #I'll use Date to estimate stuff
final_2a_BATS$Date<-as.Date(final_2a_BATS$Date,"%m/%d/%Y")

final_2a_sort_BATS<-final_2a_BATS[order(final_2a_BATS$Date),] ##sort from least recent to most recent
final_2a_sort_BATS$Date<-as.Date(final_2a_sort_BATS$Date,"%Y-%m-%d")

##Add the Julian day 
final_2a_sort_BATS$Jul<-format(final_2a_sort_BATS$Date, "%j")

###Serial Day
#Starting from
BATS_start<- as.Date('2016-01-01')
final_2a_sort_BATS$Serial<- as.numeric(-difftime(BATS_start,final_2a_sort_BATS$Date)+1)

##Days since winter solstice (21 Dic)
#Add the winter solstice (Basically you could arrange it for any date you may want to explore)
final_2a_sort_BATS$wintersols[final_2a_sort_BATS$Date>as.Date("2016-01-01") & final_2a_sort_BATS$Date<as.Date("2016-12-21")] <- "2015-12-21"
final_2a_sort_BATS$wintersols[final_2a_sort_BATS$Date>as.Date("2017-01-01") & final_2a_sort_BATS$Date<as.Date("2017-12-21")] <- "2016-12-21"
final_2a_sort_BATS$wintersols[final_2a_sort_BATS$Date>as.Date("2018-01-01") & final_2a_sort_BATS$Date<as.Date("2018-12-21")] <- "2017-12-21" #The most recent date is 2018-12-15

final_2a_sort_BATS$DaysSincewintersols<-as.numeric(-difftime(final_2a_sort_BATS$wintersols,final_2a_sort_BATS$Date)+1) #Add the numeric column

#jultonumeric
final_2a_sort_BATS$Jul<-as.numeric(final_2a_sort_BATS$Jul)

#models based on the sin and cos
final_2a_sort_BATS$cosDX1<-as.numeric(cos(2*pi*((final_2a_sort_BATS$DaysSincewintersols)/365))) #Peaking in midwinter
final_2a_sort_BATS$sinDX1<-as.numeric(sin(2*pi*((final_2a_sort_BATS$Jul)/365))) #Peaking other time

#Get the linear model for all the columns
#FIRST: loop which regresses each ASV column to cosDX and save these regressions on a list (1:501)


#[502] "Date" [503] "Jul"[504]"Serial" [505]"wintersols" [506]"DaysSincewintersols" [507] "cosDX1" [508] "senDX1"

storage_BATS <- list()
for(i in names(final_2a_sort_BATS)[1:501]){
  storage_BATS[[i]] <- lm(get(i) ~ sinDX1+cosDX1, final_2a_sort_BATS)
}

#Create a data set with the regressions
fitted.values(storage_BATS[[1]])

#Get the significant regressions (THIS MUST BE THE SAME ASVS as I did with g.fisher
pvalues_lm_BATS<-sapply(storage_BATS, function(x){
  ff <- summary(x)$fstatistic
  pf(ff[1], df1 = ff[2], df2 = ff[3], lower.tail = FALSE)
})

#Before proceding some ASVs retrieve a "NA" for the sake of our pipeline I'll remove these. 
pvalues_lm_BATS<-pvalues_lm_BATS[!is.na(pvalues_lm_BATS)]
pvalues_lm_BATS05<-pvalues_lm_BATS[pvalues_lm_BATS < 0.05] #The names here should be ~ to g.fisher HERE=223/501 = 44.5% g.fisher predicted 49.7 :O !!! This is so cool!

#All the significant 
sign_names_BATS<-str_remove(names(pvalues_lm_BATS05), ".value") #extract this from both datasets

#Number of the column that matches a list of names get the position of the names <0.05
signll_BATS<-match(sign_names_BATS,names(final_2a_sort_BATS))

#subset list of lists 
storage_BATS05<-storage_BATS[signll_BATS]

#I need the constant, sin and cos of the models, the cosDX and SinDX originals and
#Get the constant, sin and cos of the models:
coeffs05_BATS<-data.frame(coeffs=sapply(storage_BATS05, FUN=function(item){item$coefficients}))
colnames(coeffs05_BATS)<-str_remove(colnames(coeffs05_BATS), "coeffs.")

#Subset original day, sin, cos
newdf_model_BATS <- final_2a_sort_BATS[, c("Serial","cosDX1","sinDX1")]

#automatizise this :P 
newdf_model_BATS$BIOSm1<-coeffs05_BATS$BIOSm1[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1[3])
newdf_model_BATS$BIOSm15<-coeffs05_BATS$BIOSm15[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm15[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm15[3])
newdf_model_BATS$BIOSm42<-coeffs05_BATS$BIOSm42[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm42[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm42[3])
newdf_model_BATS$BIOSm267<-coeffs05_BATS$BIOSm267[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm267[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm267[3])
newdf_model_BATS$BIOSm270<-coeffs05_BATS$BIOSm270[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm270[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm270[3])
newdf_model_BATS$BIOSm18435<-coeffs05_BATS$BIOSm18435[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm18435[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm18435[3])
newdf_model_BATS$BIOSm80<-coeffs05_BATS$BIOSm80[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm80[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm80[3])
newdf_model_BATS$BIOSm12643<-coeffs05_BATS$BIOSm12643[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm12643[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm12643[3])
newdf_model_BATS$BIOSm34348<-coeffs05_BATS$BIOSm34348[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm34348[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm34348[3])
newdf_model_BATS$BIOSm28127<-coeffs05_BATS$BIOSm28127[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm28127[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm28127[3])
newdf_model_BATS$BIOSm17499<-coeffs05_BATS$BIOSm17499[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17499[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17499[3])
newdf_model_BATS$BIOSm17536<-coeffs05_BATS$BIOSm17536[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17536[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17536[3])
newdf_model_BATS$BIOSm463<-coeffs05_BATS$BIOSm463[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm463[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm463[3])
newdf_model_BATS$BIOSm12676<-coeffs05_BATS$BIOSm12676[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm12676[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm12676[3])
newdf_model_BATS$BIOSm1384<-coeffs05_BATS$BIOSm1384[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1384[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1384[3])
newdf_model_BATS$BIOSm95<-coeffs05_BATS$BIOSm95[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm95[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm95[3])
newdf_model_BATS$BIOSm995<-coeffs05_BATS$BIOSm995[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm995[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm995[3])
newdf_model_BATS$BIOSm1481<-coeffs05_BATS$BIOSm1481[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1481[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1481[3])
newdf_model_BATS$BIOSm697<-coeffs05_BATS$BIOSm697[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm697[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm697[3])
newdf_model_BATS$BIOSm17576<-coeffs05_BATS$BIOSm17576[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17576[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17576[3])
newdf_model_BATS$BIOSm26903<-coeffs05_BATS$BIOSm26903[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm26903[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm26903[3])
newdf_model_BATS$BIOSm276<-coeffs05_BATS$BIOSm276[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm276[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm276[3])
newdf_model_BATS$BIOSm18390<-coeffs05_BATS$BIOSm18390[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm18390[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm18390[3])
newdf_model_BATS$BIOSm26981<-coeffs05_BATS$BIOSm26981[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm26981[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm26981[3])
newdf_model_BATS$BIOSm17568<-coeffs05_BATS$BIOSm17568[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17568[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17568[3])
newdf_model_BATS$BIOSm851<-coeffs05_BATS$BIOSm851[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm851[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm851[3])
newdf_model_BATS$BIOSm52<-coeffs05_BATS$BIOSm52[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm52[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm52[3])
newdf_model_BATS$BIOSm201<-coeffs05_BATS$BIOSm201[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm201[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm201[3])
newdf_model_BATS$BIOSm705<-coeffs05_BATS$BIOSm705[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm705[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm705[3])
newdf_model_BATS$BIOSm468<-coeffs05_BATS$BIOSm468[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm468[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm468[3])
newdf_model_BATS$BIOSm884<-coeffs05_BATS$BIOSm884[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm884[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm884[3])
newdf_model_BATS$BIOSm199<-coeffs05_BATS$BIOSm199[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm199[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm199[3])
newdf_model_BATS$BIOSm630<-coeffs05_BATS$BIOSm630[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm630[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm630[3])
newdf_model_BATS$BIOSm1052<-coeffs05_BATS$BIOSm1052[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1052[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1052[3])
newdf_model_BATS$BIOSm1749<-coeffs05_BATS$BIOSm1749[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1749[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1749[3])
newdf_model_BATS$BIOSm5<-coeffs05_BATS$BIOSm5[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm5[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm5[3])
newdf_model_BATS$BIOSm147<-coeffs05_BATS$BIOSm147[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm147[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm147[3])
newdf_model_BATS$BIOSm17788<-coeffs05_BATS$BIOSm17788[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17788[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17788[3])
newdf_model_BATS$BIOSm26834<-coeffs05_BATS$BIOSm26834[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm26834[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm26834[3])
newdf_model_BATS$BIOSm26900<-coeffs05_BATS$BIOSm26900[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm26900[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm26900[3])
newdf_model_BATS$BIOSm22<-coeffs05_BATS$BIOSm22[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm22[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm22[3])
newdf_model_BATS$BIOSm26856<-coeffs05_BATS$BIOSm26856[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm26856[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm26856[3])
newdf_model_BATS$BIOSm1221<-coeffs05_BATS$BIOSm1221[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1221[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1221[3])
newdf_model_BATS$BIOSm7<-coeffs05_BATS$BIOSm7[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm7[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm7[3])
newdf_model_BATS$BIOSm17455<-coeffs05_BATS$BIOSm17455[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17455[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17455[3])
newdf_model_BATS$BIOSm18302<-coeffs05_BATS$BIOSm18302[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm18302[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm18302[3])
newdf_model_BATS$BIOSm17532<-coeffs05_BATS$BIOSm17532[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17532[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17532[3])
newdf_model_BATS$BIOSm1006<-coeffs05_BATS$BIOSm1006[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1006[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1006[3])
newdf_model_BATS$BIOSm17414<-coeffs05_BATS$BIOSm17414[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17414[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17414[3])
newdf_model_BATS$BIOSm17531<-coeffs05_BATS$BIOSm17531[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17531[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17531[3])
newdf_model_BATS$BIOSm243<-coeffs05_BATS$BIOSm243[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm243[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm243[3])
newdf_model_BATS$BIOSm17909<-coeffs05_BATS$BIOSm17909[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17909[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17909[3])
newdf_model_BATS$BIOSm541<-coeffs05_BATS$BIOSm541[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm541[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm541[3])
newdf_model_BATS$BIOSm59<-coeffs05_BATS$BIOSm59[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm59[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm59[3])
newdf_model_BATS$BIOSm17881<-coeffs05_BATS$BIOSm17881[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17881[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17881[3])
newdf_model_BATS$BIOSm17540<-coeffs05_BATS$BIOSm17540[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17540[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17540[3])
newdf_model_BATS$BIOSm17854<-coeffs05_BATS$BIOSm17854[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17854[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17854[3])
newdf_model_BATS$BIOSm17472<-coeffs05_BATS$BIOSm17472[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17472[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17472[3])
newdf_model_BATS$BIOSm17410<-coeffs05_BATS$BIOSm17410[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17410[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17410[3])
newdf_model_BATS$BIOSm1932<-coeffs05_BATS$BIOSm1932[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1932[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1932[3])
newdf_model_BATS$BIOSm447<-coeffs05_BATS$BIOSm447[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm447[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm447[3])
newdf_model_BATS$BIOSm2481<-coeffs05_BATS$BIOSm2481[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm2481[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm2481[3])
newdf_model_BATS$BIOSm2921<-coeffs05_BATS$BIOSm2921[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm2921[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm2921[3])
newdf_model_BATS$BIOSm313<-coeffs05_BATS$BIOSm313[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm313[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm313[3])
newdf_model_BATS$BIOSm1012<-coeffs05_BATS$BIOSm1012[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1012[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1012[3])
newdf_model_BATS$BIOSm17437<-coeffs05_BATS$BIOSm17437[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17437[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17437[3])
newdf_model_BATS$BIOSm17467<-coeffs05_BATS$BIOSm17467[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17467[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17467[3])
newdf_model_BATS$BIOSm17393<-coeffs05_BATS$BIOSm17393[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17393[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17393[3])
newdf_model_BATS$BIOSm43<-coeffs05_BATS$BIOSm43[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm43[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm43[3])
newdf_model_BATS$BIOSm17684<-coeffs05_BATS$BIOSm17684[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17684[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17684[3])
newdf_model_BATS$BIOSm12582<-coeffs05_BATS$BIOSm12582[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm12582[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm12582[3])
newdf_model_BATS$BIOSm17415<-coeffs05_BATS$BIOSm17415[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17415[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17415[3])
newdf_model_BATS$BIOSm19093<-coeffs05_BATS$BIOSm19093[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm19093[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm19093[3])
newdf_model_BATS$BIOSm1051<-coeffs05_BATS$BIOSm1051[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1051[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1051[3])
newdf_model_BATS$BIOSm17589<-coeffs05_BATS$BIOSm17589[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17589[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17589[3])
newdf_model_BATS$BIOSm17781<-coeffs05_BATS$BIOSm17781[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17781[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17781[3])
newdf_model_BATS$BIOSm18908<-coeffs05_BATS$BIOSm18908[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm18908[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm18908[3])
newdf_model_BATS$BIOSm17407<-coeffs05_BATS$BIOSm17407[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17407[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17407[3])
newdf_model_BATS$BIOSm17492<-coeffs05_BATS$BIOSm17492[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17492[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17492[3])
newdf_model_BATS$BIOSm17689<-coeffs05_BATS$BIOSm17689[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17689[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17689[3])
newdf_model_BATS$BIOSm17465<-coeffs05_BATS$BIOSm17465[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17465[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17465[3])
newdf_model_BATS$BIOSm17424<-coeffs05_BATS$BIOSm17424[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17424[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17424[3])
newdf_model_BATS$BIOSm17409<-coeffs05_BATS$BIOSm17409[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17409[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17409[3])
newdf_model_BATS$BIOSm19057<-coeffs05_BATS$BIOSm19057[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm19057[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm19057[3])
newdf_model_BATS$BIOSm17986<-coeffs05_BATS$BIOSm17986[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17986[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17986[3])
newdf_model_BATS$BIOSm402<-coeffs05_BATS$BIOSm402[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm402[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm402[3])
newdf_model_BATS$BIOSm169<-coeffs05_BATS$BIOSm169[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm169[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm169[3])
newdf_model_BATS$BIOSm423<-coeffs05_BATS$BIOSm423[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm423[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm423[3])
newdf_model_BATS$BIOSm12612<-coeffs05_BATS$BIOSm12612[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm12612[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm12612[3])
newdf_model_BATS$BIOSm122<-coeffs05_BATS$BIOSm122[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm122[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm122[3])
newdf_model_BATS$BIOSm299<-coeffs05_BATS$BIOSm299[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm299[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm299[3])
newdf_model_BATS$BIOSm12647<-coeffs05_BATS$BIOSm12647[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm12647[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm12647[3])
newdf_model_BATS$BIOSm17463<-coeffs05_BATS$BIOSm17463[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17463[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17463[3])
newdf_model_BATS$BIOSm17770<-coeffs05_BATS$BIOSm17770[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17770[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17770[3])
newdf_model_BATS$BIOSm17486<-coeffs05_BATS$BIOSm17486[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17486[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17486[3])
newdf_model_BATS$BIOSm17774<-coeffs05_BATS$BIOSm17774[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17774[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17774[3])
newdf_model_BATS$BIOSm817<-coeffs05_BATS$BIOSm817[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm817[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm817[3])
newdf_model_BATS$BIOSm1205<-coeffs05_BATS$BIOSm1205[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1205[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1205[3])
newdf_model_BATS$BIOSm1055<-coeffs05_BATS$BIOSm1055[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1055[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1055[3])
newdf_model_BATS$BIOSm1914<-coeffs05_BATS$BIOSm1914[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1914[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1914[3])
newdf_model_BATS$BIOSm17930<-coeffs05_BATS$BIOSm17930[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17930[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17930[3])
newdf_model_BATS$BIOSm1772<-coeffs05_BATS$BIOSm1772[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1772[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1772[3])
newdf_model_BATS$BIOSm577<-coeffs05_BATS$BIOSm577[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm577[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm577[3])
newdf_model_BATS$BIOSm12646<-coeffs05_BATS$BIOSm12646[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm12646[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm12646[3])
newdf_model_BATS$BIOSm18032<-coeffs05_BATS$BIOSm18032[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm18032[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm18032[3])
newdf_model_BATS$BIOSm18902<-coeffs05_BATS$BIOSm18902[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm18902[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm18902[3])
newdf_model_BATS$BIOSm75<-coeffs05_BATS$BIOSm75[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm75[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm75[3])
newdf_model_BATS$BIOSm581<-coeffs05_BATS$BIOSm581[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm581[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm581[3])
newdf_model_BATS$BIOSm1028<-coeffs05_BATS$BIOSm1028[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1028[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1028[3])
newdf_model_BATS$BIOSm1683<-coeffs05_BATS$BIOSm1683[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1683[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1683[3])
newdf_model_BATS$BIOSm17706<-coeffs05_BATS$BIOSm17706[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17706[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17706[3])
newdf_model_BATS$BIOSm253<-coeffs05_BATS$BIOSm253[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm253[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm253[3])
newdf_model_BATS$BIOSm841<-coeffs05_BATS$BIOSm841[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm841[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm841[3])
newdf_model_BATS$BIOSm1394<-coeffs05_BATS$BIOSm1394[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1394[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1394[3])
newdf_model_BATS$BIOSm607<-coeffs05_BATS$BIOSm607[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm607[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm607[3])
newdf_model_BATS$BIOSm17487<-coeffs05_BATS$BIOSm17487[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17487[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17487[3])
newdf_model_BATS$BIOSm982<-coeffs05_BATS$BIOSm982[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm982[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm982[3])
newdf_model_BATS$BIOSm275<-coeffs05_BATS$BIOSm275[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm275[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm275[3])
newdf_model_BATS$BIOSm316<-coeffs05_BATS$BIOSm316[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm316[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm316[3])
newdf_model_BATS$BIOSm414<-coeffs05_BATS$BIOSm414[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm414[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm414[3])
newdf_model_BATS$BIOSm298<-coeffs05_BATS$BIOSm298[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm298[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm298[3])
newdf_model_BATS$BIOSm17413<-coeffs05_BATS$BIOSm17413[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17413[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17413[3])
newdf_model_BATS$BIOSm517<-coeffs05_BATS$BIOSm517[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm517[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm517[3])
newdf_model_BATS$BIOSm632<-coeffs05_BATS$BIOSm632[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm632[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm632[3])
newdf_model_BATS$BIOSm19393<-coeffs05_BATS$BIOSm19393[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm19393[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm19393[3])
newdf_model_BATS$BIOSm17598<-coeffs05_BATS$BIOSm17598[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17598[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17598[3])
newdf_model_BATS$BIOSm17928<-coeffs05_BATS$BIOSm17928[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17928[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17928[3])
newdf_model_BATS$BIOSm234<-coeffs05_BATS$BIOSm234[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm234[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm234[3])
newdf_model_BATS$BIOSm17968<-coeffs05_BATS$BIOSm17968[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17968[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17968[3])
newdf_model_BATS$BIOSm18883<-coeffs05_BATS$BIOSm18883[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm18883[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm18883[3])
newdf_model_BATS$BIOSm25<-coeffs05_BATS$BIOSm25[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm25[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm25[3])
newdf_model_BATS$BIOSm18348<-coeffs05_BATS$BIOSm18348[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm18348[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm18348[3])
newdf_model_BATS$BIOSm242<-coeffs05_BATS$BIOSm242[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm242[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm242[3])
newdf_model_BATS$BIOSm27972<-coeffs05_BATS$BIOSm27972[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm27972[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm27972[3])
newdf_model_BATS$BIOSm716<-coeffs05_BATS$BIOSm716[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm716[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm716[3])
newdf_model_BATS$BIOSm1459<-coeffs05_BATS$BIOSm1459[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1459[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1459[3])
newdf_model_BATS$BIOSm1708<-coeffs05_BATS$BIOSm1708[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1708[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1708[3])
newdf_model_BATS$BIOSm18411<-coeffs05_BATS$BIOSm18411[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm18411[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm18411[3])
newdf_model_BATS$BIOSm2653<-coeffs05_BATS$BIOSm2653[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm2653[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm2653[3])
newdf_model_BATS$BIOSm723<-coeffs05_BATS$BIOSm723[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm723[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm723[3])
newdf_model_BATS$BIOSm18163<-coeffs05_BATS$BIOSm18163[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm18163[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm18163[3])
newdf_model_BATS$BIOSm466<-coeffs05_BATS$BIOSm466[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm466[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm466[3])
newdf_model_BATS$BIOSm17571<-coeffs05_BATS$BIOSm17571[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17571[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17571[3])
newdf_model_BATS$BIOSm18008<-coeffs05_BATS$BIOSm18008[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm18008[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm18008[3])
newdf_model_BATS$BIOSm12799<-coeffs05_BATS$BIOSm12799[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm12799[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm12799[3])
newdf_model_BATS$BIOSm565<-coeffs05_BATS$BIOSm565[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm565[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm565[3])
newdf_model_BATS$BIOSm520<-coeffs05_BATS$BIOSm520[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm520[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm520[3])
newdf_model_BATS$BIOSm585<-coeffs05_BATS$BIOSm585[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm585[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm585[3])
newdf_model_BATS$BIOSm377<-coeffs05_BATS$BIOSm377[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm377[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm377[3])
newdf_model_BATS$BIOSm36<-coeffs05_BATS$BIOSm36[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm36[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm36[3])
newdf_model_BATS$BIOSm162<-coeffs05_BATS$BIOSm162[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm162[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm162[3])
newdf_model_BATS$BIOSm1595<-coeffs05_BATS$BIOSm1595[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1595[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1595[3])
newdf_model_BATS$P708<-coeffs05_BATS$P708[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$P708[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$P708[3])
newdf_model_BATS$BIOSm600<-coeffs05_BATS$BIOSm600[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm600[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm600[3])
newdf_model_BATS$BIOSm750<-coeffs05_BATS$BIOSm750[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm750[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm750[3])
newdf_model_BATS$BIOSm118<-coeffs05_BATS$BIOSm118[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm118[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm118[3])
newdf_model_BATS$BIOSm126<-coeffs05_BATS$BIOSm126[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm126[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm126[3])
newdf_model_BATS$BIOSm168<-coeffs05_BATS$BIOSm168[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm168[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm168[3])
newdf_model_BATS$BIOSm219<-coeffs05_BATS$BIOSm219[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm219[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm219[3])
newdf_model_BATS$BIOSm312<-coeffs05_BATS$BIOSm312[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm312[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm312[3])
newdf_model_BATS$BIOSm443<-coeffs05_BATS$BIOSm443[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm443[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm443[3])
newdf_model_BATS$BIOSm347<-coeffs05_BATS$BIOSm347[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm347[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm347[3])
newdf_model_BATS$BIOSm334<-coeffs05_BATS$BIOSm334[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm334[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm334[3])
newdf_model_BATS$BIOSm704<-coeffs05_BATS$BIOSm704[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm704[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm704[3])
newdf_model_BATS$BIOSm208<-coeffs05_BATS$BIOSm208[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm208[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm208[3])
newdf_model_BATS$BIOSm209<-coeffs05_BATS$BIOSm209[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm209[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm209[3])
newdf_model_BATS$BIOSm462<-coeffs05_BATS$BIOSm462[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm462[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm462[3])
newdf_model_BATS$BIOSm56<-coeffs05_BATS$BIOSm56[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm56[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm56[3])
newdf_model_BATS$BIOSm102<-coeffs05_BATS$BIOSm102[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm102[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm102[3])
newdf_model_BATS$BIOSm973<-coeffs05_BATS$BIOSm973[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm973[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm973[3])
newdf_model_BATS$BIOSm575<-coeffs05_BATS$BIOSm575[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm575[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm575[3])
newdf_model_BATS$BIOSm71<-coeffs05_BATS$BIOSm71[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm71[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm71[3])
newdf_model_BATS$BIOSm164<-coeffs05_BATS$BIOSm164[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm164[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm164[3])
newdf_model_BATS$BIOSm17580<-coeffs05_BATS$BIOSm17580[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17580[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17580[3])
newdf_model_BATS$BIOSm13042<-coeffs05_BATS$BIOSm13042[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm13042[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm13042[3])
newdf_model_BATS$BIOSm678<-coeffs05_BATS$BIOSm678[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm678[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm678[3])
newdf_model_BATS$BIOSm1265<-coeffs05_BATS$BIOSm1265[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1265[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1265[3])
newdf_model_BATS$BIOSm1510<-coeffs05_BATS$BIOSm1510[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1510[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1510[3])
newdf_model_BATS$BIOSm17569<-coeffs05_BATS$BIOSm17569[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17569[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17569[3])
newdf_model_BATS$BIOSm806<-coeffs05_BATS$BIOSm806[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm806[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm806[3])
newdf_model_BATS$BIOSm108<-coeffs05_BATS$BIOSm108[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm108[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm108[3])
newdf_model_BATS$BIOSm393<-coeffs05_BATS$BIOSm393[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm393[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm393[3])
newdf_model_BATS$BIOSm1277<-coeffs05_BATS$BIOSm1277[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1277[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1277[3])
newdf_model_BATS$BIOSm120<-coeffs05_BATS$BIOSm120[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm120[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm120[3])
newdf_model_BATS$BIOSm302<-coeffs05_BATS$BIOSm302[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm302[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm302[3])
newdf_model_BATS$BIOSm307<-coeffs05_BATS$BIOSm307[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm307[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm307[3])
newdf_model_BATS$BIOSm2455<-coeffs05_BATS$BIOSm2455[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm2455[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm2455[3])
newdf_model_BATS$BIOSm172<-coeffs05_BATS$BIOSm172[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm172[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm172[3])
newdf_model_BATS$BIOSm193<-coeffs05_BATS$BIOSm193[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm193[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm193[3])
newdf_model_BATS$BIOSm1001<-coeffs05_BATS$BIOSm1001[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1001[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1001[3])
newdf_model_BATS$BIOSm124<-coeffs05_BATS$BIOSm124[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm124[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm124[3])
newdf_model_BATS$BIOSm1181<-coeffs05_BATS$BIOSm1181[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1181[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1181[3])
newdf_model_BATS$BIOSm408<-coeffs05_BATS$BIOSm408[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm408[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm408[3])
newdf_model_BATS$BIOSm1323<-coeffs05_BATS$BIOSm1323[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1323[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1323[3])
newdf_model_BATS$BIOSm634<-coeffs05_BATS$BIOSm634[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm634[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm634[3])
newdf_model_BATS$BIOSm2253<-coeffs05_BATS$BIOSm2253[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm2253[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm2253[3])
newdf_model_BATS$BIOSm1833<-coeffs05_BATS$BIOSm1833[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1833[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1833[3])
newdf_model_BATS$BIOSm946<-coeffs05_BATS$BIOSm946[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm946[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm946[3])
newdf_model_BATS$BIOSm245<-coeffs05_BATS$BIOSm245[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm245[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm245[3])
newdf_model_BATS$BIOSm2159<-coeffs05_BATS$BIOSm2159[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm2159[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm2159[3])
newdf_model_BATS$BIOSm308<-coeffs05_BATS$BIOSm308[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm308[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm308[3])
newdf_model_BATS$BIOSm420<-coeffs05_BATS$BIOSm420[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm420[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm420[3])
newdf_model_BATS$BIOSm523<-coeffs05_BATS$BIOSm523[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm523[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm523[3])
newdf_model_BATS$BIOSm930<-coeffs05_BATS$BIOSm930[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm930[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm930[3])
newdf_model_BATS$BIOSm4927<-coeffs05_BATS$BIOSm4927[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm4927[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm4927[3])
newdf_model_BATS$BIOSm871<-coeffs05_BATS$BIOSm871[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm871[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm871[3])
newdf_model_BATS$BIOSm872<-coeffs05_BATS$BIOSm872[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm872[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm872[3])
newdf_model_BATS$BIOSm160<-coeffs05_BATS$BIOSm160[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm160[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm160[3])
newdf_model_BATS$BIOSm1361<-coeffs05_BATS$BIOSm1361[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1361[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1361[3])
newdf_model_BATS$BIOSm1298<-coeffs05_BATS$BIOSm1298[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1298[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1298[3])
newdf_model_BATS$BIOSm329<-coeffs05_BATS$BIOSm329[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm329[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm329[3])
newdf_model_BATS$BIOSm1792<-coeffs05_BATS$BIOSm1792[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1792[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1792[3])
newdf_model_BATS$BIOSm2772<-coeffs05_BATS$BIOSm2772[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm2772[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm2772[3])
newdf_model_BATS$BIOSm1048<-coeffs05_BATS$BIOSm1048[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1048[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1048[3])
newdf_model_BATS$BIOSm39<-coeffs05_BATS$BIOSm39[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm39[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm39[3])
newdf_model_BATS$BIOSm17538<-coeffs05_BATS$BIOSm17538[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17538[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17538[3])
newdf_model_BATS$BIOSm99<-coeffs05_BATS$BIOSm99[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm99[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm99[3])
newdf_model_BATS$BIOSm112<-coeffs05_BATS$BIOSm112[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm112[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm112[3])
newdf_model_BATS$BIOSm272<-coeffs05_BATS$BIOSm272[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm272[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm272[3])
newdf_model_BATS$BIOSm519<-coeffs05_BATS$BIOSm519[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm519[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm519[3])
newdf_model_BATS$BIOSm17562<-coeffs05_BATS$BIOSm17562[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm17562[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm17562[3])
newdf_model_BATS$BIOSm406<-coeffs05_BATS$BIOSm406[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm406[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm406[3])
newdf_model_BATS$BIOSm1375<-coeffs05_BATS$BIOSm1375[1]+(newdf_model_BATS$sinDX1*coeffs05_BATS$BIOSm1375[2])+(newdf_model_BATS$cosDX1*coeffs05_BATS$BIOSm1375[3])


#newdf_model_BATS for melt 

Tomelt_BATS<-(newdf_model_BATS[c(1,4:226)])

#Add taxa to this melted 
newdf_model_BATS_melted<-melt(Tomelt_BATS,id.vars="Serial")

#Agregar taxa extracting sign_names variable from the Phyloseq object, thereafter extract the tax table only with genus 

##### Subset from physeq sar11.
phy_sign05_BATS = prune_taxa(sign_names_BATS, SAR11B5)
GenusSAR11sign1_BATS<-data.frame(tax_table(phy_sign05_BATS))[,8,drop=FALSE]

GenusSAR11sign1_BATS$ASV<-rownames(GenusSAR11sign1_BATS) #GenusSAR11sign1_BATS contains the neccesary to add taxa to melted

#change variable to ASV in newdf_model_BATS_melted 
names(newdf_model_BATS_melted)[2]<-"ASV"

#merge by ASV and add the column Genus 
newdf_model_BATS_melted_taxa<-merge(newdf_model_BATS_melted, GenusSAR11sign1_BATS[, c("tax", "ASV")], by="ASV")

#Subset 1 to 365 
newdf_model_BATS_melted_taxa365<-newdf_model_BATS_melted_taxa[newdf_model_BATS_melted_taxa$Serial >730 & newdf_model_BATS_melted_taxa$Serial <1095,] #2nd year
newdf_model_BATS_melted_taxa365$Serial365<-(newdf_model_BATS_melted_taxa365$Serial-730)

#unique(newdf_model_BATS_melted_taxa365$tax) => "I;Ia,Ia.3","I","I;Ia,Ia.4","I;Ib;Ib.1","I;Ib","I;Ib;Ib.2","I;Ic;Ic.2","I;Ia,Ia.1"

#Colors for SAR11 subclasses
coloresbarplot_BATS = c("I;Ia,Ia.1"="darkcyan","I;Ia,Ia.3"="blue","I"="black","I;Ia,Ia.4"="lightskyblue","I;Ib;Ib.1"="darkred","I;Ib;Ib.2"="firebrick1","II;IIa;IIa.A"="bisque4","II;IIa;IIa.B"="darkgoldenrod3 ","II;IIa;IIa_NA"="gold2","III;IIIa"="cornsilk3","II"="chocolate4","II;IIb"="olivedrab","IV"="darkgreen","I;Ic;Ic.2"="lightcoral","I;Ic;Ic.1"="hotpink","clade N2"="lavenderblush2", "I;Ib"="#d45d41")

#I will split the panels in top and bottom as the heat map. The only difference is that I´ll use 10 as value to divide. The 4 most abundant are in a larger scale

d_BATS<-colnames(SAR11_BATS_ASVt)[apply(SAR11_BATS_ASVt, 2, function(u) colSums(SAR11_BATS_ASVt)>5)] #Get the names of most abundant
d_BATS <- d_BATS[!is.na(d_BATS)] 
a_BATS<-colnames(SAR11_BATS_ASVt)[apply(SAR11_BATS_ASVt, 2, function(u) colSums(SAR11_BATS_ASVt)<=5) & colSums(SAR11_BATS_ASVt)>0.5] #Get the names of least abundant
a_BATS <- a_BATS[!is.na(a_BATS)]
t_BATS<-colnames(SAR11_BATS_ASVt)[apply(SAR11_BATS_ASVt, 2, function(u) colSums(SAR11_BATS_ASVt)<=.05)] #Get the names of least abundant
t_BATS<- t_BATS[!is.na(t_BATS)]

###Subset from the melted object

highSAR11_BATS<- newdf_model_BATS_melted_taxa365[newdf_model_BATS_melted_taxa365$ASV %in% d_BATS, ]
middleSAR11_BATS<- newdf_model_BATS_melted_taxa365[newdf_model_BATS_melted_taxa365$ASV %in% a_BATS, ]
lowSAR11_BATS<-newdf_model_BATS_melted_taxa365[newdf_model_BATS_melted_taxa365$ASV %in% t_BATS, ]

#Plot BATS
myplothigh_BATS<-ggplot(highSAR11_BATS, aes(x = Serial365, y =  value, colour=tax, group = ASV, label=ASV)) +geom_line(size=.8) + 
  directlabels::geom_dl(aes(label = ASV),  method = list(cex = .9, "smart.grid"))+ scale_color_manual(values = coloresbarplot_BATS)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=14,colour = "black"), axis.text.x=element_text(size=14,colour = "black"),legend.text=element_text(size=12))+labs(x = "day", y = "relative contribution")+ guides(color = guide_legend(override.aes = list(size = 3.5)))

myplotmid_BATS<-ggplot(middleSAR11_BATS, aes(x = Serial365, y =  value, colour=tax, group = ASV, label=ASV)) + geom_line(size=.8) + scale_color_manual(values = coloresbarplot_BATS)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=14,colour = "black"), axis.text.x=element_text(size=14,colour = "black"),legend.text=element_text(size=12))+labs(x = "day", y = "relative contribution")+ guides(color = guide_legend(override.aes = list(size = 3.5)))

myplotlow_BATS<-ggplot(lowSAR11_BATS, aes(x = Serial365, y =  value, colour=tax, group = ASV, label=ASV)) + geom_line(size=.8) + scale_color_manual(values = coloresbarplot_BATS)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_text(size=14,colour = "black"),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=14,colour = "black"), axis.text.x=element_text(size=14,colour = "black"),legend.text=element_text(size=12))+labs(x = "day", y = "relative contribution")+ guides(color = guide_legend(override.aes = list(size = 3.5)))

meds_BATS<-plot_grid(myplothigh_BATS+theme(legend.position="none"),myplotmid_BATS+theme(legend.position="none"),myplotlow_BATS+theme(legend.position="none"), ncol=1, align="hv", axis = "bl")

#Get legend
legend_BATS <- get_legend(
  # create some space to the left of the legend
  myplotlow_BATS + theme(legend.box.margin = margin(0, 0, 0, 12))
)

BATS_model_plot<-plot_grid(meds_BATS, legend_BATS, rel_widths = c(2.2, .5))


#####FINAL PLOT!####

plot_grid(BATS_model_plot,WEC_model_plot, ncol=2)

