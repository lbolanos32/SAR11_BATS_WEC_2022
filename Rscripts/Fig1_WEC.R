#This is the code to generate the Figure 1 of the paper and the WEC percentages used in the main text. To reproduce, use the provided tables to build the phyloseq object

#Figure 1: Monthly relative contribution of surface SAR11 ecotypes through a multiannual 16S rRNA survey at the station L4 from the Western English Channel. 
#Panel (a): Monthly relative contribution of the SAR11 fraction to the total of the amplicon dataset. The black line represents temperature through the sampled years. 
#Panel (b): Monthly relative contribution of ecotypes to the total SAR11 fraction.

###NOTES###

#NOTE1: Taxa I,I_basal has been relabel to I 

#Load libraries
library("phyloseq")
library("ggplot2")
library("dplyr")
library("tidyr")
library("reshape")
library("see")
library("cowplot")
library(grid)
library(gridExtra)

#create a function to to reformat the labels on the x axis of the time series (One letter code for months)

dte_formatter <- function(x) { 
  #formatter for axis labels: J12, F12, M12, etc... 
  mth <- substr(format(x, "%b"),1,1) 
  mth 
}


#Create the phyloseq object

countwec <- read.table("WEC_mar.otu", header=T, row.names=1, check.names=F)
sampleinfowec <- read.table("WEC_mar.env", header=T, row.names=1, check.names=F, sep ="\t")
taxwec <- as.matrix(read.table("WEC_marV2.tax", header=T, row.names=1, check.names=F, na.strings="", sep="\t")) #New edition with consistent Family and V2 have 1b

OTUmar = otu_table(countwec, taxa_are_rows = TRUE)
TAXmar = tax_table(taxwec)
SAMar= sample_data(sampleinfowec)

WEC<-phyloseq(OTUmar,TAXmar,SAMar)
WECcf= filter_taxa(WEC, function(x) sum(x > 1) > (0.0015*length(x)), TRUE)

#Determine "others"#

family_counts_tab <- otu_table(tax_glom(WECcf, taxrank="Family")) 

#making a vector of phyla names to set as row names

family_tax_vec <- as.vector(tax_table(tax_glom(WECcf, taxrank="Family"))[,4]) 
rownames(family_counts_tab) <- as.vector(family_tax_vec)


temp_major_taxa_counts_tab <- family_counts_tab[!row.names(family_counts_tab) %in% "SAR11_clade", ] #sacamos SAR11 y sumamos estos van a ser el rest

rest<-colSums(temp_major_taxa_counts_tab)
SAR11<-colSums(family_counts_tab)-rest
SF<-rbind("SAR11"=SAR11, "WEC_Other"=rest)

#######New OTU
SAR11_WECphy = subset_taxa(WECcf, Order=="SAR11_clade")
SAR11_ASVs_tab <- otu_table(SAR11_WECphy) 
SAR11_rest_table<-as.matrix(rbind.data.frame(SAR11_ASVs_tab, "WEC99999"=rest))
SAR11_rest_otutable<-otu_table(SAR11_rest_table, taxa_are_rows=TRUE)


######New Tax
SAR11_tax_tab <- tax_table(SAR11_WECphy) 
Othervec<-c("other","other","other","other","other","other")
SAR11_restTAX_table<-tax_table(as.matrix(rbind(SAR11_tax_tab, "WEC99999"=Othervec)))

#####New phyloseq
#SAMar set different columns for Month, Day, Year and DATE
sampleinfowec$Date<-as.Date(sampleinfowec$Date,"%m/%d/%Y") 
sampleinfowec$year = as.character(format(sampleinfowec$Date, format = "%Y"))
sampleinfowec$month = as.character(format(sampleinfowec$Date, format = "%m"))
sampleinfowec$day = as.character(format(sampleinfowec$Date, format = "%d"))
SAMar= sample_data(sampleinfowec)

SAR11_WECphy2<-phyloseq(SAR11_rest_otutable,SAR11_restTAX_table,SAMar)

#######Panel (A) SAR11 relative to the total dataset#######

glomV1filt_all<-tax_glom(SAR11_WECphy2, taxrank="Order")
WEC_1aRel<-transform_sample_counts(glomV1filt_all, function(x){x / sum(x)})

tax1a<-as.data.frame(tax_table(WEC_1aRel)[,4])

a1_ASV<-as.data.frame(otu_table(WEC_1aRel))
a1_ASV[ "Taxa" ] <- tax1a[,1]
#dim(a1_ASV)
#2 265
ASV_frame1a <- a1_ASV[,-265]
rownames(ASV_frame1a) <- a1_ASV[,265]
ASV_frw1a<-t(ASV_frame1a)
md_to_add1a<-as.data.frame(sample_data(WEC_1aRel))[,c(1,4,72,116,117)]

final_1aM<-cbind(ASV_frw1a,md_to_add1a)

ASVsToMelt1a<-colnames(final_1aM)[c(1:2)]
final_s11_melt<-melt(final_1aM,id.vars=c("ID","Date","temp_C_2.5","Chl_0m","Chl_10m"), measure.vars = c("SAR11_clade","other"))

final_s11_melt$variable <- factor(final_s11_melt$variable, levels=c("other","SAR11_clade"))

Colors1a<- c("SAR11_clade"="#2a4858", "other"="#fbfbf6")

final_s11_melt$value <- final_s11_melt$value*100

bar1a<-ggplot(final_s11_melt, aes(x = Date, y = value)) + scale_x_date(date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+ 
geom_area(aes(color = variable, fill = variable),alpha = 0.55) +scale_color_manual(values = Colors1a) +scale_fill_manual(values =Colors1a)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=12,colour = "black"), axis.text.x=element_text(size=8,colour = "black"),legend.text=element_text(size=11))+
      labs(title = "Surface SAR11 annual contribution - 2011-2018",
           y = "Relative contribution (%)")+geom_vline(xintercept=as.Date(c("2012-01-01","2013-01-01","2014-01-01","2015-01-01","2016-01-01","2017-01-01","2018-01-01")),linetype=4, colour="black")

tmp<-ggplot(final_s11_melt, aes(x = Date, y = temp_C_2.5)) + geom_line(aes(group=1))+ scale_x_date(date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+
theme(panel.background = element_blank(),legend.position = "none",strip.background = element_blank(),strip.text.x =element_text(size=18),axis.text.y=element_text(colour = 'black',size=12), axis.text.x=element_text(colour = 'black',size=8),axis.title.y=element_text(size=14)) +ylab("Temperature [°C]")+xlab(NULL)+ylim(5,25)


#######Panel (B) SAR11 subclades relative to the total SAR11#######

#subclades 
WEC_SAR11 = subset_taxa(SAR11_WECphy2, Order=="SAR11_clade") #exctract first SAR11 
WEC_SAR11rel<-transform_sample_counts(WEC_SAR11, function(x){x / sum(x)}) #Transform to % only based on the contribution to total sar11

glomV1filt<-tax_glom(WEC_SAR11rel, taxrank="Genus")


SAR11tax<-as.data.frame(tax_table(glomV1filt)[,6])

SAR11_ASV<-as.data.frame(otu_table(glomV1filt))
SAR11_ASV[ "Taxa" ] <- SAR11tax[,1]
#dim(SAR11_ASV)
#[1]  18 265
ASV_frame2 <- SAR11_ASV[,-265]
rownames(ASV_frame2) <- SAR11_ASV[,265]
ASV_frw<-t(ASV_frame2)
md_to_add<-as.data.frame(sample_data(glomV1filt))[,c(1,4,116,117)] #Date and chlorophyll 0m and 10m 
final_2a<-cbind(ASV_frw,md_to_add)

fn2a_melt<-melt(final_2a,id.vars=c("ID","Date","Chl_0m","Chl_10m"), measure.vars = c("I","I;Ia,Ia.1","I;Ia,Ia.3","I;Ia,Ia.4","I;Ib","I;Ic;Ic.1","I;Ic;Ic.2","II","II;IIa;IIa.A","II;IIa;IIa.B","II;IIa;IIa_NA","II;IIb","III;IIIa","clade N2","IV"))

#Colors for SAR11 subclasses
coloresbarplot = c("I;Ia,Ia.1"="darkcyan","I;Ia,Ia.3"="blue","I"="black","I;Ia,Ia.4"="lightskyblue","II;IIa;IIa.A"="bisque4","II;IIa;IIa.B"="darkgoldenrod3","II;IIa;IIa_NA"="gold2","III;IIIa"="cornsilk2","II"="chocolate4","II;IIb"="olivedrab","IV"="darkgreen","I;Ic;Ic.2"="lightcoral","I;Ic;Ic.1"="hotpink","clade N2"="lavenderblush2", "I;Ib"="#d45d41")

fn2a_melt$value<-fn2a_melt$value*100

subclade_barplo<-ggplot(fn2a_melt, aes(x = Date, y = value)) + scale_x_date(date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+ 
geom_area(aes(color = variable, fill = variable),alpha = 0.7) +scale_color_manual(values = coloresbarplot) +scale_fill_manual(values =coloresbarplot)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=12,colour = "black"), axis.text.x=element_text(size=8,colour = "black"),legend.text=element_text(size=11))+
      labs(title = "Surface phylotype annual contribution to SAR11 - 2011-2018",
           y = "Relative contribution (%)")+geom_vline(xintercept=as.Date(c("2012-01-01","2013-01-01","2014-01-01","2015-01-01","2016-01-01","2017-01-01","2018-01-01")),linetype=4, colour="black")


######FINAL#####
plot_grid(bar1a,tmp,subclade_barplo, align = "v",axis="lr", nrow = 3)

######NUMBERS#####

###The following section contains the code use to retrieve the percentages or other numbers from the quality control (working) dataset

####Numbers extracted for the paper:

#final_s11_melt is divided in SAR11 %s and others %. Extract SAR11 and then

#final_s11_melt is divided in SAR11 %s and others %. Extract SAR11 and then 

sar11_checktext <- subset(final_s11_melt, variable == "SAR11_clade")
min(sar11_checktext$value)	#0.04415921
max(sar11_checktext$value)	#53.96146

#distribution

sd(sar11_checktext$value)	#9.7
mean(sar11_checktext$value)	#15.69374

###Ia ecotypes -> SAR11 Ia ecotypes dominate the surface, composing X% of the total BATS (Fig 1b). 
#SUM all Ia subclades
#Extract c("I;Ia,Ia.1","I;Ia,Ia.3","I;Ia,Ia.4")
sum(colSums(otu_table(subset_taxa(WEC_SAR11, Genus=="I;Ia,Ia.3")))) #SUM= 680548
sum(colSums(otu_table(subset_taxa(WEC_SAR11, Genus=="I;Ia,Ia.1")))) #SUM= 682489
sum(colSums(otu_table(subset_taxa(WEC_SAR11, Genus=="I;Ia,Ia.4")))) #SUM= 3107

#Sum of all SAR11 reads 1366144. 
sum(colSums(otu_table(WEC_SAR11 ))) #SUM= 1901052
#1366144/1901052 = 0.7186253

