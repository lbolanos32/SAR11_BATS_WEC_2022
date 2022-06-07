#Figure 2

#This is the code to generate the Figure 2 of the paper and the BATS percentages used in the main text. To reproduce, use the provided tables to build the phyloseq object

#Figure 2: Monthly relative contribution of surface SAR11 ecotypes through multiannual 16S rRNA surveys at the Bermuda Atlantic Time-series. 
#Panel (a) Monthly relative contribution of the SAR11 fraction to the total of the amplicon dataset. The black line represents temperature through the sampled years. 
#Panel (b) Monthly relative contribution of ecotypes to the total SAR11 fraction. Samples from 1991- 1994 and 1996-2002 were sequenced using 454 FLX technology, while 2016-2018 were sequenced using the Illumina MiSeq platform. 

###NOTES###

#NOTE1: I parsed the taxa file we have been using to have the same collapsing as BATS. 
#tail -n +2 BATSBIOS_SAR11aln_mafft_trim_hd.unique.fa.aln.jplace.tab | sort -k 1,1 -n | cut -f 1,3 > 2021_tax.txt. Then I added this as the last column in the tax file as taxa

#NOTE2: To be consistent in time and station, only BATS samples were taken except for the (July) C7_N1_S112 AE1614 that was retrieved from BIOS-SCOPE collection

#NOTE3: V1V2 can define Ib.1 and Ib.2 (BATS). While V4V5 only define up to Ib (BATS)


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

#Set the time-frames to arrange based on the discontinuity of the time-series

start.end1 <- c(as.Date("1991-08-01"),as.Date("1994-02-01"))
start.end2 <- c(as.Date("1997-09-01"),as.Date("2004-01-01"))
start.end3 <- c(as.Date("2016-03-08"),as.Date("2018-12-15"))

#Create the phyloseq object

count_tabS11phyTS5m <- read.table("BIOSotu5.otu", header=T, row.names=1, check.names=F)
sample_info_tabS11phyTS5m <- read.table("BIOSenvts5.env", header=T, row.names=1, check.names=F, sep ="\t")
tax_tabS11phyTS5m <- as.matrix(read.table("BIOStax5.tax", header=T, row.names=1, check.names=F, sep="\t",na.strings = "#NA"))
sample_info_tabS11$Date <- as.Date(with(sample_info_tabS11, paste(year, month, day,sep="-")), "%Y-%m-%d")

#Use as date
sample_info_tabS11phyTS5m$Date <- as.Date(with(sample_info_tabS11phyTS5m, paste(year, month, day,sep="-")), "%Y-%m-%d")

S11OTUphyTS5m = otu_table(count_tabS11phyTS5m, taxa_are_rows = TRUE)
S11TAXphyTS5m = tax_table(tax_tabS11phyTS5m)
S11SAMphyTS5m= sample_data(sample_info_tabS11phyTS5m)

phyTS5m<-phyloseq(S11OTUphyTS5m,S11TAXphyTS5m,S11SAMphyTS5m)

phyTS5mV1 <- prune_samples(sample_sums(phyTS5m)>=1000, phyTS5m) #prune samples less than 1000 reads (eliminate 17 --> 146/163)

phyTS5mpruneV2<-filter_taxa(phyTS5mV1, function(x) sum(x >= 3) > (0.02 *length(x)), TRUE) #V1=916/2875 Again from phyTS5mprune, V2 (146 samples, .01) ~500 ASVs

#######Panel A SAR11 relative to the total dataset#######

glomV1filt_all<-tax_glom(phyTS5mpruneV2, taxrank="Order")
BATS_1aRel<-transform_sample_counts(glomV1filt_all, function(x){x / sum(x)})

tax1a<-as.data.frame(tax_table(BATS_1aRel)[,4])

a1_ASV<-as.data.frame(otu_table(BATS_1aRel))
a1_ASV[ "Taxa" ] <- tax1a[,1]
#dim(a1_ASV)
#2 147
ASV_frame1a <- a1_ASV[,-147]
rownames(ASV_frame1a) <- a1_ASV[,147]
ASV_frw1a<-t(ASV_frame1a)

md_to_add<-as.data.frame(sample_data(phyTS5mpruneV2))[,c(7,9,32)] #Date(51) and Temperature Checar esto xq no lo tengo. 
final_2a<-cbind(ASV_frw1a,md_to_add)

final_s11_melt<-melt(final_2a,id.vars=c("Date","Temp"), measure.vars = c("SAR11_clade","others"))

final_s11_melt$variable <- factor(final_s11_melt$variable, levels=c("others","SAR11_clade"))

Colors1a<- c("SAR11_clade"="#2a4858", "others"="#fbfbf6")

final_s11_melt$value <- final_s11_melt$value*100

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

S11_91_94<-ggplot(final_s11_melt, aes(x = Date, y = value)) + scale_x_date(limits=start.end1, date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+ 
geom_area(aes(color = variable, fill = variable),alpha = 0.55) +scale_color_manual(values = Colors1a) +scale_fill_manual(values =Colors1a)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=12,colour = "black"), axis.text.x=element_text(size=8,colour = "black"))+
      labs(title = "1991-1994",
           y = "Relative contribution (%)")+geom_vline(xintercept=as.Date(c("1992-01-01","1993-01-01","1994-01-01")),linetype=4, colour="black")

Forlegend<-ggplot(final_s11_melt, aes(x = Date, y = value)) + scale_x_date(limits=start.end1, date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+ 
geom_area(aes(color = variable, fill = variable),alpha = 0.55) +scale_color_manual(values = Colors1a) +scale_fill_manual(values =Colors1a)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=12,colour = "black"), axis.text.x=element_text(size=8,colour = "black"), legend.position="right",legend.text=element_text(size=11))+
      labs(title = "1991-1994",
           y = "Relative contribution (%)")+geom_vline(xintercept=as.Date(c("1992-01-01","1993-01-01","1994-01-01")),linetype=4, colour="black")

mylegend<-g_legend(Forlegend)

S11_97_2004<-ggplot(final_s11_melt, aes(x = Date, y = value)) + scale_x_date(limits=start.end2, date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+ 
geom_area(aes(color = variable, fill = variable),alpha = 0.55) +scale_color_manual(values = Colors1a) +scale_fill_manual(values =Colors1a)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y= element_blank(), axis.text.x=element_text(size=8,colour = "black"))+
      labs(title = "1997-2004",
           y = NULL)+geom_vline(xintercept=as.Date(c("1997-01-01","1998-01-01","1999-01-01","2000-01-01","2001-01-01","2002-01-01","2003-01-01")),linetype=4, colour="black")

S11_2016_2018<-ggplot(final_s11_melt, aes(x = Date, y = value)) + scale_x_date(limits=start.end3, date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+ 
geom_area(aes(color = variable, fill = variable),alpha = 0.55) +scale_color_manual(values = Colors1a) +scale_fill_manual(values =Colors1a)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y= element_blank(), axis.text.x=element_text(size=8,colour = "black"))+
      labs(title = "2016-2018",
           y = NULL)+geom_vline(xintercept=as.Date(c("2017-01-01","2018-01-01")),linetype=4, colour="black")

totalS11<-grid.arrange(arrangeGrob(S11_91_94+ theme(legend.position="none"),S11_97_2004+ theme(legend.position="none"),S11_2016_2018+ theme(legend.position="none"), widths =c(.84, 1.575,.825)),mylegend, ncol=2,widths=c(10, 1))



#######Panel B SAR11 subclades relative to the total SAR11#######

SAR11B5 = subset_taxa(phyTS5mpruneV2, Order=="SAR11_clade")

BATS_SAR11B5rel<-transform_sample_counts(SAR11B5, function(x){x / sum(x)}) #Transform to % only based on the contribution to total sar11

glomV1filtB5<-tax_glom(BATS_SAR11B5rel, taxrank="tax")

SAR11B5tax<-as.data.frame(tax_table(glomV1filtB5)[,8])

SAR11B5_ASV<-as.data.frame(otu_table(glomV1filtB5))
SAR11B5_ASV[ "Taxa" ] <- SAR11B5tax[,1]
#dim(SAR11B5_ASV)
#
ASV_frame2 <- SAR11B5_ASV[,-147]
rownames(ASV_frame2) <- SAR11B5_ASV[,147]
ASV_frw<-t(ASV_frame2)
md_to_add<-as.data.frame(sample_data(glomV1filtB5))[,c(7,9,32)] #Date,MLD and Temperature Checar esto xq no lo tengo. 
final_2a<-cbind(ASV_frw,md_to_add)

fn2a_melt<-melt(final_2a,id.vars=c("Date","Temp","MLD"),measure.vars = c("I","I;Ia,Ia.1","I;Ia,Ia.3","I;Ia,Ia.4","I;Ib","I;Ib;Ib.1","I;Ib;Ib.2","I;Ic;Ic.1","I;Ic;Ic.2","II","II;IIa;IIa.A","II;IIa;IIa.B","II;IIa;IIa_NA","II;IIb","III;IIIa","clade N2","IV"))

#Colors for SAR11 subclasses
coloresbarplot = c("I;Ia,Ia.1"="darkcyan","I;Ia,Ia.3"="blue","I"="black","I;Ia,Ia.4"="lightskyblue","I;Ib;Ib.1"="darkred","I;Ib;Ib.2"="firebrick1","II;IIa;IIa.A"="bisque4","II;IIa;IIa.B"="darkgoldenrod3 ","II;IIa;IIa_NA"="gold2","III;IIIa"="cornsilk2","II"="chocolate4","II;IIb"="olivedrab","IV"="darkgreen","I;Ic;Ic.2"="lightcoral","I;Ic;Ic.1"="hotpink","clade N2"="lavenderblush2", "I;Ib"="#d45d41")


fn2a_melt$value <- fn2a_melt$value*100

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#Rel abund scale = 100, I need to fit the max MLD at this scale and negative so the MLD follows the depth in the y axis
fn2a_melt$MLD<-fn2a_melt$MLD*.4


S11_91_94subclade<-ggplot(fn2a_melt, aes(x = Date, y = value)) + scale_x_date(limits=start.end1, date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+ 
geom_area(aes(color = variable, fill = variable),alpha = 0.7) +scale_color_manual(values = coloresbarplot) +scale_fill_manual(values =coloresbarplot)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=12,colour = "black"), axis.text.x=element_text(size=8,colour = "black"))+
      labs(title = "1991-1994",
           y = "Relative contribution (%)")+geom_vline(xintercept=as.Date(c("1992-01-01","1993-01-01","1994-01-01")),linetype=4, colour="black")

Forlegendsub<-ggplot(fn2a_melt, aes(x = Date, y = value)) + scale_x_date(limits=start.end1, date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+ 
geom_area(aes(color = variable, fill = variable),alpha = 0.7) +scale_color_manual(values = coloresbarplot) +scale_fill_manual(values =coloresbarplot)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=12,colour = "black"), axis.text.x=element_text(size=8,colour = "black"), legend.position="right",legend.text=element_text(size=11))+
      labs(title = "1991-1994",
           y = "Relative contribution (%)")+geom_vline(xintercept=as.Date(c("1992-01-01","1993-01-01","1994-01-01")),linetype=4, colour="black")

mylegendsub<-g_legend(Forlegendsub)

S11_97_2004subclade<-ggplot(fn2a_melt, aes(x = Date, y = value)) + scale_x_date(limits=start.end2, date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+ 
geom_area(aes(color = variable, fill = variable),alpha = 0.7) +scale_color_manual(values = coloresbarplot) +scale_fill_manual(values =coloresbarplot)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y= element_blank(), axis.text.x=element_text(size=8,colour = "black"))+
      labs(title = "1997-2004",
           y = NULL)+geom_vline(xintercept=as.Date(c("1997-01-01","1998-01-01","1999-01-01","2000-01-01","2001-01-01","2002-01-01","2003-01-01")),linetype=4, colour="black")

S11_2016_2018subclade<-ggplot(fn2a_melt, aes(x = Date, y = value)) + scale_x_date(limits=start.end3, date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+ 
geom_area(aes(color = variable, fill = variable),alpha = 0.7) +scale_color_manual(values = coloresbarplot) +scale_fill_manual(values =coloresbarplot)+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y= element_blank(), axis.text.x=element_text(size=8,colour = "black"))+
      labs(title = "2016-2018",
           y = NULL)+geom_vline(xintercept=as.Date(c("2017-01-01","2018-01-01")),linetype=4, colour="black")

subclade11plot<-grid.arrange(arrangeGrob(S11_91_94subclade+ theme(legend.position="none"),S11_97_2004subclade+ theme(legend.position="none"),S11_2016_2018subclade+ theme(legend.position="none"), widths = c(.84, 1.575,.825)),mylegendsub, ncol=2,widths=c(10, 1))

#######Temperature
tmp1<-ggplot(fn2a_melt, aes(x = Date, y = Temp)) + geom_line(aes(group=1))+ scale_x_date(limits=start.end1,date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_text(size=12,colour = "black"), axis.text.x=element_text(size=8,colour = "black"), legend.position="right")+
      labs(title = "1991-1994",
           y = "Temperature [°C]")+geom_vline(xintercept=as.Date(c("1992-01-01","1993-01-01","1994-01-01")),linetype=4, colour="black")+ylim(15,30)

tmp2<-ggplot(fn2a_melt, aes(x = Date, y = Temp)) + geom_line(aes(group=1))+ scale_x_date(limits=start.end2,date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_blank(), axis.text.x=element_text(size=8,colour = "black"), legend.position="right")+labs(title = "1997-2004",y = NULL)+geom_vline(xintercept=as.Date(c("1997-01-01","1998-01-01","1999-01-01","2000-01-01","2001-01-01","2002-01-01","2003-01-01")),linetype=4, colour="black")+ylim(15,30)

tmp3<-ggplot(fn2a_melt, aes(x = Date, y = Temp)) + geom_line(aes(group=1))+ scale_x_date(limits=start.end3,date_breaks = "1 month", labels = dte_formatter, expand = c(0,0))+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x =element_blank(),axis.title.y=element_text(size=14,colour = "black"),axis.text.y=element_blank(), axis.text.x=element_text(size=8,colour = "black"), legend.position="right")+
      labs(title = "2016-2018",
           y = NULL)+geom_vline(xintercept=as.Date(c("2017-01-01","2018-01-01")),linetype=4, colour="black")+ylim(15,30)

temp<-grid.arrange(arrangeGrob(tmp1+ theme(legend.position="none"),tmp2+ theme(legend.position="none"),tmp3+ theme(legend.position="none"),widths = c(.84, 1.575,.825)),mylegend, ncol=2,widths=c(10, 1))

 
######FINAL#####
plot_grid(totalS11,temp,subclade11plot, align = "v",axis="lr", nrow = 3)


####Numbers extracted for the paper:

#final_s11_melt is divided in SAR11 %s and others %. Extract SAR11 and then use the total
sar11_checktext <- subset(final_s11_melt, variable == "SAR11_clade")
min(sar11_checktext$value) # 4.406329
max(sar11_checktext$value) # 65.35184

#Get max-min of dates between 1991 to 2004 and 

S11_454_text <- subset(sar11_checktext, Date > as.Date("1991-08-01") & Date < as.Date("2004-01-01") )
max(S11_454_text$value) # 62.83762
min(S11_454_text$value) # 4.406329

#broad distribution of values 
sd(S11_454_text$value)
mean(S11_454_text$value)

S11_Illu_text <- subset(sar11_checktext, Date > as.Date("2016-03-08") & Date < as.Date("2018-12-15") )
max(S11_Illu_text$value) # 65.35184
min(S11_Illu_text$value) # 19.51858

length(S11_Illu_text$value[S11_Illu_text$value>40]) # 30 
length(S11_Illu_text$value) # 35

###Ia ecotypes -> SAR11 Ia ecotypes dominate the surface, composing X% of the total BATS (Fig 1b). 
#SUM all Ia ecotypes
#Extract c("I;Ia,Ia.1","I;Ia,Ia.3","I;Ia,Ia.4")
sum(colSums(otu_table(subset_taxa(SAR11B5, Genus=="Ia")))) #SUM=600405
#Sum of all SAR11 reads 1225300. 
sum(colSums(otu_table(SAR11B5))) #SUM=1225300
#600405/1225300 = 0.4900065

sum(colSums(otu_table(subset_taxa(SAR11B5, tax== "I;Ia,Ia.3"))))
#[1] 536531
