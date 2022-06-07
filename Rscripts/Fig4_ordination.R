##Figure 4: SAR11 ordination by season and hot-cold clusters (PANELS A,B =BATS, PANELS C,D =WEC) 

#####CREATE THE BATS
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

#####Generate the phyloseq object BATScf #####

count_tabS11phyTS5m <- read.table("BIOSotu5.otu", header=T, row.names=1, check.names=F)
SAR11df_envUPDATE<-read.table("BATS_S11_5m_envfile.txt", header=T, row.names=1, check.names=F, sep ="\t")
tax_tabS11phyTS5m <- as.matrix(read.table("BIOStax5.tax", header=T, row.names=1, check.names=F, sep="\t",na.strings = "#NA"))

#Add 2clusters
SAR11df_envUPDATE$Date<-as.Date(SAR11df_envUPDATE$Date,"%Y-%m-%d")
astroclusters <- as.integer(c("0000","0501", "0922", "1232")) #May the first, here I split the spring into cold (before) and warm (after)
astroclusters_labels <- c("Cold", "Warm", "Cold")
SAR11df_envUPDATE$clusters<-astroclusters_labels[ cut(as.integer(format(SAR11df_envUPDATE$Date, "%m%d")), astroclusters, labels = FALSE) ]

#####Check the names!!!!
rownames(SAR11df_envUPDATE) <- gsub("_", ".", rownames(SAR11df_envUPDATE))


S11SAM_UPDATE = sample_data(SAR11df_envUPDATE)
S11SAM_UPDATE_OTU<-otu_table(count_tabS11phyTS5m, taxa_are_rows = TRUE)
S11SAM_UPDATE_TAX<-tax_table(tax_tabS11phyTS5m)

###UPDATED PHYLOSEQ WITH CURATED ENVIRONMENTAL DATA **ONLY FOR 5m**
phyTS5mUPD<-phyloseq(S11SAM_UPDATE,S11SAM_UPDATE_OTU,S11SAM_UPDATE_TAX)

phyTS5mV1 <- prune_samples(sample_sums(phyTS5mUPD)>=1000, phyTS5mUPD) #prune samples less than 1000 reads (eliminate 17 --> 146/163)

phyTS5mpruneV2<-filter_taxa(phyTS5mV1, function(x) sum(x >= 3) > (0.02 *length(x)), TRUE) #V1=916/2875 Again from phyTS5mprune, V2 (146 samples, .01) =500 ASVs
colors_seasons=c("spring"="#4f9710", "summer"="#E08B1A", "autumn"="#937615", "winter"="#255C97")
colors_clusters=c("Cold"="#544a7d", "Warm"="#FABE6E")

#Curated PHYLOSEQ OBJECT 
SAR11_BATSphy = subset_taxa(phyTS5mpruneV2, Order=="SAR11_clade")

SAR11_BATSphyMiSeq<-subset_samples(SAR11_BATSphy, tech=="MiSeq")

###Rarefy 7200 ######
set.seed(717)
phyeuphminrarWSAR11_BATSphyMiSeq = rarefy_even_depth(SAR11_BATSphyMiSeq, sample.size = 7200)

###Rarefy 465 ######
set.seed(717)
phyeuphminrarWSAR11_BATSphy = rarefy_even_depth(SAR11_BATSphy, sample.size = 465) #Dependiendo de que muestras salgan aca, podemos usar mas variables quedan 211 muestras 

bray_not_na_BATSSAR11 <- phyloseq::distance(physeq = phyeuphminrarWSAR11_BATSphy, method = "bray")


####################CONSTRAINED ordination MiSeq##################
###SEASON PANEL A###
phynaBATS_MiSeq<-phyeuphminrarWSAR11_BATSphyMiSeq %>%
subset_samples(
!is.na(Date)& 
!is.na(NO3)& 
!is.na(NO2)&
!is.na(NO2.NO3)&
!is.na(Sil)&
!is.na(PO4)& 
!is.na(CTD_S)& 
!is.na(O2)& 
!is.na(Temp)&
!is.na(Turner.Chl.a)
)

bray_not_na_phynaBATSMiSeq <- phyloseq::distance(physeq = phynaBATS_MiSeq, method = "bray")

cap_ord_miseq <- ordinate(
physeq = phynaBATS_MiSeq, 
 method = "CAP",
 distance = bray_not_na_phynaBATSMiSeq,
 formula = ~ NO3+NO2+NO2.NO3+Sil+PO4+CTD_S+O2+Temp+Turner.Chl.a)


#PLOT
cap_plot_miseq <- plot_ordination(
physeq = phynaBATS_MiSeq, 
ordination = cap_ord_miseq, 
color = "season", 
 axes = c(1,2)) +  
 geom_point(colour="black", shape=21, size = 3, aes(fill = season))+ 
 scale_fill_manual(values=c(colors_seasons)) + theme(axis.text.x=element_text(angle=0,size=14), legend.title=element_blank(),  axis.text.y=element_text(vjust=0.4, hjust=1,size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+xlim(-2,2)+ylim(-2,2.25)

cap_plot2_miseq<-cap_plot_miseq + stat_ellipse(colour="black",aes(fill = season), geom="polygon",level=0.95,alpha=0.15)

#Arrows
arrowmatmiseq <- vegan::scores(cap_ord_miseq, display = "bp")

# Add labels, make a data.frame
arrowdfmiseq <- data.frame(labels = rownames(arrowmatmiseq), arrowmatmiseq)

# Define the arrow aesthetic mapping
arrow_mapmiseq <- aes(xend = CAP1, 
    yend = CAP2, 
    x = 0, 
    y = 0, 
    shape = NULL, 
    color = NULL, 
    label = labels)

label_mapmiseq <- aes(x = 1.3 * CAP1, 
    y = 1.3 * CAP2, 
    shape = NULL, 
    color = NULL, 
    label = labels)

arrowheadmiseq = arrow(length = unit(0.02, "npc"))

##PANEL A miseq
constrmiseq<-cap_plot2_miseq + 
  geom_segment(
    mapping = arrow_mapmiseq, 
    size = .6, 
    data = arrowdfmiseq, 
    color = "black", 
    arrow = arrowheadmiseq
  ) + 
  geom_text(
    mapping = label_mapmiseq, 
    size = 4,  
    data = arrowdfmiseq, 
    show.legend = FALSE
  )

#PANEL B color-coded by clusters
cap_plot_miseq_clusters <- plot_ordination(
physeq = phynaBATS_MiSeq, 
ordination = cap_ord_miseq, 
color = "clusters", 
 axes = c(1,2)) +  
 geom_point(colour="black", shape=21, size = 3, aes(fill = clusters))+ 
 scale_fill_manual(values=c(colors_clusters)) + theme(axis.text.x=element_text(angle=0,size=14), legend.title=element_blank(),  axis.text.y=element_text(vjust=0.4, hjust=1,size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+xlim(-2,2)+ylim(-2,2.25)

cap_plot2_miseq_clusters<-cap_plot_miseq_clusters + stat_ellipse(colour="black",aes(fill = clusters), geom="polygon",level=0.95,alpha=0.15)

constrmiseq_clusters<-cap_plot2_miseq_clusters + 
  geom_segment(
    mapping = arrow_mapmiseq, 
    size = .6, 
    data = arrowdfmiseq, 
    color = "black", 
    arrow = arrowheadmiseq
  ) + 
  geom_text(
    mapping = label_mapmiseq, 
    size = 4,  
    data = arrowdfmiseq, 
    show.legend = FALSE
  )

#plot_grid(constrmiseq,constrmiseq_clusters) #panels A and B

                     
##############WEC###############

#####Generate the phyloseq object WECcf #####
countwec <- read.table("WEC_mar.otu", header=T, row.names=1, check.names=F)
sampleinfowec <- read.table("WEC_mar.env", header=T, row.names=1, check.names=F, sep ="\t")
taxwec <- as.matrix(read.table("WEC_marV2.tax", header=T, row.names=1, check.names=F, na.strings="", sep="\t")) 

##Add a "warm"/ "cold" cluster label based on the dates
sampleinfowec$Date<-as.Date(sampleinfowec$Date,"%m/%d/%Y")

astroseasons <- as.integer(c("0000","0501", "0922", "1232")) #May the first, here I split the spring into cold (before) and warm (after)
astroseasons_labels <- c("Cold", "Warm", "Cold")
sampleinfowec$clusters<-astroseasons_labels[ cut(as.integer(format(sampleinfowec$Date, "%m%d")), astroseasons, labels = FALSE) ]

#I have worked since the beginning using the four canonical seasons. However, the two transition state of the pattern we are seeing in SAR11 is like autumn+winter+early spring vs late spring+summer (split spring in two). 

OTUmar = otu_table(countwec, taxa_are_rows = TRUE)
TAXmar = tax_table(taxwec)
SAMar= sample_data(sampleinfowec)

WEC<-phyloseq(OTUmar,TAXmar,SAMar)
WECcf= filter_taxa(WEC, function(x) sum(x > 1) > (0.0015*length(x)), TRUE)


####################Constrained ordination SAR11###################
family_counts_tab <- otu_table(tax_glom(WECcf, taxrank="Family")) 

#making a vector of phyla names to set as row names
family_tax_vec <- as.vector(tax_table(tax_glom(WECcf, taxrank="Family"))[,4]) 
rownames(family_counts_tab) <- as.vector(family_tax_vec)


temp_major_taxa_counts_tab <- family_counts_tab[!row.names(family_counts_tab) %in% "SAR11_clade", ] 

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
#SAMar change date column
sampleinfowec$Date<-as.Date(sampleinfowec$Date,"%m/%d/%Y") 
sampleinfowec$year = as.character(format(sampleinfowec$Date, format = "%Y"))
sampleinfowec$month = as.character(format(sampleinfowec$Date, format = "%m"))
sampleinfowec$day = as.character(format(sampleinfowec$Date, format = "%d"))
SAMar= sample_data(sampleinfowec)

SAR11_WECphy2<-phyloseq(SAR11_rest_otutable,SAR11_restTAX_table,SAMar)
SAR11_WECphy = subset_taxa(WECcf, Order=="SAR11_clade")

###Rarefy 2500 ######
set.seed(717)
phyeuphminrarWSAR11_WECphy = rarefy_even_depth(SAR11_WECphy, sample.size = 2500) 

bray_not_na_WECSAR11 <- phyloseq::distance(physeq = phyeuphminrarWSAR11_WECphy, method = "bray")

######CAP SAR11 SEASONS AND CLUSTERS 
#Environmental data available for most of the samples 

phynaWECSAR11 <- phyeuphminrarWSAR11_WECphy%>%
subset_samples(
!is.na(Date)& 
!is.na(NITRITE_0)& 
!is.na(NITRATE.NIT_0)& 
!is.na(AMMONIA_0)& 
!is.na(PHOSPHATE_0)& 
!is.na(salinity_PSU_2.5)& 
!is.na(oxygen_uM_2.5)& 
!is.na(temp_C_2.5)& 
!is.na(Chl_0m)
)



bray_not_naconstr_WECSAR11 <- phyloseq::distance(physeq =phynaWECSAR11, method = "bray")

cap_ordsar11 <- ordinate(
physeq = phynaWECSAR11, 
 method = "CAP",
 distance = bray_not_naconstr_WECSAR11,
 formula = ~ oxygen_uM_2.5+AMMONIA_0+NITRITE_0+NITRATE.NIT_0+PHOSPHATE_0+salinity_PSU_2.5+temp_C_2.5+Chl_0m)

 
#PLOT SEASONS
cap_plotsar11 <- plot_ordination(
physeq = phynaWECSAR11, 
ordination = cap_ordsar11, 
color = "season", 
 axes = c(1,2)) + 
 geom_point(colour="black", shape=21, size = 3, 
     aes(fill = season)) + 
 scale_fill_manual(values=c(colors_seasons))  + theme(axis.text.x=element_text(angle=0,size=14), legend.title=element_blank(),  axis.text.y=element_text(vjust=0.4, hjust=1,size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

cap_plotsar11_2<-cap_plotsar11 + stat_ellipse(colour="black",aes(fill = season), geom="polygon",level=0.95,alpha=0.15)

#PLOT CLUSTERS

cap_plotsar11clusters <- plot_ordination(
physeq = phynaWECSAR11, 
ordination = cap_ordsar11, 
color = "clusters", 
 axes = c(1,2)) + 
 geom_point(colour="black", shape=21, size = 3, 
     aes(fill = clusters)) + 
 scale_fill_manual(values=c(colors_clusters))  + theme(axis.text.x=element_text(angle=0,size=14), legend.title=element_blank(),  axis.text.y=element_text(vjust=0.4, hjust=1,size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))

cap_plotsar11clusters_2<-cap_plotsar11clusters + stat_ellipse(colour="black",aes(fill = clusters), geom="polygon",level=0.95,alpha=0.15)


#Arrows
arrowmatsar11 <- vegan::scores(cap_ordsar11, display = "bp")

# Add labels, make a data.frame
arrowdfsar11 <- data.frame(labels = rownames(arrowmatsar11), arrowmatsar11)

# Define the arrow aesthetic mapping
arrow_mapsar11 <- aes(xend = CAP1, 
    yend = CAP2, 
    x = 0, 
    y = 0, 
    shape = NULL, 
    color = NULL, 
    label = labels)

label_map <- aes(x = 1.3 * CAP1, 
    y = 1.3 * CAP2, 
    shape = NULL, 
    color = NULL, 
    label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

#################################
# Make a new graphic PLOT SEASONS
#################################
SAR11_2500rar<-cap_plotsar11_2 + 
  geom_segment(
    mapping = arrow_mapsar11, 
    size = .6, 
    data = arrowdfsar11, 
    color = "black", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdfsar11, 
    show.legend = FALSE
  )

#################################
# Make a new graphic PLOT CLUSTERS
#################################
SAR11_2500rar_clusters<-cap_plotsar11clusters_2 + 
  geom_segment(
    mapping = arrow_mapsar11, 
    size = .6, 
    data = arrowdfsar11, 
    color = "black", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdfsar11, 
    show.legend = FALSE
  )

###FINAL PLOT 4 panels####
plot_grid(constrmiseq,constrmiseq_clusters,SAR11_2500rar,SAR11_2500rar_clusters,ncol=2,labels = "auto",label_size = 18)

