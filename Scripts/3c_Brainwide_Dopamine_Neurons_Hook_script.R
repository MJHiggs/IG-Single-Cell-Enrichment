#### STAGE 0 - NECESSARY SETUP (RUN EVERY TIME) ###############################################################

#read in packages#
library(tidyverse)
library(liger)
library(data.table)

#set working directory/ignore if running from current directory#
setwd(...)

#"Imprinted_Gene_List.csv" - read in premade Imprinted gene file- gene name, ensmbl, chromosome order and sex_bias#
IG <- read.csv(...)

#Create an Outputs directory to save pending analysis#
ifelse(!dir.exists(file.path(getwd(), "Outputs")), dir.create(file.path(getwd(), "Outputs")), FALSE)

### STAGE 1 DATA PREPARATION AND UPREGULATION ANALYSIS ###################################################################################

#read in datafiles split by developmental age from github#
E15_Raw <- readRDS(gzcon(url("https://github.com/pwh124/sc-da-parkinsons/blob/master/data/e15.dat.filter.final.rds?raw=true")))
P7_Raw <- readRDS(gzcon(url("https://github.com/pwh124/sc-da-parkinsons/blob/master/data/P7.Fb.dat.filter.final.Rds?raw=true")))

#Select out expression data#
E15 <- E15_Raw@assayData$exprs[]
P7 <- P7_Raw@assayData$exprs[]

#Create cluster identity dictionaries for both datasets#
ClusterE15 <- as.character(E15_Raw@phenoData@data$subset.cluster)
names(ClusterE15) <- as.character(E15_Raw@phenoData@data$cell_id)
ClusterP7 <- as.character(P7_Raw@phenoData@data$subset.cluster)
names(ClusterP7) <- as.character(P7_Raw@phenoData@data$cell_id)

#Merge the dictionaries#

#Create gene name dictionaries for both datasets#
GenenamesE15 <- as.character(E15_Raw@featureData@data$gene_short_name)
names(GenenamesE15) <- as.character(E15_Raw@featureData@data$gene_id)
GenenamesP7 <- as.character(P7_Raw@featureData@data$gene_short_name)
names(GenenamesP7) <- as.character(P7_Raw@featureData@data$gene_id)

#Merge the genes to those that are present in both datasets#
E15 <- E15[names(GenenamesE15) %in% names(GenenamesP7),]

#Bind the two datasets together into one#
Final_data <- cbind(E15, P7)

#Create the 20 cell expression gene filter and filter data and gene names by the filter#
Filter <-(Matrix::rowSums(Final_data!=0) >= 20)
Final_data <- Final_data[Filter,]
GenenamesP7 <- GenenamesP7[Filter]

#Transpose the final data and turn into dataframe#
Final_data <- Matrix::t(Final_data) #396, 13095 
wil <- as.data.frame(as.matrix(Final_data))

#Create iden column with cell identities in#
wil$iden <- cluster_dict[rownames(wil)]

#Ensure only complete cases are used#
wil <- wil[complete.cases(wil), ]

#write table of number of cells per cell identity#
write.csv(table(wil$iden), paste("Outputs/Cell_Subtype_Frequencies.csv", sep = ""))

### ONE-SIDED WILCOXON TEST ###################################################################################

#Precreate dataframes for p, q, fc, pct.in, pct.rest values#
p <- data.frame("ensmbl" = colnames(wil)[-ncol(wil)], stringsAsFactors = FALSE)
fc <- data.frame("ensmbl" = colnames(wil)[-ncol(wil)], stringsAsFactors = FALSE)
percent_in <- data.frame("ensmbl" = colnames(wil)[-ncol(wil)], stringsAsFactors = FALSE)
percent_rest <- data.frame("ensmbl" = colnames(wil)[-ncol(wil)], stringsAsFactors = FALSE)
q <- data.frame("ensmbl" = p$ensmbl, stringsAsFactors = FALSE)

#Loop over every gene/column apart from the final 'iden' column#
for (b in 1:(ncol(wil)-1)){
  
  #subset the reads for this gene and the cell identities#
  dd <- wil[,c(b, ncol(wil))]
  colnames(dd) <- c("reads", "iden")
  
  #Loop through each cell identity#
  for (c in 1:length(iden)){
    #extract the reads from identity of interest#
    interest <- as.numeric(as.character(dd[dd$iden == iden[c],]$reads))
    #extract the reads from the other identities#
    rest <- as.numeric(as.character(dd[dd$iden != iden[c],]$reads))
    #Run the Wilcoxon test, one-sided#
    W <- wilcox.test(interest, rest, exact = FALSE, alternative = "greater")
    #Update the cell in the precreated database with the p value#
    p[b,c+1] <- as.numeric(W$p.value)
    #Divide Mean interest by Mean rest to create FC value#
    fc[b,c+1] <- mean(interest)/mean(rest)
    #divide number of non-zero cells by total number of cells to create percent values#
    percent_in[b,c+1] <- sum(interest != 0)/ length(interest)
    percent_rest[b,c+1] <- sum(rest != 0)/ length(rest)
  }
  #Print b as a progress counter#
  print(b)
}

#Change colnames to identity values for p, pct.in and pct.rest, change ensembl identities to gene names and save these files#
colnames(p) <- c("ensmbl", iden)
fwrite(p, paste("Outputs/Hook_p.csv", sep=""))
colnames(percent_in) <- c("ensmbl", iden)
fwrite(percent_in, paste("Outputs/Hook_percent_in.csv", sep=""))
colnames(percent_rest) <- c("ensmbl", iden)
fwrite(percent_rest, paste("Outputs/Hook_percent_rest.csv", sep=""))

#Change colnames to identity values
colnames(fc) <- c("ensmbl", iden)
#set infinite FC values to the max for that identity types#
for(x in 2:ncol(fc)){
  m = max(fc[is.finite(fc[,x]),x])
  fc[is.infinite(fc[,x]),x] = m
}
#Save this file#
fwrite(fc, paste("Outputs/Hook_fc.csv", sep=""))

# adjust p value horizontally using BH correcton and update this to the q dataframe#
for (x in 1:nrow(p)){
  q[x,2:ncol(p)] <- p.adjust(as.numeric(p[x,2:ncol(p)]), method = "BH")
}
#Change colnames to identity values and save the file#
colnames(q) <- c("ensmbl", iden)
q$ensmbl <- sapply(strsplit(as.character(q$ensmbl), "[.]"), `[`, 1)
fwrite(q, paste("Outputs/Hook_q.csv", sep=""))

### STAGE 2 ENRICHMENT ANALYSIS ########################################################################################

#define dataframes for files to later save#
Fish <- data.frame("Identity" = colnames(q[,c(-1, -ncol(q))]), stringsAsFactors = FALSE)
Genelist <- data.frame("Identity" = colnames(q[,c(-1, -ncol(q))]), stringsAsFactors = FALSE)
Marker_Genelist <- data.frame("Identity" = colnames(q[,c(-1, -ncol(q))]), stringsAsFactors = FALSE)
DEGs <- data.frame()
IGs <- data.frame()

#loop through identity groups (columns in the p/fc/q etc. files)#
for (e in 2:(ncol(q)-1)){
  ORA <- data.frame(gene = p$ensmbl, p = p[,e], q = q[,e], fc = fc[,e], pct.in = percent_in[,e], pct.rest = percent_rest[,e], stringsAsFactors = FALSE)
  
  #filter ORA by significant q values and fc limit#
  ORA <- ORA[ORA$q <= 0.05 & ORA$fc >= 1,]
  
  #create repeat column of identity name#
  ORA$Identity = colnames(q)[e]
  
  #Add rows to DEGS - all differentialy expressed genes for that tissue and the total number to Fish#
  DEGs <- rbind(DEGs, ORA)
  Fish[e-1,2] = nrow(ORA)
  
  #Filter ORA for imprinted genes only and add that to IGs file and the total number to Fish#
  ig <- ORA %>% filter(sapply(strsplit(as.character(ORA$gene), "[.]"), `[`, 1) %in% IG$ï..Ensmbl)
  ig$gene <- GenenamesP7[as.character(ig$gene)]
  IGs <- rbind(IGs, ig)
  Fish[e-1,3] = nrow(ig)
  
  #Add to the IG genelist, only if there are significant genes for that tissue#
  if (nrow(ig) > 0){
    for (z in 1:nrow(ig)){
      Genelist[e-1, z+1] <- as.character(ig$gene[z])
    }
  }
}

#add colnames to Fish and save all the files generated other than Fish#
colnames(Fish)<- c("Identity","UpReg", "IG")
fwrite(Genelist, paste("Outputs/IG_Subpopulation_List.csv", sep =""))
DEGs$gene_name <- GenenamesP7[DEGs$gene]
fwrite(DEGs, paste("Outputs/Upregulated_Genes.csv", sep =""))
fwrite(IGs, paste("Outputs/Upregulated_IGs.csv", sep =""))

### OVER REPRESENTATION ANALYSIS (ORA) ########################################################################

#Create an identity group filter based on having a minimum of 5 imprinted genes upregulated#
tissue_ORA <- Fish$Identity[Fish$IG >= (as.numeric(sum(q$gene %in% IG$Gene))/20)]

#As long as one cell identity has 5 or more IGs, run a Fisher's Exact test#  
if(length(tissue_ORA) > 0){
  #For every identity, check if it is in the tissue filter#
  for (c in 1:nrow(Fish)){
    if(Fish$Identity[c] %in% tissue_ORA){
      #Create the 2x2 matrix for the fisher's exact test - no. IG in group, no. of IG left over, no. of UpRegulated genes in group minus no. of IGs, all genes left over
      test <- data.frame(gene.interest=c(as.numeric(Fish[c,3]),as.numeric(sum(q$ensmbl %in% IG$ï..Ensmbl) - Fish[c,3])), gene.not.interest=c((as.numeric(Fish[c,2]) - as.numeric(Fish[c,3])),as.numeric(nrow(q) - Fish[c,2])- as.numeric(sum(q$ensmbl %in% IG$ï..Ensmbl) - Fish[c,3])))
      row.names(test) <- c("In_category", "not_in_category")
      f <- fisher.test(test, alternative = "greater")
      #Update Fish with 
      Fish[c,4] <- f$p.value
    }
  }
  names(Fish)[length(names(Fish))] <- "ORA_p"
  
  ## Add a new column to Fish with adjusted Fisher p values (bonferroni correction)##
  Fish$Adjust_ORA <- p.adjust(Fish$ORA_p, method = "bonferroni")
  names(Fish)[length(names(Fish))] <- "ORA_q"
}

### MEAN FOLD CHANGE #########################################################################################  

##Calculate mean fold change for upregulated imprinted genes in each identity group#
igs <- (IGs %>% group_by(Identity) %>% summarise(mean = mean(fc)))

#Input tissues without Imprinted Genes as 0 in igs#
if (length(unique(DEGs$Identity)) != length(unique(IGs$Identity))){
  leftout <- data.frame("Identity" = unique(DEGs$Identity)[!(unique(DEGs$Identity) %in% unique(IGs$Identity))])
  leftout$mean = 0
  igs <- rbind(igs, leftout)
}

#Merge Fish with igs to have a Mean Fold change imprinted gene column#
Fish <- merge(Fish, igs, by = "Identity")
names(Fish)[length(names(Fish))] <- "Mean FC IG"

##Calculate mean fold change for the rest of the upregulated genes per tissue and update a column in Fish## 
rest <- DEGs[!(sapply(strsplit(as.character(DEGs$gene), "[.]"), `[`, 1)  %in% IG$ï..Ensmbl),] %>% group_by(Identity) %>% summarise("mean" = mean(as.numeric(fc)))

#Merge Fish with rest to have a Mean Fold change Rest of genes column#
Fish <- merge(Fish, rest, by = "Identity")
names(Fish)[length(names(Fish))] <- "Mean FC Rest"

### GENE SET ENRICHMENT ANALYSIS (GSEA) ########################################################################  

#Create list of identities that fulfill the criteria for GSEA - more than 10% of imprinted genes present and a IG FC higher than Rest FC#
tissue_GSEA <- Fish$Identity[(Fish$`Mean FC IG` >= Fish$`Mean FC Rest` & Fish$IG >= 15)]

##If there are tissues in tissue_GSEA run a GSEA analysis for IGs in each of those tissues##
if(length(tissue_GSEA) > 0){
  ifelse(!dir.exists(file.path(paste(getwd(), "/Outputs/",sep = ""), "GSEA")), dir.create(file.path(paste(getwd(), "/Outputs/", sep = ""), "GSEA")), FALSE)
  for(f in 1:length(tissue_GSEA)){
    #filter DEGs for the tissue selected#
    daa <- DEGs %>% filter(DEGs$Identity == tissue_GSEA[f])
    #create log transformed fc values for DEGs#
    file <- log2(as.numeric(daa$fc))
    #replace infinity values with max values and name the file with gene names#
    file[is.infinite(file)] <- max(file[is.finite(file)])
    file <- set_names(file, daa$gene_name)
    #create imprinted gene list to use for GSEA#
    ig <- as.character(IG$Gene)
    #save the graph from the GSEA as a pdf#
    pdf(paste("Outputs/GSEA/",gsub("[^A-Za-z0-9 ]","",tissue_GSEA[f]), "_GSEA.pdf", sep=""))
    #Add GSEA p value to Fish#
    Fish[Fish$Identity == tissue_GSEA[f],8] <- gsea(file, ig)
    #Save PDF#
    dev.off()
  }
  
  names(Fish)[length(names(Fish))] <- "GSEA_p"
  
  #Adjust the GSEA p values using bonferroni as new column in Fish#
  Fish$Adjust_GSEA <- p.adjust(Fish$GSEA_p, method = "bonferroni")
  names(Fish)[length(names(Fish))] <- "GSEA_q"
}

### IMPRINTED GENE TOP CELL EXPRESSION ######################################################################    

#Filter Main data for Imprinted Genes and create Identity column#
D <- wil[,sapply(strsplit(as.character(colnames(wil)), "[.]"), `[`, 1) %in% IG$ï..Ensmbl]
colnames(D) <- GenenamesP7[colnames(D)]
D$iden <- cluster_dict[rownames(D)]
#Gather this data and group by gene, identity and get per gene per identity read averages#
D <- D %>% gather("gene", "reads", -iden)
D <- D %>% group_by(gene, iden) %>% summarise("avg" = mean(as.numeric(reads)))
#Gather fc file by iden and value and filter that for imprinted genes#
FC <- fc[,-ncol(fc)] %>% gather(iden, fc, -ensmbl) 
FC <- FC %>% filter(sapply(strsplit(as.character(FC$ensmbl), "[.]"), `[`, 1) %in% IG$ï..Ensmbl)
FC$ensmbl <- GenenamesP7[as.character(FC$ensmbl)]
colnames(FC) <- c("gene", "iden", "fc")
#Join the two datasets#
D2 <- left_join(D, FC, by = c("gene","iden"))
D2$fc <- as.numeric(D2$fc)

## create u with the tissue with max expression for each IG and save that file ##
#create u with the tissue with max expression for each IG and save that file#
u <- D2 %>% group_by(gene) %>% filter(gene %in% IG$Gene & (avg == max(avg) & max(avg) > 0))
fwrite(u, paste("Outputs/IG_Top_Expressed_Subpopulation.csv", sep =""))

#write the finished Fish file#
fwrite(Fish, paste("Outputs/Enrichment_Analysis.csv", sep =""))

#### STAGE 3 - VISUALISATION DOTPLOT ########################################################################################################

#Arrange Fish by over-representation significance, take the identity variable as an order value#
Fish <- Fish %>% arrange(Fish$ORA_p)
order <- Fish$Identity

#Filter original IG list with those in the dataset#
IG2 <- IG[IG$Gene %in% D2$gene,]
#Arrange IGs by the chromosomal order#
IG2 <- IG2 %>% arrange(IG2$Order)

#Filter IGs for Paternally expressed genes (PEGs) only and create Pat using the PEG only filter#
pat <- IG2 %>% filter(Sex == "P")
Pat <- D2 %>% filter(gene %in% pat$Gene)

#Recast Gene as a Factor and arrange by chromosomal order#
Pat$gene <- factor(Pat$gene, levels = rev(pat$Gene))
Pat <- Pat %>% arrange(rev(gene))

#Filter IGs for maternally expressed genes (MEGs) only and create Mat using the MEG only filter#
mat <- IG2 %>% filter(Sex == "M" | Sex == "I")
Mat <- D2 %>% filter(gene %in% mat$Gene)

#Recast Gene as a Factor and arrange by chromosomal order#
Mat$gene <- factor(Mat$gene, levels = rev(mat$Gene))
Mat <- Mat %>% arrange(rev(gene))

Pat$fc <- log2(Pat$fc+1)
Mat$fc <- log2(Mat$fc+1)

#Create PDF to save PEG dotplot#
pdf(paste("Outputs/", data_names[a],"/", "PEG_DOTPLOT.pdf", sep=""))

#GGplot dotplot, x = cell identity, y = gene identity, color = fc(gradated up to 5FC+), size = avg expression(0 to max expression registered)#
print(ggplot(Pat, aes(x=iden, y=gene, color=ifelse(fc == 0, NA, fc), size=ifelse(avg==0, NA, avg))) + geom_point(alpha = 0.8) +
        theme_classic() +
        scale_color_gradientn(colours = c("grey95","grey60","blue","darkblue","midnightblue"), na.value="midnightblue", values = c(0, 0.2, 0.4, 0.6, 0.8 ,1), limits = c(0,5)) +
        theme(axis.text.x = element_text(angle = 90)) +
        scale_size_continuous(limits = c(0,max(D2$avg)))+
        scale_x_discrete(limits = order) +
        labs(x = "Cell Identity", y = "Imprinted Gene", size = "Normalised Mean Expression", color = "Log2FC vs Background"))
#Save PDF#
dev.off()

#Create PDF to save MEG dotplot#  
pdf(paste("Outputs/", data_names[a],"/", "MEG_DOTPLOT.pdf", sep=""))

#GGplot dotplot, x = cell identity, y = gene identity, color = fc(gradated up to 5FC+), size = avg expression(0 to max expression registered)# 
print(ggplot(Mat, aes(x=iden, y=gene, color=ifelse(fc == 0, NA, fc), size=ifelse(avg==0, NA, avg))) + geom_point(alpha = 0.8) +
        theme_classic() +
        scale_color_gradientn(colours = c("grey95","grey60","orange","red","darkred"), na.value="darkred", values = c(0, 0.2, 0.4, 0.6, 0.8 ,1), limits = c(0,5)) +
        theme(axis.text.x = element_text(angle = 90)) +
        scale_size_continuous(limits = c(0,max(D2$avg)))+
        scale_x_discrete(limits = order) +
        labs(x = "Cell Identity", y = "Imprinted Gene", size = "Normalised Mean Expression", color = "Log2FC vs Background"))
#Save PDF#
dev.off()
