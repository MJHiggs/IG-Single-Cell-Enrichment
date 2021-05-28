#### STAGE 0 - NECESSARY SETUP (RUN EVERY TIME) ###############################################################

#read in packages#
library(tidyverse)
library(data.table)
library(liger)
library(Seurat)

#set working directory/ignore if running from current directory#
setwd(...)

#Download Glut and GABA normalised matrix and cluster mappings from#
#https://thejacksonlaboratory.ent.box.com/v/mickelsen-nat-neuro-2019-data/#
#In the current Directory#

#"Imprinted_Gene_List.csv" - read in premade Imprinted gene file- gene name, ensmbl, chromosome order and sex_bias#
IG <- read.csv(...)

#Create an Outputs directory to save pending analysis#
ifelse(!dir.exists(file.path(getwd(), "Outputs")), dir.create(file.path(getwd(), "Outputs")), FALSE)

### STAGE 1 DATA PREPARATION AND UPREGULATION ANALYSIS ###################################################################################

#Read in Expression Matrix for neurons from https://thejacksonlaboratory.ent.box.com/v/mickelsen-nat-neuro-2019-data/#
glut <- read.csv("Glut-normalized-counts.csv")
gaba <- read.csv("GABA-normalized-counts.csv")

#Read in Cluster identities for 2 datasets above from https://thejacksonlaboratory.ent.box.com/v/mickelsen-nat-neuro-2019-data/#
#and append the relevant marker identity to the cluster names e.g. glut and gaba#
cluster_glut <- read.csv("Glut-cluster-mapping.csv")
cluster_glut$dbCluster <- paste("Glut_", cluster_glut$dbCluster, sep = "")
cluster_gaba <- read.csv("GABA-cluster-mapping.csv")
cluster_gaba$dbCluster <- paste("Gaba_", cluster_gaba$dbCluster, sep = "")


#Set up dictionary to name neurons based on cell barcodes using gaba and glut datasets#
dict <- c(cluster_gaba$dbCluster, cluster_glut$dbCluster)
names(dict) <- c(as.character(cluster_gaba$X), as.character(cluster_glut$X))
names(dict)<-str_replace(names(dict), "-", ".")

#Set up Gaba Dataset and add iden column#
colnames <- gaba$X
gaba$X = NULL
Gaba <- data.frame(t(gaba))
colnames(Gaba) <- colnames

#Set up Glut Dataset and add iden column#
colnames <- glut$X
glut$X = NULL
Glut <- data.frame(t(glut))
colnames(Glut) <- colnames

#Acquire Gene names from Ensmebl IDs#
library('biomaRt')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- unique(c(colnames(Gaba), colnames(Glut)))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mart)

#Create dictionary of gene names to map onto the ensembl IDs
genes_dict <- G_list$mgi_symbol
names(genes_dict) <- G_list$ensembl_gene_id

#Establish two datasets to loop over#
data_names <- c("Gaba", "Glut")

for (a in 1:length(data_names)){
  #Read in data#
  wil <- eval(parse(text = data_names[a]))
  #Remove 2 summary columns at the end#
  wil <- wil[,c(-ncol(wil)+1, -ncol(wil))]
  #Create logical vector for genes expressed in 20 or more cells#
  percell <- colSums(wil != 0, na.rm = TRUE) >= 20
  wil <- wil[,percell]
  wil <- wil[!(is.na(genes_dict[colnames(wil)])),]
  wil$iden <- dict[rownames(wil)]
  wil <- wil[wil$iden != "Gaba_0" & wil$iden != "Glut_0",]
  
  #Ensure only complete cases are used#
  wil <- wil[complete.cases(wil), ]
  
  #Create iden - list of unique cell identities to loop over#
  iden <- sort(unique(wil$iden))
  
  #Create a specific folder in Outputs to save this data type's files to#
  ifelse(!dir.exists(file.path(paste(getwd(), "/Outputs", sep = ""), data_names[a])), dir.create(file.path(paste(getwd(), "/Outputs", sep = ""), data_names[a])), FALSE)
  
  #write table of number of cells per cell identity#
  write.csv(table(wil$iden), paste("Outputs/", data_names[a], "/Cell_Subtypes_Frequencies.csv", sep = ""))
  
### ONE-SIDED WILCOXON TEST ###################################################################################
  
  #Precreate dataframes for p, q, fc, pct.in, pct.rest values#
  p <- data.frame("gene" = colnames(wil)[-ncol(wil)], stringsAsFactors = FALSE)
  fc <- data.frame("gene" = colnames(wil)[-ncol(wil)], stringsAsFactors = FALSE)
  percent_in <- data.frame("gene" = colnames(wil)[-ncol(wil)], stringsAsFactors = FALSE)
  percent_rest <- data.frame("gene" = colnames(wil)[-ncol(wil)], stringsAsFactors = FALSE)
  q <- data.frame("gene" = colnames(wil)[-ncol(wil)], stringsAsFactors = FALSE)
  
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
  
  #Change colnames to identity values for p, pct.in and pct.rest and save these files#
  colnames(p) <- c("gene", iden)
  fwrite(p, paste("Outputs/", data_names[a], "/LHA_p.csv", sep=""))
  colnames(percent_in) <- c("gene", iden)
  fwrite(percent_in, paste("Outputs/", data_names[a], "/LHA_percent_in.csv", sep=""))
  colnames(percent_rest) <- c("gene", iden)
  fwrite(percent_rest, paste("Outputs/", data_names[a], "/LHA_percent_rest.csv", sep=""))
  
  #Change colnames to identity values#
  colnames(fc) <- c("gene", iden)
  #set infinite FC values to the max for that identity types#
  for(x in 2:ncol(fc)){
    m = max(fc[is.finite(fc[,x]),x])
    fc[is.infinite(fc[,x]),x] = m
  }
  
  #Save this file#
  fwrite(fc, paste("Outputs/", data_names[a], "/LHA_fc.csv", sep=""))
  
  # adjust p value horizontally using BH correcton and update this to the q dataframe#
  for (x in 1:nrow(p)){
    q[x,2:ncol(p)] <- p.adjust(as.numeric(p[x,2:ncol(p)]), method = "BH")
  }
  
  #Change colnames to identity values and save the file#
  colnames(q) <- c("gene", iden)
  fwrite(q, paste("Outputs/", data_names[a], "/LHA_q.csv", sep=""))
  
### STAGE 2 ENRICHMENT ANALYSIS ########################################################################################

  #define dataframes for files to later save#
  Fish <- data.frame("Identity" = colnames(q[,-1]), stringsAsFactors = FALSE)
  Genelist <- data.frame("Identity" = colnames(q[,-1]), stringsAsFactors = FALSE)
  DEGs <- data.frame()
  IGs <- data.frame()
  
  #loop through identity groups (columns in the p/fc/q etc. files)#
  for (e in 2:(ncol(q))){
    ORA <- data.frame(gene = q$gene, p = p[,e], q = q[,e], fc = fc[,e], pct.in = percent_in[,e], pct.rest = percent_rest[,e], stringsAsFactors = FALSE)
    
    #filter ORA by significant q values and fc limit#
    ORA <- ORA[ORA$q <= 0.05 & ORA$fc >= 1,]
    
    #create repeat column of identity name#
    ORA$Identity = colnames(q)[e]
    
    #Add rows to DEGS - all differentialy expressed genes for that tissue and the total number to Fish#
    DEGs <- rbind(DEGs, ORA)
    Fish[e-1,2] = nrow(ORA)
    
    #Filter ORA for imprinted genes only and add that to IGs file and the total number to Fish#
    ig <- ORA %>% filter(gene %in% IG$ï..Ensembl)
    ig$gene <- genes_dict[as.character(ig$gene)]
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
  fwrite(Genelist, paste("Outputs/", data_names[a],"/IG_Subpopulation_List.csv", sep =""))
  DEGs$gene_name <- genes_dict[as.character(DEGs$gene)]
  fwrite(DEGs, paste("Outputs/", data_names[a], "/Upregulated_Genes.csv", sep =""))
  fwrite(IGs, paste("Outputs/", data_names[a], "/Upregulated_IGs.csv", sep =""))
  
  ### OVER REPRESENTATION ANALYSIS (ORA) ########################################################################
  
  #Create an identity group filter based on having a minimum of 5 imprinted genes upregulated#
  tissue_ORA <- Fish$Identity[Fish$IG >= (as.numeric(sum(q$gene %in% IG$ï..Ensembl))/20)]
  
  #As long as one cell identity has 5 or more IGs, run a Fisher's Exact test#  
  if(length(tissue_ORA) > 0){
    #For every identity, check if it is in the tissue filter#
    for (c in 1:nrow(Fish)){
      if(Fish$Identity[c] %in% tissue_ORA){
        #Create the 2x2 matrix for the fisher's exact test - no. IG in group, no. of IG left over, no. of UpRegulated genes in group minus no. of IGs, all genes left over
        test <- data.frame(gene.interest=c(as.numeric(Fish[c,3]),as.numeric(sum(q$gene %in% IG$ï..Ensembl) - Fish[c,3])), gene.not.interest=c((as.numeric(Fish[c,2]) - as.numeric(Fish[c,3])),as.numeric(nrow(q) - Fish[c,2])- as.numeric(sum(q$gene %in% IG$ï..Ensembl) - Fish[c,3])))
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
  rest <- DEGs[!(DEGs$gene %in% IG$ï..Ensembl),] %>% group_by(Identity) %>% summarise("mean" = mean(fc))
  
  #Merge Fish with rest to have a Mean Fold change Rest of genes column#
  Fish <- merge(Fish, rest, by = "Identity")
  names(Fish)[length(names(Fish))] <- "Mean FC Rest"
  
  ### GENE SET ENRICHMENT ANALYSIS (GSEA) ########################################################################  
  
  #Create list of identities that fulfill the criteria for GSEA - more than 10% of imprinted genes present and a IG FC higher than Rest FC#
  tissue_GSEA <- Fish$Identity[(Fish$`Mean FC IG` >= Fish$`Mean FC Rest` & Fish$IG >= 15)]
  
  ##If there are tissues in tissue_GSEA run a GSEA analysis for IGs in each of those tissues##
  if(length(tissue_GSEA) > 0){
    ifelse(!dir.exists(file.path(paste(getwd(), "/Outputs/", data_names[a], sep = ""), "GSEA")), dir.create(file.path(paste(getwd(), "/Outputs/", data_names[a], sep = ""), "GSEA")), FALSE)
    for(f in 1:length(tissue_GSEA)){
      #filter DEGs for the tissue selected#
      daa <- DEGs %>% filter(DEGs$Identity == tissue_GSEA[f])
      #create log transformed fc values for DEGs#
      file <- log2(daa$fc)
      #replace infinity values with max values and name the file with gene names#
      file[is.infinite(file)] <- max(file[is.finite(file)])
      file <- set_names(file, daa$gene_name)
      #create imprinted gene list to use for GSEA#
      ig <- as.character(IG$Gene)
      #save the graph from the GSEA as a pdf#
      pdf(paste("Outputs/", data_names[a],"/GSEA/",gsub("[^A-Za-z0-9 ]","",tissue_GSEA[f]), "_GSEA.pdf", sep=""))
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
  
  ##Create D - file of average expression of imprinted genes for each tissue##
  
  #Filter Main data for Imprinted Genes and create Identity column#
  D <- wil[,colnames(wil) %in% c(as.character(IG$ï..Ensembl), "iden")] %>% gather("gene", "reads", -iden)
  D$gene <- genes_dict[D$gene]
  #Gather this data and group by gene, identity and get per gene per identity read averages#
  D <- D %>% group_by(gene, iden) %>% summarise("avg" = mean(as.numeric(reads)))
  #Gather fc file by iden and value and filter that for imprinted genes#
  FC <- fc %>% gather(iden, fc, -gene) %>% filter(gene %in% IG$ï..Ensembl)
  FC$gene <- genes_dict[as.character(FC$gene)]
  #Join the two datasets#
  D2 <- left_join(D, FC, by = c("gene","iden"))
  
  #create u with the tissue with max expression for each IG and save that file#
  u <- D2 %>% group_by(gene) %>% filter(gene %in% IG$Gene & (avg == max(avg) & max(avg) > 0))
  fwrite(u, paste("Outputs/", data_names[a], "/IG_Top_Expressed_Subpopulation.csv", sep =""))

  #write the finished Fish file#
  fwrite(Fish, paste("Outputs/", data_names[a], "/Enrichment_Analysis.csv", sep =""))  
}