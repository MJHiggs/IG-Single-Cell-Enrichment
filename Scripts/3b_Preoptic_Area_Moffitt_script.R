#### STAGE 0 - NECESSARY SETUP (RUN EVERY TIME) ###############################################################

#read in packages#
library(tidyverse)
library(data.table)
library(liger)
library(Seurat)
library(Matrix)
BiocManager::install("GEOquery")
library(GEOquery)

#set working directory/ignore if running from current directory#
setwd(...)

#"Imprinted_Gene_List.csv" - read in premade Imprinted gene file- gene name, ensmbl, chromosome order and sex_bias#
IG <- read.csv(...)

#Create an Outputs directory to save pending analysis#
ifelse(!dir.exists(file.path(getwd(), "Outputs")), dir.create(file.path(getwd(), "Outputs")), FALSE)

#Read in the supplementary file from GSE87544#
getGEOSuppFiles("GSE113576")

#Rename the files to remove the prefix#
setwd("GSE113576/")
file.rename(list.files(), gsub("GSE113576_", "", list.files()))
setwd("../")

#Download cell metadata from supplementary table 1 of Moffitt -#
#https://science.sciencemag.org/highwire/filestream/717711/field_highwire_adjunct_files/1/aau5324_Moffitt_Table-S1.xlsx#

#Read in the cell identity table#
ID <- readxl::read_xlsx("aau5324_Moffitt_Table-S1.xlsx")
colnames(ID) <- ID[1,]
ID <- ID[-1,]

### STAGE 1 DATA PREPARATION AND UPREGULATION ANALYSIS ###################################################################################

#Create dictionaries for neuronal identities#
ID$`Cell name` <- gsub("[-]",".",ID$`Cell name`)
ID$`Non-neuronal cluster (determined from clustering of all cells)` <- ifelse(ID$`Non-neuronal cluster (determined from clustering of all cells)` == "", ID$`Neuronal cluster (determined from clustering of inhibitory or excitatory neurons)`, ID$`Non-neuronal cluster (determined from clustering of all cells)`)
Neurons <- ID %>% filter(`Cell class (determined from clustering of all cells)` == "Excitatory" | `Cell class (determined from clustering of all cells)`== "Inhibitory")
Neurons$`Cell name` <- gsub("[-]",".",Neurons$`Cell name`)
dict <- Neurons$`Neuronal cluster (determined from clustering of inhibitory or excitatory neurons)`
names(dict) <- Neurons$`Cell name`

#Read in the Seurat Files and convert to Seurat Object#
counts <- Read10X(data.dir = "GSE113576/")
seurat<-CreateSeuratObject(counts = counts, min.cells = 0, min.features = 500, project = "10X_MOF", meta.data = ID, assay = "RNA")

#Create the variables for cell and gene filters#
counts_per_cell <- Matrix::colSums(counts)
counts_per_gene <- Matrix::rowSums(counts)
genes_per_cell <- Matrix::colSums(counts > 0) # count gene only if it has non-zero reads mapped.
cell_per_gene <- Matrix::rowSums(counts > 0)

#Create gene filter - 50 reads or 20 individual cells#
Filter <- (counts_per_gene >= 50) + (cell_per_gene >= 20)
Filter <- Filter > 0

#Calculate percentage mito genes for each cell and set a logical variable for anything less than 10%#
mito.genes <- grep(pattern = "^mt", x = rownames(seurat@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(seurat@assays[["RNA"]][mito.genes, ])/Matrix::colSums(seurat@assays[["RNA"]])
mito <- percent.mito < 0.1

#Convert Seurat object into a matrix and then transformed data frame#
x <- GetAssayData(object = seurat, assay.type = "RNA", slot = "data")
x <- as(Class = 'matrix', object = x)
data <- data.frame(t(x))

#Create an RAW directory to save data so far#
ifelse(!dir.exists(file.path(getwd(), "RAW")), dir.create(file.path(getwd(), "RAW")), FALSE)

#Save data so far#
fwrite(data, "RAW/prefiltered_matrix.csv", quote = FALSE, row.names = TRUE)

#filter out mitochondrial genes#
data <- data[mito,]

#Create size factor vector for each cell to 10000 reads#
MOF_cell_Totals <- rowSums(data)
umi <- 10000 / MOF_cell_Totals

#Create col means vector#
MOF_gene_Averages <- colMeans(data)

#Make new dataframe#
Data <- data

#Loop through each cell and normalise each cell to 10000 reads and log transform the data#
for (j in 1:nrow(Data)){
  Data[j,] = Data[j,]*umi[j]
  Data[j,] = log(Data[j,] + 1)
}

#Apply gene filter to data#
Data <- Data[,Filter]

# deal with any duplicated gene names by adding _1 after duplicates#
if(sum(duplicated(rownames(Data))) != 0){ 
  dup <- which(duplicated(colnames(Data)))
  for(i in 1:length(dup)){
    colnames(Data)[dup[i]] <- paste(colnames(Data)[dup[i]], "_2", sep = "")
  }
}

#Create Neuron only data and genelist by filtering out non-neuron cells#
filt <- rownames(Data) %in% Neurons$`Cell name`
N <- Data[filt,]

#Read in the raw matrix to calculate neuron specific gene filters#
data <- fread("RAW/prefiltered_matrix.csv")
data$V1 <- str_replace(data$V1, "-", ".")
filt <- colnames(data) %in% colnames(N)
filt2 <- data$V1 %in% rownames(N)
data_neurons <- data[filt2, ..filt]

#Calculate the gene filter for the neuron only dataset#
counts_per_gene <- colSums(data_neurons)
cell_per_gene <- colSums(data_neurons>0)
cell = cell_per_gene >= 20
umi <- counts_per_gene >= 50
Filter <- umi + cell
Filter <- Filter > 0

#Apply the filter, create the cell column and save the neuronal dataset N#
N <- N[,..Filter]
N$iden <- dict[rownames(N)]
fwrite(N, "N.csv", row.names = TRUE)

#Remove the non-POA neuron types and save the POA specific neuron type as N2#
N2 <- N[!(N$iden %in% c("h3:Slc32a1/Gsc", "i14:Avp/Cck", "i6:Avp/Nms", "e20:Crh", "e21:Glut/Rxpf3", "i27:Th/Trh", "i28:Gaba/Six6", "i30: Vip", "i31:Calca")),]
fwrite(N2, "N2.csv", row.names = TRUE)

##Create the list of data to loop through and the names of the data to be used when creating folders#
data_list <- c("N", "N2")
data_names <- c("just_neurons", "just_POA_neurons")

#Create an Outputs directory to save pending analysis#
ifelse(!dir.exists(file.path(getwd(), "Outputs")), dir.create(file.path(getwd(), "Outputs")), FALSE)

#Loop through the data configuration available#
for(a in 1:length(data_names)){
  #Create wil object of the specific data to be used on this run
  wil <- data.frame(fread(paste(data_list[a], ".csv", sep = "")))
  
  #Name rownames by column V1 and delete it#
  rownames(wil) <- wil$V1
  wil$V1 = NULL

  #Ensure only complete cases are used#
  wil <- wil[complete.cases(wil), ]
  
  #Ensure cell identity colname is called "iden"#
  colnames(wil)[length(wil)] = "iden"
  
  #Create iden - a list of unique sorted identity types to use to organise analysis#
  iden <- sort(unique(wil$iden))
  
  #Create a specific folder in Outputs to save this data type's files to, and write a file with the number of cells for each identity type# 
  ifelse(!dir.exists(file.path(paste(getwd(), "/Outputs", sep = ""), data_names[a])), dir.create(file.path(paste(getwd(), "/Outputs", sep = ""), data_names[a])), FALSE)
  
  #write table of number of cells per cell identity#
  write.csv(table(wil$iden), paste("Outputs/", data_names[a], "/cell_numbers_Global.csv", sep = ""))
  
  #set dataframe p,q,fc,percent_in,percent_rest with each row a gene and each column an identity type#
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
  fwrite(p, paste("Outputs/", data_names[a], "/Moffitt_p.csv", sep=""))
  colnames(percent_in) <- c("gene", iden)
  fwrite(percent_in, paste("Outputs/", data_names[a], "/Moffitt_percent_in.csv", sep=""))
  colnames(percent_rest) <- c("gene", iden)
  fwrite(percent_rest, paste("Outputs/", data_names[a], "/Moffitt_percent_rest.csv", sep=""))
  
  #Change colnames to identity values#
  colnames(fc) <- c("gene", iden)
  #set infinite FC values to the max for that identity types#
  for(x in 2:ncol(fc)){
    m = max(fc[is.finite(fc[,x]),x])
    fc[is.infinite(fc[,x]),x] = m
  }
  
  #Save this file#
  fwrite(fc, paste("Outputs/", data_names[a], "/Moffitt_fc.csv", sep=""))
  
  # adjust p value horizontally using BH correcton and update this to the q dataframe#
  for (x in 1:nrow(p)){
    q[x,2:ncol(p)] <- p.adjust(as.numeric(p[x,2:ncol(p)]), method = "BH")
  }
  
  #Change colnames to identity values and save the file#
  colnames(q) <- c("gene", iden)
  fwrite(q, paste("Outputs/", data_names[a], "/Moffitt_q.csv", sep=""))
  
### STAGE 2 ENRICHMENT ANALYSIS ########################################################################################

  #Set two fold change limits for upregulated genes#  
  limit <- c(2,1)
  
  #loop over the number of FC limits#  
  for (r in 1:length(limit)){
    #create folder for that fc limit#
    ifelse(!dir.exists(file.path(paste(getwd(), "/Outputs/", data_names[a],"/",sep = ""), limit[r])), dir.create(file.path(paste(getwd(), "/Outputs/",data_names[a],"/", sep = ""), limit[r])), FALSE)
    
    #define dataframes for files to later save#
    Fish <- data.frame("Identity" = colnames(q[,-1]), stringsAsFactors = FALSE)
    Genelist <- data.frame("Identity" = colnames(q[,-1]), stringsAsFactors = FALSE)
    DEGs <- data.frame()
    IGs <- data.frame()
    
    #loop through identity groups (columns in the p/fc/q etc. files)#
    for (e in 2:(ncol(q))){
      ORA <- data.frame(gene = q$gene, p = p[,e], q = q[,e], fc = fc[,e], pct.in = percent_in[,e], pct.rest = percent_rest[,e], stringsAsFactors = FALSE)
      
      #filter ORA by significant q values and fc limit#
      ORA <- ORA[ORA$q <= 0.05 & ORA$fc >= limit[r],]
      
      #create repeat column of identity name#
      ORA$Identity = colnames(q)[e]
      
      #Add rows to DEGS - all differentialy expressed genes for that tissue and the total number to Fish#
      DEGs <- rbind(DEGs, ORA)
      Fish[e-1,2] = nrow(ORA)
      
      #Filter ORA for imprinted genes only and add that to IGs file and the total number to Fish#
      ig <- ORA %>% filter(gene %in% IG$Gene)
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
    fwrite(Genelist, paste("Outputs/",data_names[a],"/", limit[r], "/IG_Subpopulation_List.csv", sep =""))
    fwrite(DEGs, paste("Outputs/", data_names[a],"/", limit[r], "/Upregulated_Genes.csv", sep =""))
    fwrite(IGs, paste("Outputs/", data_names[a],"/", limit[r], "/Upregulated_IGs.csv", sep =""))
    
### OVER REPRESENTATION ANALYSIS (ORA) ########################################################################
    
    #Create an identity group filter based on having a minimum of 5 imprinted genes upregulated#
    tissue_ORA <- Fish$Identity[Fish$IG >= (as.numeric(sum(q$gene %in% IG$Gene))/20)]
    
    #As long as one cell identity has 5 or more IGs, run a Fisher's Exact test#  
    if(length(tissue_ORA) > 0){
      #For every identity, check if it is in the tissue filter#
      for (c in 1:nrow(Fish)){
        if(Fish$Identity[c] %in% tissue_ORA){
          #Create the 2x2 matrix for the fisher's exact test - no. IG in group, no. of IG left over, no. of UpRegulated genes in group minus no. of IGs, all genes left over
          test <- data.frame(gene.interest=c(as.numeric(Fish[c,3]),as.numeric(sum(q$gene %in% IG$Gene) - Fish[c,3])), gene.not.interest=c((as.numeric(Fish[c,2]) - as.numeric(Fish[c,3])),as.numeric(nrow(q) - Fish[c,2])- as.numeric(sum(q$gene %in% IG$Gene) - Fish[c,3])))
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
    
    #Calculate mean fold change for upregulated imprinted genes in each identity group#
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
    rest <- DEGs[!(DEGs$gene %in% IG$Gene),] %>% group_by(Identity) %>% summarise("mean" = mean(fc))
    
    #Merge Fish with rest to have a Mean Fold change Rest of genes column#
    Fish <- merge(Fish, rest, by = "Identity")
    names(Fish)[length(names(Fish))] <- "Mean FC Rest"
    
### GENE SET ENRICHMENT ANALYSIS (GSEA) ########################################################################  
    
    #Create list of identities that fulfill the criteria for GSEA - more than 10% of imprinted genes present and a IG FC higher than Rest FC#
    tissue_GSEA <- Fish$Identity[(Fish$`Mean FC IG` >= Fish$`Mean FC Rest` & Fish$IG >= 15)]
    
    ##If there are tissues in tissue_GSEA run a GSEA analysis for IGs in each of those tissues##
    if(length(tissue_GSEA) > 0){
      ifelse(!dir.exists(file.path(paste(getwd(), "/Outputs/",data_names[a],"/", limit[r], sep = ""), "GSEA")), dir.create(file.path(paste(getwd(), "/Outputs/",data_names[a],"/", limit[r], sep = ""), "GSEA")), FALSE)
      for(f in 1:length(tissue_GSEA)){
        #filter DEGs for the tissue selected#
        daa <- DEGs %>% filter(DEGs$Identity == tissue_GSEA[f])
        #create log transformed fc values for DEGs#
        file <- log2(daa$fc)
        #replace infinity values with max values and name the file with gene names#
        file[is.infinite(file)] <- max(file[is.finite(file)])
        file <- set_names(file, daa$gene)
        #create imprinted gene list to use for GSEA#
        ig <- as.character(IG$Gene)
        #save the graph from the GSEA as a pdf#
        pdf(paste("Outputs/",data_names[a],"/", limit[r],"/GSEA/",gsub("[^A-Za-z0-9 ]","",tissue_GSEA[f]), "_GSEA.pdf", sep=""))
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
    
### IMPRINTED GENE TOP TISSUE EXPRESSION ######################################################################    
    
    #Filter Main data for Imprinted Genes and create Identity column#
    D <- wil[,colnames(wil) %in% c(as.character(IG$Gene), "iden")] %>% gather("gene", "reads", -iden)
    #Gather this data and group by gene, identity and get per gene per identity read averages#
    D <- D %>% group_by(gene, iden) %>% summarise("avg" = mean(as.numeric(reads)))
    #Create matching file of fc's with just imprinted genes#
    FC <- fc %>% gather(iden, fc, -gene) %>% filter(gene %in% IG$Gene)
    #Create D2 - combination of D and FC ready for making dotplots#
    D2 <- left_join(D, FC, by = c("gene","iden"))
    
    #Create u with the tissue with max expression for each IG and save that file ##
    u <- D2 %>% group_by(gene) %>% filter(gene %in% IG$Gene & avg == max(avg))
    fwrite(u, paste("Outputs/", data_names[a], "/", limit[r], "/IG_Top_Expressed_Subpopulation.csv", sep =""))
    #update a Fish column with the number of IGs with top expression in that tissue#
    u <- u %>% group_by(iden) %>% summarise("Top_Tissue_IGs" = n())
    Fish <- merge(Fish, u, by.x = "Identity", by.y = "iden", all = TRUE)
    
    #write the finished Fish file#
    fwrite(Fish, paste("Outputs/", data_names[a], "/", limit[r], "/Enrichment_Analysis.csv", sep =""))
  }
  
  
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
  Pat$gene <- factor(Pat$gene, levels = rev(IG$Gene))
  Pat <- Pat %>% arrange(rev(gene))
  
  #Filter IGs for maternally expressed genes (MEGs) only and create Mat using the MEG only filter#
  mat <- IG2 %>% filter(Sex == "M" | Sex == "I")
  Mat <- D2 %>% filter(gene %in% mat$Gene)
  
  #Recast Gene as a Factor and arrange by chromosomal order#
  Mat$gene <- factor(Mat$gene, levels = rev(IG$Gene))
  Mat <- Mat %>% arrange(rev(gene))
  
  #Create PDF to save PEG dotplot#
  pdf(paste("Outputs/", data_names[a],"/", "PEG_DOTPLOT.pdf", sep=""))
  
  #GGplot dotplot, x = cell identity, y = gene identity, color = fc(gradated up to 5FC+), size = avg expression(0 to max expression registered)#
  print(ggplot(Pat, aes(x=iden, y=gene, color=ifelse(fc == 0, NA, fc), size=ifelse(avg==0, NA, avg))) + geom_point(alpha = 0.8) +
          theme_classic() +
          scale_color_gradientn(colours = c("grey95","grey60","blue","darkblue","midnightblue"), na.value="midnightblue", values = c(0, 0.2, 0.4, 0.6, 0.8 ,1), limits = c(0,5)) +
          theme(axis.text.x = element_text(angle = 90)) +
          scale_size_continuous(limits = c(0,max(D2$avg)))+
          scale_x_discrete(limits = order) +
          labs(size = "Normalised_Mean_Expression", color = "Proportion_Expression_vs_Mean"))
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
          labs(size = "Normalised_Mean_Expression", color = "Proportion_Expression_vs_Mean"))
  #Save PDF#
  dev.off()
}