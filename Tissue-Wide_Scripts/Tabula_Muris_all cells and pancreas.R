#### STAGE 0 - NECESSARY SETUP (RUN EVERY TIME) ######################################

#Prior to running - download 20 FACS Robjects from #
#https://figshare.com/articles/dataset/Robject_files_for_tissues_processed_by_Seurat/5821263#
#Place in Directory - "Facs_seurat"

#read in packages#library(Seurat)
library(data.table)
library(tidyverse)
library(liger)
library(cowplot)
library(ggnewscale)
library(Seurat)

#set directory#
setwd("...")

#"Imprinted_Gene_List.csv" - read in premade Imprinted gene file- gene name, ensmbl, chromosome order and sex_bias#
IG <- read.csv("...")

#Create an Outputs directory to save pending analysis#
ifelse(!dir.exists(file.path(getwd(), "Outputs")), dir.create(file.path(getwd(), "Outputs")), FALSE)

#Redirect to folder with data files and list them
setwd("Facs_seurat/")
files <- list.files()

#Load the first RAW COUNTS dataset and create a dataframe from it with gene names as columns#
load(files[1])
total <- data.frame(t(as.matrix(tiss@raw.data)))

#Loop over the subsequent datasets and transform them in the same way before row binding to make the large total dataset#
for(a in 2:length(files)){
  load(files[a])
  x <- data.frame(t(as.matrix(tiss@raw.data)))
  total <- rbind(total, x)
}

#Create the gene filters, genes expressed in 20 unique cells#
filter <- colSums(total != 0) >= 20

#Load the first file again and this time create the dataframe with the transformed PROCESSED dataset#
load(files[1])
total <- data.frame(t(as.matrix(tiss@data)))

#Loop over the subsequent PROCESSED datasets and transform them in the same way before row binding to make the large total dataset#
for(a in 2:length(files)){
  load(files[a])
  x <- data.frame(t(as.matrix(tiss@data)))
  total <- rbind(total, x)
}

#Apply the previously generated gene filter on the total processed dataset#
total <- total[,filter]

#Save this datafile to retrieve later#
fwrite(total, "../Outputs/merged.csv", row.names = TRUE)

#Leave the raw data directory#
setwd("../")

#Read in the cell annotations#
ID <- read.csv("https://ndownloader.figshare.com/files/13088129")

#Create id dataframe with all the levels of cell identification on it#
id <- data.frame(iden = ID$cell,cell = as.character(ID$cell_ontology_class), free = as.character(ID$free_annotation), subtissue = as.character(ID$subtissue), tissue = as.character(ID$tissue), stringsAsFactors = FALSE)

#Loop through each cell and assign stand in identities where a cell doesnt have one e.g. the whole tissue identity to the sub tissue identity if that tissue doesnt have a subtissue one#
for(a in 1:nrow(id)){
  if(id$cell[a] == "" & id$free[a] == ""){
    id$cell[a] = id$subtissue[a]
  }else if(id$cell[a] == ""){
    id$cell[a] = id$free[a]
  }
  id$subtissue[a] <- paste(id$tissue[a], id$subtissue[a])
}

#Create the 2 dictionaries to assign cell barcodes with tissue and unique cell identities#
dict_tissue <- id$tissue
dict_cell <- id$cell
names(dict_tissue) <- id$iden
names(dict_cell) <- id$iden

#Create the lists of data and relevant dictionaries to loop through#
data_names <- c("cells", "Pancreas")
limit <- c("2", "1")

#Loop through each dataset#
for (a in 1:length(data_names)){
  #Assign the main data to wil
  wil <- total
  #Give that data an iden column with the relevant cell identities#
  wil$iden <-  dict_cell[rownames(wil)]
  
  if(data_names[a] != "cells"){
    wil <- wil[dict_tissue[rownames(wil)] == data_names[a],]
    wil <- wil[,colSums(wil != 0) >= 20]
  }
  
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
  fwrite(p, paste("Outputs/", data_names[a], "/tabulamuris_p.csv", sep=""))
  colnames(percent_in) <- c("gene", iden)
  fwrite(percent_in, paste("Outputs/", data_names[a], "/tabulamuris_percent_in.csv", sep=""))
  colnames(percent_rest) <- c("gene", iden)
  fwrite(percent_rest, paste("Outputs/", data_names[a], "/tabulamuris_percent_rest.csv", sep=""))
  
  #Change colnames to identity values#
  colnames(fc) <- c("gene", iden)
  #set infinite FC values to the max for that identity types#
  for(x in 2:ncol(fc)){
    m = max(fc[is.finite(fc[,x]),x])
    fc[is.infinite(fc[,x]),x] = m
  }
  
  #Save this file#
  fwrite(fc, paste("Outputs/", data_names[a], "/tabulamuris_fc.csv", sep=""))
  
  # adjust p value horizontally using BH correcton and update this to the q dataframe#
  for (x in 1:nrow(p)){
    q[x,2:ncol(p)] <- p.adjust(as.numeric(p[x,2:ncol(p)]), method = "BH")
  }
  
  #Change colnames to identity values and save the file#
  colnames(q) <- c("gene", iden)
  fwrite(q, paste("Outputs/", data_names[a], "/tabulamuris_q.csv", sep=""))
  
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
    ORA <- ORA[ORA$q <= 0.05 & ORA$fc >= limit[a],]
    
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
  fwrite(Genelist, paste("Outputs/", data_names[a],"/IG_Subpopulation_List.csv", sep =""))
  fwrite(DEGs, paste("Outputs/", data_names[a], "/Upregulated_Genes.csv", sep =""))
  fwrite(IGs, paste("Outputs/", data_names[a], "/Upregulated_IGs.csv", sep =""))
  
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
  rest <- DEGs[!(DEGs$gene %in% IG$Gene),] %>% group_by(Identity) %>% summarise("mean" = mean(fc))
  
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
      file <- set_names(file, daa$gene)
      #create imprinted gene list to use for GSEA#
      ig <- as.character(IG$Gene)
      #save the graph from the GSEA as a pdf#
      #pdf(paste("Outputs/", data_names[a], "/GSEA/",gsub("[^A-Za-z0-9 ]","",tissue_GSEA[f]), "_GSEA.pdf", sep=""))
      #Add GSEA p value to Fish#
      Fish[Fish$Identity == tissue_GSEA[f],8] <- gsea(file, ig)
      #Save PDF#
      #dev.off()
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
  
  fwrite(Fish, paste("Outputs/", data_names[a], "/Enrichment_Analysis.csv", sep =""))


  #### STAGE 3 - VISUALISATION DOTPLOT ########################################################################################################
  
  #Arrange Fish by over-representation significance, take the identity variable as an order value#
  Fish <- Fish %>% arrange(Fish$ORA_p)
  order <- c("pancreatic D cell", "pancreatic A cell", "type B pancreatic cell", "pancreatic PP cell", "pancreatic stellate cell", "endothelial cell", "pancreatic ductal cell", "pancreatic acinar cell", "leukocyte")
  
  #Filter original IG list with those in the dataset#
  IG2 <- IG[IG$Gene %in%  IGs$gene[IGs$Identity %in% c("pancreatic D cell", "pancreatic A cell", "type B pancreatic cell	", "pancreatic PP cell")],]
  #Arrange IGs by the chromosomal order#
  IG2$Gene <- as.character(IG2$Gene)
  IG2 <- IG2 %>% arrange(IG2$Order)
  
  #Filter IGs for Paternally expressed genes (PEGs) only and create Pat using the PEG only filter#
  pat <- IG2 %>% filter(Sex == "P")
  Pat <- D2 %>% filter(gene %in% pat$Gene)
  Pat$fc <- log2(Pat$fc+1)
  
  #Filter IGs for maternally expressed genes (MEGs) only and create Mat using the MEG only filter#
  mat <- IG2 %>% filter(Sex == "M" | Sex == "I")
  Mat <- D2 %>% filter(gene %in% mat$Gene)
  Mat$fc <- log2(Mat$fc+1)
  
  ### Create Function to align the legend of the dotplot to the centre ###
  align_legend <- function(p, hjust = 0.5)
  {
    # extract legend
    g <- cowplot::plot_to_gtable(p)
    grobs <- g$grobs
    legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
    legend <- grobs[[legend_index]]
    
    # extract guides table
    guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
    
    # there can be multiple guides within one legend box  
    for (gi in guides_index) {
      guides <- legend$grobs[[gi]]
      
      # add extra column for spacing
      # guides$width[5] is the extra spacing from the end of the legend text
      # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
      # both sides, we get an aligned legend
      spacing <- guides$width[5]
      guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
      guides$widths[6] <- (1-hjust)*spacing
      title_index <- guides$layout$name == "title"
      guides$layout$l[title_index] <- 2
      
      # reconstruct guides and write back
      legend$grobs[[gi]] <- guides
    }
    
    # reconstruct legend and write back
    g$grobs[[legend_index]] <- legend
    g
  }
  
  dot <- ggplot(D2, aes(x=iden, y=gene, color=ifelse(fc == 0, NA, fc), size=ifelse(avg==0, NA, avg))) + 
    geom_point(data = Pat,aes(color = ifelse(fc == 0, NA, fc)), alpha = 0.8) +
    scale_color_gradientn(colours = c("grey95","grey90","blue","darkblue","midnightblue"), na.value="midnightblue", values = c(0, 0.2, 0.4, 0.6, 0.8 ,1), limits = c(0,4)) +
    labs(color = "Log2FC vs\nBackground\n(PEGs)")+
    new_scale_color()+
    geom_point(data = Mat, aes(color = ifelse(fc == 0, NA, fc)), alpha = 0.8) +
    scale_color_gradientn(colours = c("grey95","grey90","orange","red","darkred"), na.value="darkred", values = c(0, 0.2, 0.4, 0.6, 0.8 ,1), limits = c(0,4)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, face = "italic"), axis.title.x = element_text(angle = 180), axis.title.y.right = element_text(angle = 90), axis.text.y = element_text(angle = 180), legend.title.align=0.5) +
    scale_size_continuous(limits = c(0,max(D2$avg)))+
    scale_x_discrete(limits = order, position = "top") +
    scale_y_discrete(limits = (IG2$Gene))+
    guides(size = guide_legend(order = 1))+
    coord_flip()+
    labs(x = "Cell Identity (Pancreas - TM)", y = "Imprinted Gene", size = "Normalised\nMean\nExpression", color = "Log2FC vs\nBackground\n(MEGs)")
  
  pdf("Outputs/TM_Pancreas_Cell_lineage_DOTPLOT.pdf", paper= "a4r",  width = 28, height = 18)
  ggdraw(align_legend(dot))
  dev.off()
}    
