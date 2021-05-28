#### STAGE 0 - NECESSARY SETUP (RUN EVERY TIME) ###############################################################

#read in packages#
library(loomR)
library(data.table)
library(cowplot)
library(tidyverse)
library(data.table)
library(liger)

#set working directory/ignore if running from current directory#
setwd(...)

#"Imprinted_Gene_List.csv" - read in premade Imprinted gene file- gene name, ensmbl, chromosome order and sex_bias#
IG <- read.csv(...)

#Download "l5.all.loom file" from https://storage.googleapis.com/linnarsson-lab-loom/l5_all.loom#
#And save into "Data/" directory#

#Create an Outputs directory to save pending analysis#
ifelse(!dir.exists(file.path(getwd(), "Outputs")), dir.create(file.path(getwd(), "Outputs")), FALSE)

### STAGE 1 DATA PREPARATION AND UPREGULATION ANALYSIS ###################################################################################

## Connect to the lfile with ALL CELLS in ##
lfile <- connect(filename = "Data/l5_all.loom", mode = "r+")

#Create various objects with cell identity labels to use later#
class_whole <- lfile$col.attrs$Class[]
tissue_whole <- lfile$col.attrs$Tissue[]
subclass_whole <- lfile$col.attrs$ClusterName[]
tissue_neurons <- lfile$col.attrs$Tissue[lfile$col.attrs$Class[] == "Neurons"]
subclass_neurons <- lfile$col.attrs$ClusterName[lfile$col.attrs$Class[] == "Neurons"]

#Extract gene short names to relabel ensembl IDs with#
genes <- lfile$row.attrs$Accession[]
tag <- lfile$row.attrs$Gene[]
names(tag) <- genes

#Preset vectors to record no. of transcript and no. of cells per gene for ALL CELLS and for ONLY NEURONS#
transcripts = vector()
percell = vector()
transcripts_neuron = vector()
percell_neuron = vector()

#Loop over the datamatrix per gene and append to the vectors above with the specific values#
for (a in 1:length(lfile$row.attrs$`_Total`[])){
  #get data for one gene only#
  data.gene <- lfile[["matrix"]][,a]
  #Append total transcript number#
  transcripts <- append(transcripts, sum(data.gene))
  #Append total cell number#
  percell <- append(percell, sum(data.gene != 0))
  #Limit the data to only neurons#
  data.gene <- lfile[["matrix"]][lfile$col.attrs$Class[] == "Neurons",a]
  #Append total transcript number#
  transcripts_neuron <- append(transcripts_neuron, sum(data.gene))
  #Append total cell number#
  percell_neuron <- append(percell_neuron, sum(data.gene != 0))
  #Print loop number as tracker#
  print(a)
}

#Create logical gene filter for all cells and for just neurons
filter <- percell >= 20
filter_neuron <- percell_neuron >= 20

#Create numeric vector to indicate valid genes#
Fil <- which(filter, arr.ind = FALSE)
Fil_Neuron <- which(filter_neuron, arr.ind = FALSE)

#Create scaling factor to 5000 reads#
Scale <- 5000/lfile$col.attrs$`_Total`[]

#Create list of datasets to put through basis workflow#
data_names <- c("class_whole", "tissue_neurons", "subclass_neurons")

#Loop through the data configuration available#
for (a in 1:length(data_names)){
  #If data includes all available cells#
  if(a == 1){
    #create gene list with valid genes
    gene <- genes[filter]
    #else if the data is neuronal use the neuronally valid genes#
  }else{gene <- genes[filter_neuron]}
  
  #Create a specific folder in Outputs to save this data type's files to#
  ifelse(!dir.exists(file.path(paste(getwd(), "/Outputs", sep = ""), data_names[a])), dir.create(file.path(paste(getwd(), "/Outputs", sep = ""), data_names[a])), FALSE)
  
  #Create iden - a list of unique sorted identity types to use to organise analysis#
  iden <- sort(unique(eval(parse(text = data_names[a]))))
  
  #write table of number of cells per cell identity#
  write.csv(table(eval(parse(text = data_names[a]))), paste("Outputs/", data_names[a], "/cell_numbers_Global.csv", sep = ""))
  
  #set dataframe p,q,fc,percent_in,percent_rest with each row a gene and each column an identity type#
  p <- data.frame("ensmbl" =gene, stringsAsFactors = FALSE)
  fc <- data.frame("ensmbl" = gene, stringsAsFactors = FALSE)
  percent_in <- data.frame("ensmbl" = gene, stringsAsFactors = FALSE)
  percent_rest <- data.frame("ensmbl" = gene, stringsAsFactors = FALSE)
  q <- data.frame("ensmbl" = gene, stringsAsFactors = FALSE)
  
  #Loop over every gene/column#
  for (b in 1:length(gene)){
    #If data is not neuron specific
    if(a == 1){
      #Extract the data for just that gene across all cells#
      data.gene <- lfile[["matrix"]][,lfile$row.attrs$Accession[] == gene[b]]
      #Scale and normalise data#
      data.gene <- log((data.gene * Scale)+1)
    } else{
      #Extract the data for just that gene and restrict cells to only neurons#
      data.gene <- lfile[["matrix"]][lfile$col.attrs$Class[] == "Neurons" ,lfile$row.attrs$Accession[] == gene[b]]
      #Scale and normalise data#
      data.gene <- log((data.gene * Scale[lfile$col.attrs$Class[] == "Neurons"])+1)
    }
    
    #Convert that data into a dataframe#
    dd <- data.frame("reads" = data.gene[])
    
    #Create iden column of cell identities based on relevant dictionary#
    dd$iden <- eval(parse(text = data_names[a]))
    
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
  colnames(p) <- c("ensmbl", iden)
  p$gene <- tag[as.character(p$ensmbl)]
  fwrite(p, paste("Outputs/", data_names[a], "/Zeisel_p.csv", sep=""))
  colnames(fc) <- c("ensmbl", iden)
  for(x in 2:ncol(fc)){
    m = max(fc[is.finite(fc[,x]),x])
    fc[is.infinite(fc[,x]),x] = m
  }
  fc$gene <- tag[as.character(fc$ensmbl)]
  fwrite(fc, paste("Outputs/", data_names[a], "/Zeisel_fc.csv", sep=""))
  colnames(percent_in) <- c("ensmbl", iden)
  percent_in$gene <- tag[as.character(percent_in$ensmbl)]
  fwrite(percent_in, paste("Outputs/", data_names[a], "/Zeisel_percent_in.csv", sep=""))
  colnames(percent_rest) <- c("ensmbl", iden)
  percent_rest$gene <- tag[as.character(percent_rest$ensmbl)]
  fwrite(percent_rest, paste("Outputs/", data_names[a], "/Zeisel_percent_rest.csv", sep=""))
  
  # adjust p value horizontally #
  

  for (x in 1:nrow(p)){
    q[x,2:(ncol(p)-1)] <- p.adjust(as.numeric(p[x,2:(ncol(p)-1),]), method = "BH")
  }
  colnames(q) <- c("ensmbl", iden)
  q$gene <- tag[as.character(q$ensmbl)]
  fwrite(q, paste("Outputs/", data_names[a], "/Zeisel_q.csv", sep=""))
  
  #### ORA and GSEA ANALYSIS ################################################################################################  
  
  #Set up the paretal origin nature of the IGs and the FC limits to loop over#
  if(a %in% c(1,2)){
    sex <- c("P","M","ALL")
  }else{ sex <- c("ALL")}

  #Loop over PEGs MEGS and all IGs#
  for (s in c(1:length(sex))){
    #set a second imprinted gene list to subset#
    IIGG <- IG
    #For Pegs and Megs reduce the list#
    if(sex[s] %in% c("P","M")){
      IIGG <- IIGG[IIGG$Sex == sex[s],]
    }
    
    #Set up a directory for the sex subset#
    ifelse(!dir.exists(file.path(paste(getwd(), "/Outputs/",data_names[a], sep = ""), sex[s])), dir.create(file.path(paste(getwd(), "/Outputs/",data_names[a], sep = ""), sex[s])), FALSE)
    
    #define dataframes for files to later save#
    Fish <- data.frame("Identity" = colnames(q[,c(-1, -ncol(q))]), stringsAsFactors = FALSE)
    Genelist <- data.frame("Identity" = colnames(q[,c(-1, -ncol(q))]), stringsAsFactors = FALSE)
    DEGs <- data.frame()
    IGs <- data.frame()
    
    #loop through identity groups (columns in the p/fc/q etc. files)#
    for (e in 2:(ncol(q)-1)){
      ORA <- data.frame(gene = q$ensmbl, p = p[,e], q = q[,e], fc = fc[,e], pct.in = percent_in[,e], pct.rest = percent_rest[,e], stringsAsFactors = FALSE)
      
      #filter ORA by significant q values and fc limit#
      ORA <- ORA[ORA$q <= 0.05 & ORA$fc >= 2,]
      
      #create repeat column of identity name#
      ORA$Identity = colnames(q)[e]
      
      #Add rows to DEGS - all differentialy expressed genes for that tissue and the total number to Fish#
      DEGs <- rbind(DEGs, ORA)
      Fish[e-1,2] = nrow(ORA)
      
      #Filter ORA for imprinted genes only and add that to IGs file and the total number to Fish#
      ig <- ORA %>% filter(gene %in%  IIGG$ï..Ensembl)
      ig$gene <- tag[as.character(ig$gene)]
      IGs <- rbind(IGs, ig)
      Fish[e-1,3] = nrow(ig)
      
      #Add to the genelists for IGs and Marker genes only if there are significant genes for that tissue#
      if (nrow(ig) > 0){
        for (z in 1:nrow(ig)){
          Genelist[e-1, z+1] <- as.character(ig$gene[z])
        }
      }
    }
    
    #add colnames to Fish and save all the files generated other than Fish#
    colnames(Fish)<- c("Identity","UpReg", "IG")
    fwrite(Genelist, paste("Outputs/", data_names[a],"/", sex[s], "/IG_Subpopulation_List.csv", sep =""))
    DEGs$gene_name <- tag[as.character(DEGs$gene)]
    fwrite(DEGs, paste("Outputs/", data_names[a],"/", sex[s], "/Upregulated_Genes.csv", sep =""))
    fwrite(IGs, paste("Outputs/", data_names[a], "/", sex[s], "/Upregulated_IGs.csv", sep =""))
    
    ### OVER REPRESENTATION ANALYSIS (ORA) ########################################################################
    
    #Create an identity group filter based on having a minimum of 5 imprinted genes upregulated#
    tissue_ORA <- Fish$Identity[Fish$IG >= (as.numeric(sum(q$ensmbl %in% IIGG$ï..Ensembl))/20)]
    
    #As long as one cell identity has 5 or more IGs, run a Fisher's Exact test#  
    if(length(tissue_ORA) > 0){
      #For every identity, check if it is in the tissue filter#
      for (c in 1:nrow(Fish)){
        if(Fish$Identity[c] %in% tissue_ORA){
          #Create the 2x2 matrix for the fisher's exact test - no. IG in group, no. of IG left over, no. of UpRegulated genes in group minus no. of IGs, all genes left over
          test <- data.frame(gene.interest=c(as.numeric(Fish[c,3]),as.numeric(sum(q$ensmbl %in% IIGG$ï..Ensembl) - Fish[c,3])), gene.not.interest=c((as.numeric(Fish[c,2]) - as.numeric(Fish[c,3])),as.numeric(nrow(q) - Fish[c,2])- as.numeric(sum(q$ensmbl %in% IIGG$ï..Ensembl) - Fish[c,3])))
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
    rest <- DEGs[!(DEGs$gene %in% IIGG$ï..Ensembl),] %>% group_by(Identity) %>% summarise("mean" = mean(fc))
    
    #Merge Fish with rest to have a Mean Fold change Rest of genes column#
    Fish <- merge(Fish, rest, by = "Identity")
    names(Fish)[length(names(Fish))] <- "Mean FC Rest"
    
    ### GENE SET ENRICHMENT ANALYSIS (GSEA) ########################################################################  
    
    #Create list of identities that fulfill the criteria for GSEA - more than 10% of imprinted genes present and a IG FC higher than Rest FC#
    tissue_GSEA <- Fish$Identity[(Fish$`Mean FC IG` >= Fish$`Mean FC Rest` & Fish$IG >= 15)]
    
    ##If there are tissues in tissue_GSEA run a GSEA analysis for IGs in each of those tissues##
    if(length(tissue_GSEA) > 0){
      ifelse(!dir.exists(file.path(paste(getwd(), "/Outputs/", data_names[a], "/", sex[s], sep = ""), "GSEA")), dir.create(file.path(paste(getwd(), "/Outputs/", data_names[a], "/", sex[s], sep = ""), "GSEA")), FALSE)
      for(f in 1:length(tissue_GSEA)){
        #filter DEGs for the tissue selected#
        daa <- DEGs %>% filter(DEGs$Identity == tissue_GSEA[f])
        #create log transformed fc values for DEGs#
        file <- log2(daa$fc)
        #replace infinity values with max values and name the file with gene names#
        file[is.infinite(file)] <- max(file[is.finite(file)])
        file <- set_names(file, daa$gene)
        #create imprinted gene list to use for GSEA#
        ig <- as.character(IIGG$ï..Ensembl)
        #save the graph from the GSEA as a pdf#
        pdf(paste("Outputs/", data_names[a], "/", sex[s], "/GSEA/",gsub("[^A-Za-z0-9 ]","",tissue_GSEA[f]), "_GSEA.pdf", sep=""))
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
    if(a == 1){
      D <-  data.frame(lfile[["matrix"]][, lfile$row.attrs$Accession[] %in% q$ensmbl[q$ensmbl %in% IIGG$ï..Ensembl]])
      D <- log((D * Scale)+1)
    }else{ 
      D <- data.frame(lfile[["matrix"]][lfile$col.attrs$Class[] == "Neurons", lfile$row.attrs$Accession[] %in% q$ensmbl[q$ensmbl %in% IIGG$ï..Ensembl]])
      D <- log((D * Scale[lfile$col.attrs$Class[] == "Neurons"])+1)
    }
    
    #Rename colnames with gene names and ad the iden column#
    colnames(D) <- lfile$row.attrs$Gene[lfile$row.attrs$Accession[] %in% q$ensmbl[q$ensmbl %in% IIGG$ï..Ensembl]]
    D$iden <- eval(parse(text = data_names[a]))
    #Gather this data and group by gene, identity and get per gene per identity read averages#
    D <- D  %>% gather("gene", "reads", -iden)
    D <- D %>% group_by(gene, iden) %>% summarise("avg" = mean(as.numeric(reads)))
    #Create matching file of fc's with just imprinted genes#
    FC <- fc %>% filter(fc$ensmbl %in% IIGG$ï..Ensembl)
    FC <- FC[,-1] %>% gather(iden, fc, -gene) 
    #Create D2 - combination of D and FC ready for making dotplots#
    D2 <- left_join(D, FC, by = c("gene","iden"))
    
    #Create u with the tissue with max expression for each IG and save that file ##
    u <- D2 %>% group_by(gene) %>% filter(gene %in% IIGG$Gene & avg == max(avg))
    fwrite(u, paste("Outputs/", data_names[a], "/", sex[s], "/IG_Top_Expressed_Subpopulation.csv", sep =""))
    #update a Fish column with the number of IGs with top expression in that tissue#
    u <- u %>% group_by(iden) %>% summarise("Top_Tissue_IGs" = n())
    Fish <- merge(Fish, u, by.x = "Identity", by.y = "iden", all = TRUE)
    
    #write the finished Fish file#
    fwrite(Fish, paste("Outputs/", data_names[a], "/", sex[s], "/Enrichment_Analysis.csv", sep =""))
  }  
}     
#### STAGE 3 - VISUALISATION DOTPLOT ########################################################################################################
  
#Arrange Fish by GSEA significance, take the identity variable as an order value#

#Filter original IG list with those in the dataset#
if(a == 1){
  Fish <- Fish %>% arrange(Fish$ORA_p)
  order <- Fish$Identity
  IG2 <- IIGG[IIGG$Gene %in% IGs$gene[IGs$Identity == "Neurons"],]
}else if(a == 2){
  Fish <- Fish %>% arrange(Fish$ORA_p)
  order <- Fish$Identity
  IG2 <- IIGG[IIGG$Gene %in% IGs$gene[IGs$Identity %in% c("Hypoth", "Medulla", "Pons", "Mbv")],]
}else{
  Fish <- Fish %>% arrange(Fish$GSEA_p)
  order <- Fish$Identity[Fish$ORA_q <= 0.05]
  order <- order[!is.na(order)]
  IG2 <- IIGG[IIGG$Gene %in% IGs$gene[IGs$Identity %in% c("HBSER2", "HBSER4", "HBSER5")],]
}
#Arrange IGs by the chromosomal order#
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
  labs(x = "Cell Subpopulation (Mouse Brain Atlas)", y = "Imprinted Gene", size = "Normalised\nMean\nExpression", color = "Log2FC vs\nBackground\n(MEGs)")

pdf(paste("Outputs/", data_names[a], "/MBA_DOTPLOT.pdf", sep=""), paper= "a4r",  width = 28, height = 18)
ggdraw(align_legend(dot))
dev.off()
