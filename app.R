######## To make a PCA on which I can highlight specified gene expression
### Use Salmon-limma logCPM data and limma-generated MDS dimensional data for this exercise
library(ggplot2)
library(tidyverse)
library(lubridate)
library(shiny)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(DT)

# Load a list of human and mouse genes
list <- read.delim("genelist_human_mouse.txt", header=FALSE)
list <- as.vector(list$V1)

# Load and prep in-house human GC datasets
atac_vs_cpm <- read_excel("inhouse/logatac_vs_logcpm.xlsx")
mds <- read.delim("inhouse/reduced_salmon_limma_mds.xls") # keep only desired columns for plotting
lcpm <- read.delim("inhouse/salmon_limma_lcpm.xls")
lcpm$gene <- row.names(lcpm)
kgn_count <- read.delim("inhouse/kgn_logtpm.xls", header=TRUE)
kgn_volcano <- read_excel("inhouse/kgn_volcano.xlsx")

# Load in-house mouse GC datasets
mouse <- read_excel("inhouse/mouse_rnaseq.xlsx")

# Load timecourse datasets
GSE133868 <- read.delim("timecourse/GSE133868_rma_box.xls", header=TRUE)
GSE178314 <- read.delim("timecourse/GSE178314_logCPM_box.xls", header=TRUE)
GSE15315 <- read.delim("timecourse/GSE15315_rma_dot.xls", header=TRUE)
GSE44651 <- read.delim("timecourse/GSE44651_rma_box.xls", header=TRUE)
GSE20466 <- read.delim("timecourse/GSE20466_rma_box.xls", header=TRUE)
GSE119508 <- read.delim("timecourse/GSE119508_logCPM_box.xls", header=TRUE)
GSE140371 <- read.delim("timecourse/GSE140371_logCPM_box.xls", header=TRUE)
GSE167939 <- read.delim("timecourse/GSE167939_FKPM_dot.xls", header=TRUE)
GSE23084 <- read.delim("KO/GSE23084_rma_dot.xls", header=TRUE)

# Load KO datasets
GSE140922 <- read.delim("KO/GSE140922_logCPM_box.xls", header=TRUE)
GSE168213 <- read.delim("KO/GSE168213_logCPM_box.xls", header=TRUE)

# Load treatment datasets
GSE137608 <- read.delim("treatment/GSE137608_logCPM_dot.xls", header=TRUE)
GSE152727 <- read.delim("treatment/GSE152727_logCPM_box.xls", header=TRUE)
GSE158218_human <- read.delim("treatment/GSE158218_KGN_logCPM_box.xls", header=TRUE)
GSE158218_mouse <- read.delim("treatment/GSE158218_mouse_logCPM_box.xls", header=TRUE)
GSE166443 <- read.delim("treatment/GSE166443_logCPM_box.xls", header=TRUE)

# Load PCOS datasets
GSE34526 <- read.delim("PCOS/GSE34526_rma_pca.xls", header=TRUE)
GSE62093 <- read.delim("PCOS/GSE62093_logFPKM_pca.xls", header=TRUE)
GSE80432 <- read.delim("PCOS/GSE80432_rma_pca.xls", header=TRUE)
GSE110924 <- read.delim("PCOS/GSE110924_normalised_pca.xls", header=TRUE)
GSE114419 <- read.delim("PCOS/GSE114419_rma_pca.xls", header=TRUE)
GSE138518 <- read.delim("PCOS/GSE138518_logCPM_pca.xls", header=TRUE)
GSE168214 <- read.delim("PCOS/GSE168214_logCPM_box_pca.xls", header=TRUE)
GSE168404 <- read.delim("PCOS/GSE168404_logcount_pca.xls", header=TRUE)
GSE193123 <- read.delim("PCOS/GSE193123_logCPM_pca.xls", header=TRUE)

# Load metadata tables
table1 <- read_excel("table1.xlsx")
table2 <- read_excel("table2.xlsx")

# UI updated for faster gene lookup!
ui <- fluidPage(
  titlePanel("Gene expression in ovulatory granulosa cells"),
  sidebarLayout(
    sidebarPanel(
      selectizeInput("pickme", "Give me a (human format) gene", choices = NULL),
      br(),
      p("This is a compilation of gene expression data from human and mouse granulosa cells from the Russell lab and publicly available data. The aim is to have all the available granulosa cell data in a single user-friendly space."),
      p("This version contains bulk RNA-seq and some microarray datasets. scRNA-seq to follow once I've got the hang of it. Further details on the included datasets are available in the final tab."),
      br(),
      p("Mini Q&A:"),
      p("Can't find your gene? - The genelist was curated in favour of human nomenclature, so please make sure that the genename is accurate for human."),
      p("Gene not found in human or mouse datasets? - The gene likely has no human-mouse ortholog, or has been removed due to low expression."),
      p("Why are there different versions of the same gene? - In RNA-seq datasets, a gene can be represented by more than one transcript (isoform). In microarray, a gene can be targeted by multiple probes. Whenever possible, all transcripts or probes annotated to the gene are shown."),
      p("Where are all the stats? - This app has been designed to give a quick overview of the expression profiles across multiple different contexts, thus stats have not been added. That being said, all possible pairwise comparisons have been performed and data is avail upon request."),
      p("One of the plot isn't showing while the rest are fine? - Most likely because the selected gene is not found in the dataset for the missing plot. If this persists through different genes (esp commonly expressed genes), please notify Thao!"),
      p("Why are the units different between plots? logCPM/TPM/FPKM/RMA count are obtained through different normalisation formulas. These units are not interchangable, so do not compare the results ACROSS different plots."),
      p("See Thao for any other enquiries, bug reporting or data request :)"),
      p("Have fun exploring!")
    ),
    mainPanel(tabsetPanel(type = "tabs",
                          tabPanel("Human GC (inhouse)", 
                                   p("The following data comes from in-house human granulosa cell RNA-seq data."),
                                   br(),
                                   span(textOutput(outputId = "text_human1"), style = "color:red; font-size:20px; font-family:arial; font-style:italic"),
                                   column(6, plotOutput("plot_inhouse1")),
                                   column(6, plotOutput("plot_inhouse2")),
                                   p("MDS plot (left) shows the unsupervised clustering of all human GC samples and with the average gene expression level (logCPM)."),
                                   p("Scatter plot (left) shows the average gene expression level (logCPM) vs promoter chromatin accessibility (logATAC) in human GC across all samples."),
                                   br(),
                                   plotOutput("plot_inhouse3",height = 680),
                                   p("Bar graph shows the ranked gene expression (logCPM) for each library. Associated clinical data and fertility outcomes are listed in colour banners.")
                          ),
                          tabPanel("KGN (inhouse)",
                                   p("The following data comes from in-house human granulosa cell RNA-seq data."),
                                   br(),
                                   column(12, plotOutput("plot_KGN1")),
                                   p("Boxplot of gene expression (logTPM) for KGN overexpressed with PRA or PRB, treated with R5020 for 0h/8h/24h."),
                                   br(),
                                   column(6, plotOutput("plot_KGN2")),
                                   column(6, plotOutput("plot_KGN3")),
                                   column(6, plotOutput("plot_KGN4")),
                                   column(6, plotOutput("plot_KGN5")),
                                   p("Volcano plot with highlighted gene, from left to right: PRA 8h vs 0h, PRA 24h vs 0h (top plots), PRB 8h vs 0h, PRB 24h vs 0h (bottom plots).")
                          ),
                          tabPanel("Mouse GC LH / PRKO (inhouse)",
                                   p("The following data comes from in-house mouse granulosa cell RNA-seq."),
                                   br(),
                                   span(textOutput(outputId = "text_mouse1"), style = "color:red; font-size:20px; font-family:arial; font-style:italic"),
                                   column(6, plotOutput("plot_inhouse4")),
                                   column(6, plotOutput("plot_inhouse5")),
                                   p("Volcano plot with highlighted gene, from left to right: GC in vivo 8h vs 0h post-hCG, PGRKO vs PGRWT GC in vivo 8h post-hCG.")
                          ),
                          tabPanel("Ovulatory timecourse",
                                   p("See info tab for details of these datasets"),
                                   p("Boxplot (multiple replicates) or dotplot (single replicate) of gene expression (RMA count for microarray, logCPM or FPKM for RNA-seq)."),
                                   br(),
                                   column(3, plotOutput("plot_timecourse1")),
                                   column(3, plotOutput("plot_timecourse2")),
                                   column(3, plotOutput("plot_timecourse3")),
                                   column(3, plotOutput("plot_timecourse4")),
                                   column(3, plotOutput("plot_timecourse5")),
                                   column(3, plotOutput("plot_timecourse6")),
                                   column(3, plotOutput("plot_timecourse7")),
                                   column(3, plotOutput("plot_timecourse8")),
                                   column(3, plotOutput("plot_timecourse9"))
                                   ),
                          tabPanel("KO mouse models",
                                   p("See info tab for details of these datasets"),
                                   p("Boxplot (multiple replicates) or dotplot (single replicate) of gene expression (RMA count for microarray, logCPM for RNA-seq)."),
                                   br(),
                                   span(textOutput(outputId = "text_mouse_timecourse"), style = "color:red; font-size:20px; font-family:arial; font-style:italic"),
                                   column(3, plotOutput("plot_KO1")),
                                   column(3, plotOutput("plot_KO2")),
                                   column(3, plotOutput("plot_KO3")),
                                   column(3, plotOutput("plot_KO4")),
                                   column(3, plotOutput("plot_KO5")),
                                   column(3, plotOutput("plot_KO6")),
                                   column(3, plotOutput("plot_KO7"))
                          ),
                          tabPanel("Cell culture treatment",
                                   p("See info tab for details of these datasets"),
                                   p("Boxplot (multiple replicates) or dotplot (single replicate) of gene expression (logCPM)."),
                                   br(),
                                   column(4, plotOutput("plot_treatment1")),
                                   column(4, plotOutput("plot_treatment2")),
                                   column(4, plotOutput("plot_treatment3")),
                                   column(4, plotOutput("plot_treatment4")),
                                   column(4, plotOutput("plot_treatment5"))
                          ),
                          tabPanel("PCOS / other human",
                                   p("See info tab for details of these datasets"),
                                   br(),
                                   column(4, plotOutput("plot_PCOS1_box")),
                                   column(4, plotOutput("plot_PCOS2_box")),
                                   column(4, plotOutput("plot_PCOS3_box")),
                                   column(4, plotOutput("plot_PCOS4_box")),
                                   column(4, plotOutput("plot_PCOS6_box")),
                                   column(4, plotOutput("plot_PCOS7_box")),
                                   p("Boxplot (multiple replicates) of gene expression (RMA count for microarray, logCPM or FPKM for RNA-seq)."),
                                   column(6, plotOutput("plot_PCOS8_pca")),
                                   column(6, plotOutput("plot_PCOS9_pca")),
                                   p("MDS plot with gene expression level (logCPM) for un-labelled human GC RNA-seq samples.")
                          ),
                          tabPanel("Datasets",
                                   h2(strong("The following datasets are included in this app:"), style = "font-size:15px;"),
                                   br(),
                                   DT::dataTableOutput("table1"),
                                   br(),
                                   h2(strong("The following datasets are not included in this app due to technical issues but might also be of interest:"), style = "font-size:15px;"),
                                   DT::dataTableOutput("table2"),
                                   br(),
                                   HTML("<p>In-house ChIP/ATAC/other genomic data can also be viewed at the UCSC Genome Browser for <a href='http://genome.ucsc.edu/s/thao.dinh/complete_GC_GR_AR'>mouse</a> and <a href='http://genome.ucsc.edu/s/thao.dinh/human_ATAC_PGR_GSE86189'>human</a>.</p>")
                          ),
                          tabPanel("Analysis details",
                                   h2(strong("In-house RNA-seq datasets:"), style = "font-size:25px;"),
                                   p("- Details on in-house mouse RNA-seq datasets can be found at https://doi.org/10.1101/2021.06.17.448908[to be updated]."),
                                   p("- In-house human GC RNA-seq datasets:"),
                                   p("Human granulosa cells were purified from follicular fluid of patients seeking ART, with red blood cells removed using RBC lysis buffer. Cells were counted and 200k live cells were used for RNA extraction and total RNA library preparation."),
                                   p("Libraries were paired-end sequenced on the NovaSeq platform at ~38M reads/library. Bioinformatics analysis includes adaptor trimming, transcriptomic alignment using Salmon and differential analysis following the limma/voom pipeline. LogCPM is used as normalised count. Various clinical and fertility outcome data is incorporated with each library. No libraries were exluded from this analysis."),
                                   br(),
                                   p("- In-house KGN RNA-seq datasets:"),
                                   p("KGN were overexpressed with dox-inducible PGR-A or PGR-B and treated with R5020 for indicated time. Cells were harvested for RNA extraction and mRNA library prepared. Libraries were sequenced on the NovaSeq platform at ~30M reads/library. Bioinformatics analysis includes adaptor trimming, transcriptomic alignment using Salmon and differential analysis following the Deseq2 pipeline. LogTPM is used as normalised count. No libraries were exluded from this analysis."),
                                   br(),
                                   h3(strong("Public RNA-seq datasets:"), style = "font-size:25px;"),
                                   p("Raw read count data was obtained from GEO without re-alignment. Data was filtered, normalised and logCPM was calculated using edgeR/limma packages. When raw data is unavailable, normalised count (FPKM) was used in place of logCPM. No libraries were exluded from any datasets."),
                                   p(""),
                                   br(),
                                   h4(strong("Public microarray datasets:"), style = "font-size:25px;"),
                                   p("CEL files were obtained from microarray datasets and normalised using the RMA method from oligo package. No samples were excluded from any datasets.")
                          )
    )
    )
  )
)

# Server side
server <- function(input, output, session) 
{updateSelectizeInput(session, 'pickme', choices = list, server = TRUE)
  
  # this next chunk is to show a text and not a plot when the specified gene is not in the human gene list or mouse gene list
  human_text<-function(i){
    output[[paste0("text_human",i)]] <- renderText({
      if (!(input$pickme %in% lcpm$gene)) {
        paste("This gene either has no human ortholog or has been filtered out due to low expression. If this is a mouse gene, you can check the other tabs for data in mouse.")}
    })}
  lapply(1:3, human_text)
  
  mouse_text<-function(i){
    output[[paste0("text_mouse",i)]] <- renderText({
      if (!(input$pickme %in% mouse$adj_name)) {
        paste("This gene likely has no mouse ortholog. If this is a human gene, you can check the other tabs for data in human.")}
    })}
  lapply(1:3, mouse_text)
  
  output$text_human_timecourse <- renderText({
    if (!(input$pickme %in% GSE133868$gene)) {
      paste("This gene either has no human ortholog or has been filtered out due to low expression. If this is a mouse gene, you can check the other tabs for data in mouse.")
    }
  })
  
  # this next bit is to create plots and accompanied text
  output$plot_inhouse1 <- renderPlot({
    if (input$pickme %in% lcpm$gene) {
      filtered <- filter(lcpm, gene == input$pickme)
      mds_t <- as.data.frame(t(mds))
      colnames(mds_t) <- mds_t[4, ]
      mds_t <- mds_t[-4 ,]
      filtered <- subset(filtered, select = -c(gene) )
      merge <- rbind(mds_t,filtered)
      merge_t <- as.data.frame(t(merge))
      colnames(merge_t)[21] <- "gene"
      merge_t$file_ID <- row.names(merge_t)
      cols <- colnames(merge_t)
      cols <- cols[! cols %in% c('group', 'lib_batch', 'file_ID')]
      merge_t[cols] <- sapply(merge_t[cols],as.numeric)
      merge_t %>% ggplot(aes(Dim1, Dim2, color = gene, shape = group)) +
        geom_point(size=3)+
        scale_colour_gradient(low = "black", high = "red")}
  })
  output$plot_inhouse2 <- renderPlot({
    if (input$pickme %in% lcpm$gene) {
      g1 <- subset(atac_vs_cpm, gene == input$pickme)
      ggplot(data = atac_vs_cpm, aes(x = logCPM, y = logATAC)) + geom_point(col = "darkgrey") + 
        geom_point(data = g1, col = "red", size=5) + 
        geom_text(data = g1, label = g1$gene, vjust = 1, size=5) + 
        labs(title = "logATAC vs logCPM", x = "logATAC", y = "logCPM") + 
        theme_minimal() 
    }
  })
  output$plot_inhouse3 <- renderPlot({
    if (input$pickme %in% lcpm$gene) {
      filtered <- filter(lcpm, gene == input$pickme)
      mds_t <- as.data.frame(t(mds))
      colnames(mds_t) <- mds_t[4, ]
      mds_t <- mds_t[-4 ,]
      filtered <- subset(filtered, select = -c(gene) )
      merge <- rbind(mds_t,filtered)
      merge_t <- as.data.frame(t(merge))
      colnames(merge_t)[21] <- "gene"
      merge_t$file_ID <- row.names(merge_t)
      cols <- colnames(merge_t)
      cols <- cols[! cols %in% c('group', 'lib_batch', 'file_ID')]
      merge_t[cols] <- sapply(merge_t[cols],as.numeric)
      merge_t <- merge_t[order(merge_t$gene),]
      plot = HeatmapAnnotation(logCPM = anno_barplot(merge_t$gene, height = unit(10, "cm")),
                               text = anno_text(merge_t$file_ID, rot = 45, gp = gpar(fontsize = 10)),
                               annotation_name_side = "left",
                               annotation_name_gp= gpar(fontsize = 8)) %v% 
        HeatmapAnnotation(group = anno_simple(merge_t$group, col = c("ML" = "black", "GOLD" = "gold")),
                          batch = anno_simple(merge_t$lib_batch, col = c("first" = "#1B1464", "second" = "#65C8D0", "trial" = "#FF9289")),
                          BMI = anno_simple(merge_t$BMI, col = colorRamp2(c(17, 51), c("#ABC3DD", "#8856A7"))),
                          age = anno_simple(merge_t$age, col = colorRamp2(c(19, 47), c("#ffffe0", "#94003a"))),
                          repeated = anno_simple(merge_t$rep_patient_ID, col = c("0" = "#240e8b", "1" = "#f04393", "2" = "#f9c449", "3" = "#ce6c47", "4" = "#e8a49c")),
                          eggs = anno_simple(merge_t$N_EGGS, col = colorRamp2(c(2, 33), c("#EDF8FB", "#3BAA71"))),
                          first_stim = anno_simple(merge_t$STIM_FIRST, col = c("0" = "#fbf6b3", "1" = "#DD1C77")),
                          tubal = anno_simple(merge_t$CI_TUBE, col = c("0" = "#240e8b", "1" = "#f9c449")),
                          endometriosis = anno_simple(merge_t$CI_ENDO, col = c("0" = "#240e8b", "1" = "#f9c449")),
                          PCOS = anno_simple(merge_t$PCOS, col = c("0" = "#240e8b", "1" = "#f9c449")),
                          other_female = anno_simple(merge_t$CI_OTH, col = c("0" = "#240e8b", "1" = "#f9c449")),
                          male = anno_simple(merge_t$CI_MALE, col = c("0" = "#240e8b", "1" = "#f9c449")),
                          unexplained = anno_simple(merge_t$CI_UNEX, col = c("0" = "#240e8b", "1" = "#f9c449")),
                          fertilisation_rate = anno_simple(merge_t$FERT_RATE, colorRamp2(c(0, 1), c("#fbf6b3", "#D95F0E"))),
                          ontime_day3_rate = anno_simple(merge_t$DAY3_RATE, colorRamp2(c(0, 1), c("#fbf6b3", "#D95F0E"))),
                          blast_rate = anno_simple(merge_t$BLAST_RATE, colorRamp2(c(0, 1), c("#fbf6b3", "#D95F0E"))),
                          preg_rate = anno_simple(merge_t$CUM_PR_RATE, colorRamp2(c(0, 1), c("#fbf6b3", "#D95F0E"))),
                          foetal_heart = anno_simple(merge_t$CUM_FH_RATE, colorRamp2(c(0, 1), c("#fbf6b3", "#D95F0E"))))
      lgd1 = Legend(labels = c("ML", "GOLD"), title = "group", legend_gp = gpar(fill = c("black", "gold")))
      lgd2 = Legend(labels = c("first", "second", "trial"), title = "batch", legend_gp = gpar(fill = c("#1B1464", "#65C8D0", "#FF9289")))
      lgd3 = Legend(col_fun = colorRamp2(c(17, 51), c("#ABC3DD", "#8856A7")), title = "BMI", at = c(17,51))
      lgd4 = Legend(col_fun = colorRamp2(c(17, 51), c("#ffffe0", "#94003a")), title = "age", at = c(19,47))
      lgd5 = Legend(labels = c("non", "repeated#1", "repeated#2","repeated#3","repeated#4"), title = "repeated", legend_gp = gpar(fill = c("#240e8b", "#f04393", "#f9c449", "#ce6c47", "#e8a49c")))
      lgd6 = Legend(col_fun = colorRamp2(c(2, 33), c("#EDF8FB", "#3BAA71")), title = "eggs", at = c(2,33))
      lgd7 = Legend(labels = c("non-first", "first"), title = "first stim", legend_gp = gpar(fill = c("#FFFFCC", "#DD1C77")))
      lgd8 = Legend(labels = c("none", "diagnosed"), title = "infertility cause", legend_gp = gpar(fill = c("#240e8b", "#f9c449")))
      lgd9 = Legend(col_fun = colorRamp2(c(0,1), c("#fbf6b3", "#D95F0E")), title = "fertility outcome rate", at = c(0,1))
      pd = packLegend(lgd1, lgd2, lgd3, lgd4, lgd5, lgd6, lgd7, lgd8, lgd9, direction = "horizontal")
      draw(plot, annotation_legend_list = pd, annotation_legend_side = "bottom")
    }
  })
  output$plot_inhouse4 <- renderPlot({
    if (input$pickme %in% mouse$adj_name) {
      g1 <- subset(mouse, adj_name == input$pickme)
      ggplot(data = mouse, aes(x = LH_FC, y = LH_padj)) + geom_point(col = "darkgrey") + 
        geom_point(data = g1, col = "red", size=5) + 
        geom_text(data = g1, label = g1$adj_name, vjust = 1, size=5) + 
        geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") + 
        geom_hline(yintercept = 2, col = "black", linetype = "dashed") + 
        labs(x = "logFC", y = "-log(FDR)") + theme_minimal()+
        ggtitle("LH RNA-seq - 8h vs 0h")
    }
  })
  output$plot_inhouse5 <- renderPlot({
    if (input$pickme %in% mouse$adj_name) {
      g1 <- subset(mouse, adj_name == input$pickme)
      ggplot(data = mouse, aes(x = PGRKO_FC, y = PGRKO_padj)) + geom_point(col = "darkgrey") + 
        geom_point(data = g1, col = "red", size=5) + 
        geom_text(data = g1, label = g1$adj_name, vjust = 1, size=5) + 
        geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") + 
        geom_hline(yintercept = 2, col = "black", linetype = "dashed") + 
        labs(x = "logFC", y = "-log(FDR)") + theme_minimal()+
        ggtitle("PGRKO vs PGRWT 8h hCG")
    }
  })
  
  output$plot_KGN1 <- renderPlot({
    if (input$pickme %in% kgn_count$gene) {
      new <- filter(kgn_count, gene == input$pickme)
      new %>% ggplot(aes(x = factor(second_column,levels =c('PRA_0h','PRA_8h','PRA_24h','PRB_0h','PRB_8h','PRB_24h')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "group", y = "logTPM", fill="gene/transcript") + 
        theme_bw() +
        theme(legend.position = "bottom", legend.direction = "horizontal",
              text = element_text(size = 15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
  })
  output$plot_KGN2 <- renderPlot({
    if (input$pickme %in% kgn_volcano$gene) {
      g1 <- subset(kgn_volcano, gene == input$pickme)
      ggplot(data = kgn_volcano, aes(x = PRA_8hvs0h_logFC, y = PRA_8hvs0h_logpvaj)) + geom_point(col = "darkgrey") + 
        geom_point(data = g1, col = "red", size=5) + 
        geom_text(data = g1, label = g1$gene, vjust = 1, size=5) + 
        geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") + 
        geom_hline(yintercept = 2, col = "black", linetype = "dashed") + 
        labs(x = "logFC", y = "-log(FDR)") + theme_minimal() +
        ggtitle("PRA 8h vs 0h")
    }
  })
  output$plot_KGN3 <- renderPlot({
    if (input$pickme %in% kgn_volcano$gene) {
      g1 <- subset(kgn_volcano, gene == input$pickme)
      ggplot(data = kgn_volcano, aes(x = PRA_24hvs0h_logFC, y = PRA_24hvs0h_logpvaj)) + geom_point(col = "darkgrey") + 
        geom_point(data = g1, col = "red", size=5) + 
        geom_text(data = g1, label = g1$gene, vjust = 1, size=5) + 
        geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") + 
        geom_hline(yintercept = 2, col = "black", linetype = "dashed") + 
        labs(x = "logFC", y = "-log(FDR)") + theme_minimal() +
        ggtitle("PRA 24h vs 0h")
    }
  })
  output$plot_KGN4 <- renderPlot({
    if (input$pickme %in% kgn_volcano$gene) {
      g1 <- subset(kgn_volcano, gene == input$pickme)
      ggplot(data = kgn_volcano, aes(x = PRB_8hvs0h_logFC, y = PRB_8hvs0h_logpvaj)) + geom_point(col = "darkgrey") + 
        geom_point(data = g1, col = "red", size=5) + 
        geom_text(data = g1, label = g1$gene, vjust = 1, size=5) + 
        geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") + 
        geom_hline(yintercept = 2, col = "black", linetype = "dashed") + 
        labs(x = "logFC", y = "-log(FDR)") + theme_minimal() +
        ggtitle("PRB 8h vs 0h")
    }
  })
  output$plot_KGN5 <- renderPlot({
    if (input$pickme %in% kgn_volcano$gene) {
      g1 <- subset(kgn_volcano, gene == input$pickme)
      ggplot(data = kgn_volcano, aes(x = PRB_24hvs0h_logFC, y = PRB_24hvs0h_logpvaj)) + geom_point(col = "darkgrey") + 
        geom_point(data = g1, col = "red", size=5) + 
        geom_text(data = g1, label = g1$gene, vjust = 1, size=5) + 
        geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") + 
        geom_hline(yintercept = 2, col = "black", linetype = "dashed") + 
        labs(x = "logFC", y = "-log(FDR)") + theme_minimal() +
        ggtitle("PRB 24h vs 0h")
    }
  })
  
  output$plot_timecourse1 <- renderPlot({
    if (input$pickme %in% GSE140371$original) {
      new <- filter(GSE140371, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = second_column, y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "time post PMSG", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("GSE140371") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_timecourse2 <- renderPlot({
    if (input$pickme %in% GSE15315$original) {
      new <- filter(GSE15315, original == input$pickme)
      new <- subset(new, third_column == "WT") #for geno/time matrices
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = second_column, y = value)) + 
        geom_dotplot(aes(fill = gene),  binaxis='y', stackdir='center', 
                     position = position_jitter(width = 0.025, height = 0)) +
        labs(x = "time post hCG", y = "RMA count") + 
        theme_bw() +
        ggtitle("GSE15315") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_timecourse3 <- renderPlot({
    if (input$pickme %in% GSE44651$original) {
      new <- filter(GSE44651, original == input$pickme)
      new <- subset(new, third_column == "WT") #for geno/time matrices
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = second_column, y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "time post hCG", y = "RMA count", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("GSE44651") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_timecourse4 <- renderPlot({
    if (input$pickme %in% GSE119508$original) {
      new <- filter(GSE119508, original == input$pickme)
      new <- subset(new, third_column == "WT") #for geno/time matrices
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = second_column, y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "time post hCG", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("GSE119508") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_timecourse5 <- renderPlot({
    if (input$pickme %in% GSE167939$original) {
      new <- filter(GSE167939, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(second_column, level = c('0h hCG', '4h hCG', '12h hCG')), y = value)) + 
        geom_dotplot(aes(fill = gene),  binaxis='y', stackdir='center', 
                     position = position_jitter(width = 0.025, height = 0)) +
        labs(x = "time post hCG", y = "FPKM") + 
        theme_bw() +
        ggtitle("GSE167939") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_timecourse6 <- renderPlot({
    if (input$pickme %in% GSE20466$original) {
      new <- filter(GSE20466, original == input$pickme)
      new <- subset(new, third_column == "WT") #for geno/time matrices
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = second_column, y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "time post hCG", y = "RMA count", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("GSE20466") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_timecourse7 <- renderPlot({
    if (input$pickme %in% GSE178314$original) {
      new <- filter(GSE178314, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = second_column, y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "time post hCG", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("GSE178314") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_timecourse8 <- renderPlot({
    if (input$pickme %in% GSE23084$original) {
      new <- filter(GSE23084, original == input$pickme)
      new <- subset(new, third_column == "WT") 
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(second_column, level = c('8h hCG', '24h hCG')), y = value)) + 
        geom_dotplot(aes(fill = gene),  binaxis='y', stackdir='center', 
                     position = position_jitter(width = 0.025, height = 0)) +
        labs(x = "time post hCG", y = "RMA count") + 
        theme_bw() +
        ggtitle("GSE23084") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_timecourse9 <- renderPlot({
    if (input$pickme %in% GSE133868$gene) {
      new <- filter(GSE133868, gene == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = timepoint, y = count)) + 
        geom_dotplot(aes(fill = gene),  binaxis='y', stackdir='center', 
                     position = position_jitter(width = 0.025, height = 0)) +
        labs(x = "time post hCG", y = "RMA count") + 
        theme_bw() +
        ggtitle("GSE133868") +
        theme(text = element_text(size = 15))
    }
  })
  
  output$plot_KO1 <- renderPlot({
    if (input$pickme %in% GSE168213$original) {
      new <- filter(GSE168213, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(third_column,levels =c('WT','KO')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "genotype", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("PRWT vs KO, 8h post-hCG") +
        theme(legend.position = "bottom", legend.direction = "horizontal",
              text = element_text(size = 15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
  })
  output$plot_KO2 <- renderPlot({
    if (input$pickme %in% GSE140922$original) {
      new <- filter(GSE140922, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(third_column,levels =c('WT','KO')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "genotype", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("CBFB/RUNX2 WT vs KO, 11h post-hCG") +
        theme(legend.position = "bottom", legend.direction = "horizontal",
              text = element_text(size = 15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
  })
  output$plot_KO3 <- renderPlot({
    if (input$pickme %in% GSE119508$original) {
      new <- filter(GSE119508, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(interaction(third_column, second_column, sep=':'), levels=c('WT:0h hCG', 'KO:0h hCG', 'WT:4h hCG', 'KO:4h hCG')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "genotype:timepoint", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("LRH1 WT vs KO, 0h/4h post-hCG") +
        theme(legend.position = "bottom", legend.direction = "horizontal",
              text = element_text(size = 15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
  })
  output$plot_KO4 <- renderPlot({
    if (input$pickme %in% GSE23084$original) {
      new <- filter(GSE23084, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(interaction(third_column, second_column, sep=':'), levels=c('WT:8h hCG', 'KO:8h hCG', 'WT:24h hCG', 'KO:24h hCG')), y = value)) + 
        geom_dotplot(aes(fill = gene),  binaxis='y', stackdir='center', 
                     position = position_jitter(width = 0.025, height = 0)) +
        labs(x = "time post hCG", y = "RMA count") + 
        theme_bw() +
        ggtitle("CEBPB WT vs KO, 8h/24h post-hCG") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
  })
  output$plot_KO5 <- renderPlot({
    if (input$pickme %in% GSE44651$original) {
      new <- filter(GSE44651, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(interaction(third_column, second_column, sep=':'), levels=c('WT:0h hCG', 'KO:0h hCG', 'WT:4h hCG', 'KO:4h hCG')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "genotype:timepoint", y = "RMA count", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("ERb WT vs KO, 0h/4h post-hCG") +
        theme(legend.position = "bottom", legend.direction = "horizontal",
              text = element_text(size = 15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
  })
  output$plot_KO6 <- renderPlot({
    if (input$pickme %in% GSE15315$original) {
      new <- filter(GSE15315, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(interaction(third_column, second_column, sep=':'), levels=c('WT:0h hCG', 'KO:0h hCG', 'WT:2.5h hCG', 'KO:2.5h hCG','WT:4h hCG', 'KO:4h hCG')), y = value)) + 
        geom_dotplot(aes(fill = gene),  binaxis='y', stackdir='center', position = position_jitter(width = 0.025, height = 0)) +
        labs(x = "genotype:timepoint", y = "logCPM") + 
        theme_bw() +
        ggtitle("ERK1/2 WT vs KO, 0h/2.5h/4h post-hCG") +
        theme(legend.position = "bottom", legend.direction = "horizontal",
              text = element_text(size = 15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
  })
  output$plot_KO7 <- renderPlot({
    if (input$pickme %in% GSE20466$original) {
      new <- filter(GSE20466, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(interaction(third_column, second_column, sep=':'), levels=c('WT:0h hCG', 'KO:0h hCG', 'WT:6h hCG', 'KO:6h hCG')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "genotype:timepoint", y = "RMA count", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("INHA WT vs KO, 0h/6h post-hCG") +
        theme(legend.position = "bottom", legend.direction = "horizontal",
              text = element_text(size = 15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
  })
  
  output$plot_treatment1 <- renderPlot({
    if (input$pickme %in% GSE152727$original) {
      new <- filter(GSE152727, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(second_column, levels=c('VEH','DHT')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "treatment", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("Mouse GC, 24h DHT or vehicle") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_treatment2 <- renderPlot({
    if (input$pickme %in% GSE158218_mouse$original) {
      new <- filter(GSE158218_mouse, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(second_column, levels=c('VEH','DHT')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "treatment", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("Mouse in vivo, 18h DHT or vehicle") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_treatment3 <- renderPlot({
    if (input$pickme %in% GSE158218_human$original) {
      new <- filter(GSE158218_human, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(second_column, levels=c('VEH','DHT')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "treatment", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("KGN, 12h DHT or vehicle") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_treatment4 <- renderPlot({
    if (input$pickme %in% GSE166443$original) {
      new <- filter(GSE166443, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(second_column, levels=c('hCG','hCG+T5224')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "treatment", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("human GC, 12h hCG+FOS inhibitor or hCG") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_treatment5 <- renderPlot({
    if (input$pickme %in% GSE137608$original) {
      new <- filter(GSE137608, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(second_column, levels=c('control','FSH', 'FSK', 'PMA', 'SC79')), y = value)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center',
                     position = position_jitter(width = 0.025, height = 0)) +
        labs(x = "treatment", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("KGN, 24h treatment of protein kinase activators") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  
  output$plot_PCOS1_box <- renderPlot({
    if (input$pickme %in% GSE114419$original) {
      new <- filter(GSE114419, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(second_column, levels=c('control','PCOS')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "group", y = "RMA count", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("GSE114419") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_PCOS2_box <- renderPlot({
    if (input$pickme %in% GSE138518$gene) {
      new <- filter(GSE138518, gene == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(second_column, levels=c('control','PCOS')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "group", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("GSE138518") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_PCOS3_box <- renderPlot({
    if (input$pickme %in% GSE168404$gene) {
      new <- filter(GSE168404, gene == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(second_column, levels=c('control','PCOS')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "group", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("GSE168404") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_PCOS4_box <- renderPlot({
    if (input$pickme %in% GSE193123$gene) {
      new <- filter(GSE193123, gene == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = factor(second_column, levels=c('control','PCOS')), y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "group", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("GSE193123") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  output$plot_PCOS6_box <- renderPlot({
    if (input$pickme %in% GSE80432$original) {
      new <- filter(GSE80432, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = second_column, y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "group", y = "RMA count", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("GSE80432") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
  })
  output$plot_PCOS7_box <- renderPlot({
    if (input$pickme %in% GSE168214$original) {
      new <- filter(GSE168214, original == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(x = second_column, y = value)) + 
        geom_boxplot(aes(fill = gene), width = 0.1, size = 0.4, position = position_nudge(x = -0.1)) +
        geom_dotplot(aes(fill = gene), binaxis='y', stackdir='center') +
        labs(x = "group", y = "logCPM", fill="gene/transcript") + 
        theme_bw() +
        ggtitle("GSE168214") +
        theme(legend.position = "bottom", 
              legend.direction = "horizontal",
              text = element_text(size = 15))
    }
  })
  
  output$plot_PCOS8_pca <- renderPlot({
    if (input$pickme %in% GSE110924$gene) {
      new <- filter(GSE110924, gene == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(Dim1, Dim2, color = value)) +
        geom_point(size=3)+
        scale_colour_gradient(low = "black", high = "red") +
        labs(color="logCPM") +
        theme_bw() + 
        ggtitle("GSE110924")
    }
  })
  output$plot_PCOS9_pca <- renderPlot({
    if (input$pickme %in% GSE62093$gene) {
      new <- filter(GSE62093, gene == input$pickme)
      new <- as.data.frame(new) 
      new %>% ggplot(aes(Dim1, Dim2, color = value)) +
        geom_point(size=3)+
        scale_colour_gradient(low = "black", high = "red") +
        labs(color="logCPM") +
        theme_bw() + 
        ggtitle("GSE62093")
    }
  })
  
  output$table1 = DT::renderDataTable({table1})
  output$table2 = DT::renderDataTable({table2})
  
  session$onSessionEnded(function() { #this is to stop R from crashing when closing app window rather than Stop
    stopApp()
  })
}


shinyApp(ui = ui, server = server)

