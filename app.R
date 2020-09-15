## app.R ##
#options(repos = BiocInstaller::biocinstallRepos())
#getOption("repos")
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(rsconnect))
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(shinythemes))
suppressPackageStartupMessages(library(rintrojs))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(dyno))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(enrichplot))

ui <- dashboardPage(
  dashboardHeader(
    title = "Protein and transcriptional responses of blood derived human immune cells to diverse stimuli at single cell resolution",
    titleWidth = 1200
  ),
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(id = "tabs",
      # dashboard should display tutorial or instructions of how to use the app
      menuItem("Overview", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Cell Proportions", tabName = "Flow", icon = icon("chart-bar")),
      menuItem("Gene Expression", tabName = "Single", icon = icon("chart-line")),
      menuItem("ADT Expression", tabName = "ADT", icon = icon("chart-line")),
      menuItem("Response Pathways", tabName = "Response", icon = icon("reply")),
      menuItem("Trajectories", tabName = "Trajectory", icon = icon("route")),
      menuItem("Source Code/Questions?", icon = icon("file-code-o"), 
               href = "https://github.com/nlawlor/PBMC_CITEseq")
    )
  ),
  
  ## Body content
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard", 
              
              h2("The data presented here is associated with the following publication/preprint:"),
              fluidRow(
                # image of journal and link to the preprint/etc
                box()
              ),
              
              h2("This study was funded and supported by:"),
              fluidRow(
                box(width = 2, tags$a("",
                            href="https://chanzuckerberg.com", target="_blank",
                            tags$img(src = "czi.png", height = "125", width = "200")
                ))),
              h2("This study is a collaboration of the following groups:"),
              fluidRow(
                box(width = 2, tags$a("", 
                            href="https://www.humancellatlas.org", target="_blank",
                            tags$img(src = "human.cell.atlas.png", height = "150", width = "200")
                )),
                box(width = 2, tags$a("", 
                           href="https://www.jax.org", target="_blank",
                           tags$img(src = "jax.png", height = "125", width = "200")
                )),
                box(width = 2, tags$a("", 
                           href="https://www.nygenome.org", target="_blank",
                           tags$img(src = "nygc.png", height = "150", width = "200")
                ))
              )
      ),
      
      # Second tab content
      tabItem(tabName = "Flow",
              h2("Cell Type Proportions from Flow Cytometry and Single Cell Data"),
              # make ggplot interactive points showing the # cell type proportions 
              fluidRow(
                box(plotlyOutput("flow_plot1", height = 500)),
                box(plotlyOutput("flow_plot2", height = 500)),
                box(plotlyOutput("flow_plot", height = 500))
              )
              
      ),
      
      tabItem(tabName = "Single",
              h2("Single Cell RNA-seq Gene Expression"),
              # plot of gene expression
              fluidRow(
                box(
                  title = "Select Gene", 
                  selectInput(inputId = "gen_sym", label = "Gene Symbol: ",
                              choices = NULL, 
                              selected = NULL,
                              multiple = FALSE),
                width = 2),
                box(
                  title = "Plot By",
                  selectInput(inputId = "gen_plot_type", label = "Facet: ",
                              choices = c("Treatment", "Cell Type"), 
                              selected = "Treatment",
                              multiple = FALSE),
                  width = 2),
                box(plotOutput("rna_plot", height = 700), width = 8)
              )
      ),
      # end of gene expression tab
      
      tabItem(tabName = "ADT",
              h2("Single Cell Antibody Derived Tag (ADT) Expression"),
              fluidRow(
                box(title = "Select ADT", 
                    selectInput(inputId = "adt_sym", label = "ADT Symbol: ",
                                choices = NULL, 
                                selected = NULL,
                                multiple = FALSE),
                    width = 2),
                box(title = "Plot By",
                    selectInput(inputId = "adt_plot_type", label = "Facet: ",
                              choices = c("Treatment", "Cell Type"), 
                              selected = "Treatment",
                              multiple = FALSE),
                    width = 2),
                box(plotOutput("adt_plot", height = 700), width = 8)
              )
      ),
      # end of ADT expression tab
      
      tabItem(tabName = "Response",
              h2("Cell Type Response Genes and Pathways associated with Stimulation"),
              # plotting space
              fluidRow(
                column(width = 12,
                  # specify response and cell type
                  box(
                    title = "Select Response and Cell Type", 
                    selectInput(inputId = "cell_type", label = "Choose: ",
                                choices = c("Bcell_Anti_CD3_CD28",
                                            "Monocyte_LPS_induced", "Monocyte_LPS_reduced",
                                            "NK_Anti_CD3_CD28",
                                            "CD4T_naive_Anti_CD3_CD28", "CD4T_memory_Anti_CD3_CD28",
                                            "CD8T_naive_Anti_CD3_CD28", "CD8T_memory_Anti_CD3_CD28"), 
                                selected = NULL,
                                multiple = FALSE),
                    width = 5),
                  box(
                    title = "Select Pathway Type", 
                    selectInput(inputId = "path_type", label = "Choose: ",
                                choices = c("BP", "MF", "CC", "KEGG", "Wikipathways", "DO", "DGN", "Tmod"), 
                                selected = NULL,
                                multiple = FALSE),
                    width = 2),
                  box(
                    title = "Number of Terms to Display",
                    numericInput(inputId = "term_num", label = "# of Terms",
                                 value = 10, min = 0, max = NA, step = 5),
                    width = 3
                  ),
                  box(
                    title = "Select Plot Type", 
                    selectInput(inputId = "view_type", label = "Choose: ",
                                choices = c("Dotplot", "Barplot", "Map"), 
                                selected = NULL,
                                multiple = FALSE),
                    width = 2)
                ),
                
                box(plotOutput("path_plot", height = 700), width = 8),
                
                # table of genes 
                box(
                  h3("Table of Cell Type Response Genes"),
                  downloadButton("Download_Pathways", "Download Table"),
                  DT::dataTableOutput("path_table"), width = 12
                )
              )
      ), # end of response tab
      
      tabItem(tabName = "Trajectory",
              h2("Pseudo-temporal ordering and trajectory inference"),
              # plotting space make 2 side by side plots of trajectories (to color points by metadata, and color by gene expression)
              fluidRow(
                column(width = 12,
                       # specify response and cell type
                       box(
                         title = "Select Trajectory", 
                         selectInput(inputId = "traj_name", label = "Choose: ",
                                     choices = c("Bcell_Anti_CD3",
                                                 "CD4T_Anti_CD3", "CD8T_Anti_CD3",
                                                 "NK_Anti_CD3",
                                                 "Monocyte_LPS"), 
                                     selected = "Monocyte_LPS",
                                     multiple = FALSE),
                         width = 4),
                       box(
                         title = "Color Trajectory By",
                         selectInput(inputId = "traj_plot_type", label = "Metadata: ",
                                     choices = c("Stimulation", "Pseudotime"), 
                                     selected = NULL,
                                     multiple = FALSE),
                         width = 2),
                       box(
                         title = "Select Gene/ADT to View on Trajectory", 
                         selectInput(inputId = "traj_marker_id", label = "Choose: ",
                                     choices = NULL, 
                                     selected = NULL,
                                     multiple = FALSE),
                         width = 4)
                ) # end of first column
              ), # end of first fluid row
              
              # second row will have plots
              fluidRow(
                column(width = 12,
                # first traj plot
                box(plotOutput("traj_plot", height = 600), width = 6),
                # second traj plot
                box(plotOutput("traj_exp_plot", height = 600), width = 6)
              )
            ) # end of second fluid row
              
      ) # end of trajectory tab
    )
  )
)

# code and functions to actually do stuff
server <- shinyServer(function(input, output, session) {
  #set.seed(122)
  
  # table of event reactive values (for loading files just once)
  dataTables <- reactiveValues(
    flow_cytom_prop = NULL,
    rna_cell_prop = NULL,
    flow_and_rna_prop = NULL,
    adt_exp_data = NULL,
    gene_exp_data = NULL,
    gene_exp_sel = NULL,
    gene_exp_names = NULL,
    gene_letter_choice = NULL,
    dat_obj = NULL,
    model_obj = NULL
  )
  
  # tab for cell proportion info
  observeEvent(input$tabs,{
    if (input$tabs == "Flow") {
      # barplot of flow cytom proportions
      if (is.null(dataTables$flow_cytom_prop)) {
        output$flow_plot1 <- renderPlotly({
          # read in data
          withProgress(expr = flow_df <- readRDS("Data/Flow.cytometry.cell.type.proportions.per.donor.Rds"),
                                  message = "Loading flow cytometry cell proportions, please wait")
          f1 <- ggplot(flow_df, aes(fill=CellType, y=Cell_Type_Proportion, x=Donor, label = Sex)) + 
            geom_bar(stat="identity", position = "fill") +
            labs(y = "Cell % Per Donor", x = "Donor Number") + scale_fill_manual(values = c("#fdbf6f", "red2", "#2171b5", "#c2a5cf")) +
            ggtitle("Flow Cytometry Cell Type Proportions")
          ggplotly(f1, tooltip = c("y", "x", "Sex", "Donor", "CellType"))
        })
        # add reactive value
        dataTables$flow_cytom_prop <- 1
      }
      
      # barplot of scrna-seq cell type proportions
      # order same as flow cytom proportions
      if (is.null(dataTables$rna_cell_prop)) {
        output$flow_plot2 <- renderPlotly({
          # read in data
          withProgress(expr = cell_data <- readRDS("Data/scRNAseq.cell.proportions.per.donor.Rds"),
                                    message = "Loading scRNAseq cell proportions, please wait")
          # color pallete
          cols <- c( "#fdbf6f", "#f7fcb9", "#41ab5d", "lightskyblue1", "#2171b5", "red2", "#c2a5cf")
          # celltype info
          celltype <- c("B", "CD4T_Naive", "CD4T_Mem", "CD8T_Naive", "CD8T_Mem",
                        "CD14_Mono", "NK")
          
          f2 <- ggplot(cell_data, aes(fill=CellType, y=Percent, x=Donor, label = Cell_Number)) + 
            geom_bar(stat="identity", position="fill") +
            labs(y = "Cell % Per Donor", x = "Donor") + scale_fill_manual(values = cols, labels = celltype) +
            ggtitle("scRNA-seq Cell Type Proportions") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
          ggplotly(f2, tooltip = c("y", "x", "Cell_Number", "Donor", "CellType"))
          })
        
        # add reactive value
        dataTables$rna_cell_prop <- 1
      }
      
      # plot for flow cytom vs single cell rna-seq proportions
      if (is.null(dataTables$flow_and_rna_prop)) {
        output$flow_plot <- renderPlotly({
          
          # read in data
          withProgress(expr = flow_dat <- readRDS("Data/scRNA_flow_baseline_proportions.Rds"),
                                   message = "Loading scRNA and flow cytometry cell proportions, please wait")
          p1 <- ggplot(flow_dat, aes(y=Control, x=Single_Cell_Prop, group = CellType, color = CellType)) + 
            geom_point(aes(color=CellType, pch = Sex, labels = Donor_Number)) + 
            scale_y_continuous(limits = c(0,100)) +
            scale_x_continuous(limits = c(0,100)) + 
            geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
            labs(y = "Flow Cytometry Cell %", x = "scRNA-seq Cell %") +
            scale_color_manual(values = c("#fdbf6f", "red2", "#2171b5", "#c2a5cf")) +
            ggtitle("Cell Type Proportions (Baseline)")
          ggplotly(p1, tooltip = c("y", "x", "group", "Sex", "Donor_Number"))
          })
        # add reactive value
        dataTables$flow_and_rna_prop <- 1
      }
      
    }
  })
  
  # tab for gene expression data
  observeEvent(input$tabs,{
    if (input$tabs == "Single") {
      
      # load in gene expression names
      if (is.null(dataTables$gene_exp_names)) {
        withProgress(message = "Loading gene symbols, please wait",
                                   expr = dataTables$gene_exp_names <- readRDS("Data/gene.expression.names.no.IgG.IgM.Rds"))
        
        # update gene choices
        updateSelectInput(session = session, inputId = "gen_sym", label = "Gene Symbol: ",
                          choices = dataTables$gene_exp_names)
        
        # load first expression matrix
        dataTables$gene_letter_choice <- substr(x = dataTables$gene_exp_names[1], start = 1, stop = 1)
        withProgress(message = "Loading gene expression data, please wait",
                                                 expr = dataTables$gene_exp_data <- readRDS("Data/gene.expression.formatted.no.IgG.IgM.A.Rds"))
      }
      
      # load a gene expression matrix depending on the users choice of gene
      observeEvent(c(input$gen_sym, input$gen_plot_type), {
        if (input$gen_sym == "") {
          # do nothing if no input gene symbol
        } else {
          # get first letter of input gene
          gene_let <- substr(x = input$gen_sym, start = 1, stop = 1)
          
          # for DEBUGGING
          # print(paste("symbol", input$gen_sym))
          # print(paste("Gene letter:", gene_let))
          # print(paste("data table gene letter:", dataTables$gene_letter_choice))
          
          # if no gene input, default to genes that start with "A"
          if (dataTables$gene_letter_choice == "A") {
            # do nothing
          }
          if (dataTables$gene_letter_choice != gene_let) {
            # determine which expression matrix need to load
            gene_list <- list.files(path = "Data/", pattern = "gene.expression*")
            gen_idx <- which(grepl(x = gene_list, pattern = paste(".", gene_let, ".Rds", sep = "")))
            # update gene letter choice
            dataTables$gene_letter_choice <- gene_let
            # read in gene expression data
            withProgress(message = "Loading gene expression data, please wait",
                         expr = dataTables$gene_exp_data <- readRDS(paste("Data/", gene_list[gen_idx], sep = "")))
          } else {}
          exp_sel <- dataTables$gene_exp_data[dataTables$gene_exp_data$Symbol == input$gen_sym, ]
          # FOR DEBUGGING
          # print(paste("gene matrix:", dim(exp_sel)))
          # print(paste("symbol", input$gen_sym))
          # print(paste("Data/", gene_list[gen_idx], sep = ""))
          if (input$gen_plot_type == "Treatment") {
            output$rna_plot <- renderPlot({
              r1 <- ggplot(exp_sel, aes(x=group_name, y=LogEx, color=Treatment)) +
                geom_boxplot() +
                ylab("Log Normalized Gene Expression") +
                xlab("") +
                facet_wrap(~Treatment, nrow = 1) +
                coord_flip() +
                ggtitle(exp_sel$Symbol[1]) +
                theme(plot.title = element_text(face = "italic"))
              plot(r1)
            })
          }
          else if (input$gen_plot_type == "Cell Type") {
            output$rna_plot <- renderPlot({
              r1 <- ggplot(exp_sel, aes(x=Treatment, y=LogEx, color=group_name)) +
                geom_boxplot() +
                ylab("Log Normalized Gene Expression") +
                xlab("") +
                facet_wrap(~group_name, nrow = 3) +
                coord_flip() +
                ggtitle(exp_sel$Symbol[1]) +
                theme(plot.title = element_text(face = "italic"))
              plot(r1)
            })
          }
        }
        
        
      }) # end of observevent
      
    } # end of tab
  }) # end of observe to click tab
  
  # tab for ADT expression data
  observeEvent(input$tabs,{
    if (input$tabs == "ADT") {
      
      if (is.null(dataTables$adt_exp_data)) {
        # read in gene expression data
        withProgress(expr = adt_df <- readRDS("Data/ADT.expression.formatted.no.IgG.IgM.Rds"),
                               message = "Loading ADT expression data, please wait")
        # update gene choices
        updateSelectInput(session = session, inputId = "adt_sym", label = "ADT Symbol: ",
                          choices = unique(adt_df$Symbol))
        
        # boxplot of gene expression
        observeEvent(c(input$adt_sym, input$adt_plot_type), {
          if (input$adt_sym == "") {
            # do nothing
          } else {
            adt_sel <- adt_df[adt_df$Symbol == input$adt_sym, ]
            if (input$adt_plot_type == "Treatment") {
              output$adt_plot <- renderPlot({
                a1 <- ggplot(adt_sel, aes(x=group_name, y=LogEx, color=Treatment)) + 
                  geom_boxplot() +
                  ylab("Log Normalized ADT Expression") +
                  xlab("") +
                  facet_wrap(~Treatment, nrow = 1) +
                  coord_flip() +
                  ggtitle(adt_sel$Symbol[1])
                plot(a1)
              })
            } else if (input$adt_plot_type == "Cell Type") {
              output$adt_plot <- renderPlot({
                a1 <- ggplot(adt_sel, aes(x=Treatment, y=LogEx, color=group_name)) + 
                  geom_boxplot() +
                  ylab("Log Normalized ADT Expression") +
                  xlab("") +
                  facet_wrap(~group_name, nrow = 3) +
                  coord_flip() +
                  ggtitle(adt_sel$Symbol[1])
                plot(a1)
              })
            }
          }
         
        }) # end of plotting observe event
        # add to reactive value
        dataTables$adt_exp_data <- 1
      }
    }
  })
  
  # tab for pathway response data
  observeEvent(c(input$cell_type, input$path_type, input$view_type), {
    # load appropriate file
    path_files <- list.files("Data/Pathway_Files/", pattern = ".Rds")
    id_path <- which(grepl(x = path_files, pattern = input$path_type) & grepl(x = path_files, pattern = input$cell_type))
    path_res <- readRDS(paste("Data/Pathway_Files/", path_files[id_path], sep = ""))
    
    # plot of enriched pathways 
    if (isolate(input$view_type == "Dotplot")) {
      output$path_plot <- renderPlot({
        clusterProfiler::dotplot(object = path_res, 
                                 showCategory = input$term_num) + 
          ggtitle(input$path_type)
      })
    } else if (isolate(input$view_type == "Barplot")) {
      output$path_plot <- renderPlot({
        barplot(path_res, showCategory = input$term_num) + 
          ggtitle(input$path_type)
      })
    } else if (isolate(input$view_type == "Map")) {
      output$path_plot <- renderPlot({
        emapplot(path_res, showCategory = input$term_num, color="p.adjust") + 
          ggtitle(input$path_type)
      })
    }
    
    # table of enriched pathways
    output$path_table <- DT::renderDataTable({
      DT::datatable(path_res@result)
    })
    
    # download table
    output$Download_Pathways <- downloadHandler(
      filename = function() {
        new_file <- path_files[id_path]
        new_nam <- gsub(x = new_file, pattern = ".Rds", replacement = ".csv")
        paste(new_nam, sep = "")
      },
      content = function(file) {
        write.csv(path_res@result, file)
      }
    )
  })
  
  # tab for trajectory data
  observeEvent(input$tabs,{
    if (input$tabs == "Trajectory") {
        if (is.null(dataTables$dat_obj)) {
          observeEvent(c(input$traj_name), {
            # read in trajectory and expression object
            traj_dat <- list.files("Data/Trajectory/", pattern = ".dataset.Rds")
            traj_mod <- list.files("Data/Trajectory/", pattern = ".model.Rds")
            id_trj <- which(grepl(x = traj_dat, pattern = input$traj_name))
            withProgress(expr = dataTables$dat_obj <- readRDS(paste("Data/Trajectory/", traj_dat[id_trj], sep = "")),
                         message = "Loading in trajectory dataset, please wait...")
            withProgress(expr = dataTables$model_obj <- readRDS(paste("Data/Trajectory/", traj_mod[id_trj], sep = "")),
                         message = "Loading in trajectory model, please wait...")
            
            # update gene entries
            updateSelectInput(session = session, inputId = "traj_marker_id", label = "Choose: ",
                              choices = colnames(dataTables$dat_obj$expression))
          })
          
          observeEvent(c(input$traj_plot_type), {
            # first trajectory plot
            if (isolate(input$traj_plot_type == "Stimulation")) {
              output$traj_plot <- renderPlot({
                withProgress(expr = plot_dimred(
                  dataTables$model_obj, 
                  color_density = "grouping",
                  grouping = dataTables$dat_obj$cell_info$Stimulation,
                  label_milestones = F,
                  alpha_cells = 0.5
                ),
                message = "Updating trajectory plot, please wait...")
              })
            } else if (isolate(input$traj_plot_type == "Pseudotime")) {
              output$traj_plot <- renderPlot({
                withProgress(expr = plot_dimred(dataTables$model_obj, "pseudotime", pseudotime = calculate_pseudotime(dataTables$model_obj)) + ggtitle("Pseudotime"),
                             message = "Calculating and plotting pseudotime, please wait...")
              })
            }
          }) # end of first observe event in traj tab
          
          # second event for second trajectory plot with markers
          observeEvent(c(input$traj_marker_id), {
            output$traj_exp_plot <- renderPlot({
              withProgress(expr = plot_dimred(dataTables$model_obj, feature_oi = isolate(input$traj_marker_id), expression_source = dataTables$dat_obj) + ggtitle(isolate(input$traj_marker_id)),
                           message = "Plotting expression of gene/ADT on trajectory, please wait...")
            })
          })
        }
      
    } # end of if tab
  }) # end of trajectory tab
      
  
})

# run the app
shinyApp(ui, server)
