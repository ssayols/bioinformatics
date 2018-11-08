##
## Interactive visualization of gene expression on the t-SNE plot
##
## What do you need? An RData file with 3 objects called like this:
##   - sce: a single cell experiment object, like the one created in the Rmd file.
##   - genes: a data frame with gene annotation, consisting of 2 columns: gene_id and gene_name.
##   -        The gene_id is the key field of the sce object (rownames(sce))
##   - top.hvg: subset of the top variable genes (gene_id) to be used in the t-SNE.
##
library(shinydashboard)
library(DT)
library(plotly)
library(ggplot2)
library(scater)

##
## Define UI for application
##
ui <- dashboardPage(
  dashboardHeader(title="Gut-Stomach scRNA-seq"),
  dashboardSidebar(disable=TRUE), 
  dashboardBody(
    fluidRow(
      column(width=9,
        box(width=NULL, solidHeader=TRUE, plotlyOutput("tsnePlot"))
      ),
      column(width=3,
        box(width=NULL, status="warning",
            fluidRow(
              column(6, selectInput("dataset", label="dataset", choices=list())),
              column(6, numericInput("perp", label="perplexity", min=1, max=500, value=5))
            ),
            DT::dataTableOutput("genes"),
            hr(),
            downloadLink('downloadPlot', 'Download figure'),
            fluidRow(
              column(6, numericInput("width", label="width", min=1, max=20, value=7)),
              column(6, numericInput("height", label="height", min=1, max=20, value=7))
            )
        )
      )
    )
  )
)

##
## Define server logic
##
server <- function(input, output, session) {
 
  session$onSessionEnded(stopApp)
  
  # show the available datasets to the interface
  observe({
    updateSelectInput(session, "dataset", label="dataset", choices=list.files(pattern="\\.RData$"))
  })
  
  # load pre-rendered data
  data <- reactiveValues()
  observe({
    if(file.exists(input$dataset)) {
      withProgress(message="Reading data...", value=0, {
        x <- local({
          load(input$dataset)
          list(sce=sce, genes=genes, top.hvg=top.hvg)
        })
        data$sce <- x$sce
        data$genes <- x$genes
        data$top.hvg <- x$top.hvg
        rm(x)
      })
    }
  })
  
  # load the table with gene names
  output$genes <- DT::renderDataTable({
    if(!is.null(data$genes))
      DT::datatable(data$genes, selection="single", options=list(pageLength=5))
  })
  
  # the plot
  tSNEplot <- reactiveValues(p=NULL)
  observe({
    if(length(input$genes_rows_selected) == 1) {
      withProgress(message="Calculating tSNE plot", value=0, {
        p <- plotTSNE(data$sce, exprs_values="norm_exprs",
                 colour_by=data$genes$gene_id[input$genes_rows_selected],
                 shape_by="sort", perplexity=input$perp, rand_seed=100, feature_set=data$top.hvg) +
          ggtitle(data$genes$gene_name[input$genes_rows_selected]) +
          theme(legend.title=element_blank())
        p$data$label <- rownames(p$data)
        p$mapping$label <- as.name("label")
        tSNEplot$p <- p
      })
    } else {
      tSNEplot$p <- NULL
    }
  })
  
  # display the plot
  output$tsnePlot <- renderPlotly({
    if(!is.null(tSNEplot$p))
      ggplotly(tSNEplot$p, tooltip=c("colour_by", "shape_by", "label"))
    else
      ggplotly(ggplot() + geom_blank() + theme_minimal())
  })
  
  # download the plot
  output$downloadPlot <- downloadHandler(
    filename=function() {
      if(length(input$genes_rows_selected) == 1)
        paste0(data$genes$gene_name[input$genes_rows_selected], ".pdf")
      else
        "empty.pdf"
    },
    content=function(con) {
      pdf(con, width=input$width, height=input$height)
      if(!is.null(tSNEplot$p))
        p <- tSNEplot$p
      else
        p <- ggplot() + geom_blank() + theme_minimal()
      print(p)
      dev.off()
    }
  ) 
}

# Run the application 
shinyApp(ui = ui, server = server)
