#' ---
#' title: "Call GeneMania and get the network for that gene"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float:
#'       collapsed: false
#'       smooth_scroll: false
#' ---
# Note: call with Rscript ./get_gm_net.R targets=genes.txt \
#                                        outdir=get_gm_net \
#                                        org=hs \
#                                        limit=20 \
#                                        cores=4 \
#                                        gm_path=~/bin//genemania/genemania-cytoscape-plugin-2.18.jar \
#                                        gm_data_path=~/bin//genemania/gmdata-2017-07-13-core/
library(XML)
library(igraph)
library(parallel)

#' get arguments from the command line
parseArgs <- function(args, string, default=NULL, convert="as.character") {
  if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  if(!is.null(default)) default else do.call(convert,list(NA))
}

args    <- commandArgs(T)
TARGETS <- parseArgs(args, "targets=", default="targets.txt")
OUTDIR  <- parseArgs(args, "outdir=" , default="get_gm_net")
ORG     <- parseArgs(args, "org="    , default="hs")
LIMIT   <- parseArgs(args, "limit="  , default=20, convert="as.numeric")
CORES   <- parseArgs(args, "cores="  , default=4, convert="as.numeric")
GM_PATH <- parseArgs(args, "gm_path=", default="/home/sayolspu/src/toxpred/src/networks/genemania-cytoscape-plugin-2.18.jar")
GM_DATA_PATH <- parseArgs(args, "gm_data_path=", default="/home/sayolspu/src/toxpred/src/networks/genemania-data/gmdata-2017-07-13-core/")

runstr <- paste("Call with:",
                "Rscript ./get_gm_net.R targets=genes.txt \\", 
                "                       outdir=get_gm_net \\",
                "                       org=hs \\",
                "                       limit=20 \\",
                "                       cores=4 \\",
                "                       gm_path=~/bin//genemania/genemania-cytoscape-plugin-2.18.jar \\",
                "                       gm_data_path=~/bin//genemania/gmdata-2017-07-13-core/",
                sep="\n")
if(!file.exists(TARGETS)) stop("targets file not found.\n", runstr)
if(!file.exists(GM_PATH)) stop("GeneMania binary genemania-cytoscape-plugin-2.18.jar not found.\n", runstr)
if(!dir.exists(GM_DATA_PATH)) stop("GM_DATA_PATH not found.\n", runstr)
if(!ORG %in% c("hs", "mm")) stop("org must be 'hs' or 'mm'.\n", runstr)
ORG <- switch(ORG,
              "hs"="H. Sapiens",
              "mm"="M. Musculus")

targets <- read.csv(TARGETS)

#' ## GENEMANIA PART
#' predict target's network with genemania (external jar)
gmn <- mclapply(targets[, 1], function(x) {
  # temp file with the query:
  # S. Cerevisiae
  # CDC27 APC11 APC4  XRS2  RAD54 APC2  RAD52 RAD10 MRE11 APC5
  # preferred
  # 150
  # bp
  input   <- tempfile()
  output  <- tempfile()
  writeLines(paste(ORG,
                   paste(toupper(x), collapse="\t"),
                   "preferred",  # coexp (Co-expr), gi (Genetic int), pi (Physical int)
                   LIMIT,  # related-gene-limit
#                   "bp",  # Networks are weighted in an attempt to reproduce GO BP co-annotation patterns
                   sep="\n"),
             input)  

  # run command
  cmd <- paste("java -Xmx1800M",
               paste("-cp", GM_PATH, "org.genemania.plugin.apps.QueryRunner"),
               paste("--data", GM_DATA_PATH),
               "--out xml",
               "--threads", CORES,
               "--results", output,
               input)
  dir.create(output)
  try(system(cmd))

  # parse output
  f <- file.path(output, paste0(basename(input), "-results.report.xml"))
  nw <- matrix(nrow=0, ncol=0)
  if(file.exists(f)) {
    tryCatch({
      gm <- xmlToList(xmlParse(f))
      nw <- do.call(rbind, gm[["results"]][["interactions"]])
      nw <- as_adjacency_matrix(graph_from_edgelist(nw[, 1:2, drop=FALSE]), sparse=FALSE)
    }, error=function(e) invisible(e))
  }

  nw
}, mc.cores=CORES)

names(gmn) <- targets[, 1]
saveRDS(gmn, file=file.path(OUTDIR, "get_gm_net.RDS"))

#' Extract adjacency matrixes and graphml to be open in cytoscape
try(dir.create(OUTDIR))
invisible(
  Map(function(x, name) {
    write.csv(x, file=paste0(file.path(OUTDIR, name), ".csv"))
  }, gmn, names(gmn))
)

#' Save as graphml files to open in Cytoscape
invisible(
  Map(function(x, name) {
    write_graph(graph_from_adjacency_matrix(x), file=paste0(file.path(OUTDIR, name), ".graphml"), "graphml")
  }, gmn, names(gmn))
)

