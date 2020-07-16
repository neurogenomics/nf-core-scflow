#!/usr/bin/env Rscript
# Perform differential gene expression on a SingleCellExperiment Object
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = future::availableCores())

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)
library(cli)

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

# specify options
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")

required$add_argument(
  "--de_table",
  help = "path to de_table.qs",
  metavar = "/dir/de_table.qs",
  required = TRUE
)

required$add_argument(
  "--fc_threshold",
  type = "double",
  default = 1.1,
  metavar = "number",
  help = "Absolute fold-change cutoff for DE [default %(default)s]"
)

required$add_argument(
  "--pval_cutoff",
  type = "double",
  default = 0.05,
  metavar = "number",
  help = "adjusted p-value cutoff for DE report [default %(default)s]"
)

required$add_argument(
  "--n_label",
  help = "Number of top genes to label on volcano plot",
  metavar = "integer",
  required = TRUE
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
#   _____________________________________________________________________

#Reading the de_table output from process scflow_perform_de

res <- qs::qread(args$de_table)

#make output dirs

new_dirs <- c(
  "de_report",
  "de_plot",
  "de_plot_data")

purrr::walk(new_dirs, ~ dir.create(file.path(getwd(), .)))


#generate file name

celltype <- attr(res, "de_params")["celltype"]
de_method <- attr(res, "de_params")["de_method"]
pseudobulk <- attr(res, "de_params")["pseudobulk"]

if (pseudobulk == "Yes") {
  pb_str <- "_pb"
} else if (pseudobulk == "No"){
  pb_str <- ""
}

file_name <- paste0(celltype, "_",
                    de_method, pb_str, "_")

result <- pseudobulk <- attr(res, "de_params")["contrast_name"]

#Generating the report
report_de(res = res,
          fc_threshold = args$fc_threshold,
          pval_cutoff = args$pval_cutoff,
          n_label = args$n_label,
          report_folder_path = file.path(getwd(), "de_report"),
          report_file = paste0(file_name, result, "_scflow_de_report"))
    
#Saving the high resolution plot
p <- .volcano_plot(
  dt = res,
  fc_threshold = args$fc_threshold,
  pval_cutoff = args$pval_cutoff,
  n_label = args$n_label
)

ggplot2::ggsave(file.path(getwd(), "de_plot", 
                          paste0(file_name, result, "_volcano_plot.png")), 
                plot = p,
                width = 247, height = 170, units = "mm", dpi = 600)


#Saving the plot data
plot_data <- p$data

write.table(p$data,
            file = file.path(getwd(), "de_plot_data",
                             paste0(file_name, result, ".tsv")),
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)



