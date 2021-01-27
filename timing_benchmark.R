
###########################
## Function to run roary ##
###########################
run_roary <- function(gffs,    # The gffs to use.
                      p = 1,   # The number of cpus to use.
                      f = ".") # The path to the output directory.
  {
  gffs <- normalizePath(gffs)
  cmd <- paste(
    "roary",
    "-p", p,
    "-f", f,
    paste(gffs, collapse = " "),
    sep = " "
    )
  system(cmd)
  list.files(path = f, full.names = TRUE)
}


####################################################
## Function to generate pagoo raw in-memory input ##
####################################################
#' This function reads roary's output and the gff files
#' and generates the input tables to feed the pagoo()
#' function. 
generate_input <- function(gene_presence_absence_csv, # Roary's csv file. 
                           gffs,                      # List of gff files. 
                           sep = "__")                # Separator.
  {
  
  # Load Biostrings package
	library(Biostrings)

  # Read csv file and parse it to be pagoo compatible.
  df <- read.csv(gene_presence_absence_csv, 
                 header = TRUE,
                 sep = ",", 
                 stringsAsFactors = FALSE, 
                 check.names = FALSE)
  cluster_meta <- df[, c("Gene", "Annotation")]
  colnames(cluster_meta) <- c("cluster", "Annotation")
  dims <- dim(df)
  df <- df[, 15:dims[2]]
  lp <- lapply(df, strsplit, "\t")
  lp <- lapply(lp, setNames, cluster_meta$cluster)
  lths <- lapply(lp, sapply, length)
  sums <- sapply(lths, sum)
  mm <- cbind(unlist(lp, use.names = FALSE), unlist(lapply(lths,
                                                           function(x) {
                                                             unlist(rep(names(x), x), use.names = FALSE)
                                                           }), use.names = FALSE), rep(names(sums), sums))
  mm <- as.data.frame(mm)
  colnames(mm) <- c("gene", "cluster", "org")
  
  # Read gffs.
  names(gffs) <- sub("[.]gff$", "", basename(gffs))
  seqs <- lapply(gffs, function(x) {
    sq <- pagoo:::read_gff(x)
    sq <- sq[which(names(sq) %in% mm$gene)]
  })
  mcls <- lapply(seqs, mcols)
  mcls <- lapply(names(mcls), function(x) {
    mcls[[x]]$org <- x
    mcls[[x]]
  })
  cols <- c("seqid", "type", "start", "end", "strand",
            "product", "org", "locus_tag")
  mcls <- lapply(mcls, function(x) x[, cols])
  mcls <- do.call(rbind, mcls)
  ma <- match(paste(mm$org, mm$gene, sep = sep), paste(mcls$org,
                                                       mcls$locus_tag, sep = sep))
  mm <- cbind(mm, mcls[ma, c("seqid", "type", "start",
                             "end", "strand", "product")])
  
  # Generate pagoo's input list
  ret <- list(data = mm, 
              cluster_meta = cluster_meta, 
              sequences = DNAStringSetList(seqs), 
              sep = sep)

  # Return
  return(ret)
}


######################
## Download dataset ##
######################
#' Download Decano et al. 2019 dataset. 
#' About ~ 7Gb (compressed).

## Download and decompress data
gffdir <- "decano_downing_2019_dataset_gffs/"
if ( !file.exists(gffdir) ){
  system(paste("zenodo_get -d 10.5281/zenodo.3341535 -o", gffdir))
  system(paste("gunzip ", gffdir, "*.gz", sep = ""))
}

# List files
allgffs <- list.files(path = "decano_downing_2019_dataset_gffs/", 
                      pattern = "[.]gff$", 
                      full.names = TRUE)

library(pagoo)
library(parallel)

# Set sample sizes
ns <- c(10, 100, 500, 1000)
fout <- paste("N", ns, sep = "_")

# Seed and sample gff files
set.seed(123)
gffs <- lapply(ns, function(x) sample(allgffs, x) )

# At most 5 processes with 6 cores each
p <- 6
mc.cores <- 5

timings_dfs <- mcmapply(function(gffs, n, f, p) { 

  # Run roary
  roary_out <- run_roary(gffs, p = p, f = f)

  # grep roary's csv file
  gpa <- grep("gene_presence_absence[.]csv$", roary_out, value = TRUE)

  # repeat 10 times each operation
  dfout <- lapply(1:10, function(x){
    
    #' Operation 1: time to load a pagoo object using the 
    #' roary_2_pagoo() function with only the roary's csv 
    #' as input.
    st1 <- system.time(p <- roary_2_pagoo(gpa))
    df1 <- data.frame(Operation = "roary_csv",
                      N = n,
                      Dir = f,
                      User = st1[["user.self"]],
                      System = st1[["sys.self"]],
                      Elapsed = st1[["elapsed"]])

    #' Operation 2: time to generate a pangenome barplot.
    st2 <- system.time(p$gg_barplot())
    df2 <- data.frame(Operation = "barplot",
                      N = n,
                      Dir = f,
                      User = st2[["user.self"]],
                      System = st2[["sys.self"]],
                      Elapsed = st2[["elapsed"]])

    #' Operation 3: time to load a pagoo object using the 
    #' roary_2_pagoo() function including sequences and 
    #' data stored in the gff files.
    st3 <- system.time(p <- roary_2_pagoo(gpa, gffs))
    df3 <- data.frame(Operation = "roary_csv_gff",
                      N = n,
                      Dir = f,
                      User = st3[["user.self"]],
                      System = st3[["sys.self"]],
                      Elapsed = st3[["elapsed"]])

    #' Operation 4: time to return core sequences.
    st4 <- system.time(p$core_seqs_4_phylo())
    df4 <- data.frame(Operation = "core_seq",
                      N = n,
                      Dir = f,
                      User = st4[["user.self"]],
                      System = st4[["sys.self"]],
                      Elapsed = st4[["elapsed"]])
    rm(p)



    #' Operation 5: time to load a pagoo object using the
    #' pagoo() function. Similar to Operation 3, but all
    #' the information is already in memory (in R objects) 
    #' instead of in files. 
    input <- generate_input(gpa, gffs)
    st5 <- system.time(pagoo(data = input$data, cluster_meta = input$cluster_meta,
                                  sequences = input$sequences, sep = input$sep))
    df5 <- data.frame(Operation = "preloaded_seq",
                      N = n,
                      Dir = f,
                      User = st5[["user.self"]],
                      System = st5[["sys.self"]],
                      Elapsed = st5[["elapsed"]])

    # Bind data.frames and return
    rbind(df1, df2, df3, df4, df5)
  })
  return(dfout)

}, gffs = gffs, 
   n = ns, 
   f = fout, 
   MoreArgs = list(p = p), 
   SIMPLIFY = FALSE, 
   mc.cores = mc.cores)

# Save list of data.frames
saveRDS(timings_dfs, file = "timings.RDS")

# Write tsv
x <- do.call(rbind, unlist(timings_dfs, recursive = FALSE))
write.table(x, file = "timings.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

# Plot
library(ggplot2)
g <- ggplot(x, aes(y=Elapsed, x=N, color=Operation, group=N)) +
  geom_point() + 
  facet_wrap(~Type, scales="free_y", nrow=2) + 
  stat_summary(fun=median, geom="line", aes(group=Operation)) + 
  theme_light() + 
  scale_color_discrete(labels=c("$gg_barplot()", 
                                "$core_seqs_4_phylo()", 
                                "pagoo()", 
                                "roary_2_pagoo()", 
                                "roary_2_pagoo() + gff"), 
                       guide = guide_legend(reverse=T)) + 
  ylab("Elapsed time (s)") + 
  xlab("Number of organisms")

# Save plot
ggsave("bench_pagoo.pdf", plot=g)
