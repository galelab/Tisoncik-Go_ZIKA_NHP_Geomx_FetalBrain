
library(GSEABase)
library(fgsea)

generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}
results_folder<- "2.Network"
generate_folder(results_folder)
pathsofinterest <- c("GO_NEURON_PROJECTION_GUIDANCE", "GO_NEURON_DEVELOPMENT",
    "GO_NEURON_DIFFERENTIATION", "GO_AXON_DEVELOPMENT", 
    "GO_CELL_MORPHOGENESIS_INVOLVED_IN_NEURON_DIFFERENTIATION", "GO_NEURON_MIGRATION",
    "GO_SYNAPTIC_SIGNALING")

GMTFILE <- "/share/lwhitmo/GeneEnrichmentDBs/msigdb_v7.0_files_to_download_locally/msigdb_v7.0_GMTs/c5.bp.v7.0.symbols.gmt"
pathways <- gmtPathways(GMTFILE)
DEgenes <- read.csv("./1.de_results_full/sigLFCs.csv", header=T)
write("Pathway\tGenes", file.path(results_folder, "Pathwayswithgenesfullpathway.tsv"))
write("Pathway\tGenes", file.path(results_folder, "PathwayswithgenesSiglpathway.tsv"))

sigGenes <- list()
totalgenes <- c()
for (path in pathsofinterest) {
    genes <- pathways[[path]]
    print(length(genes))
    write(paste(path, paste(genes, collapse =","), sep="\t"), 
        file.path(results_folder, "Pathwayswithgenesfullpathway.tsv"), append=T)
    genesint <- intersect(DEgenes$X, genes)
    message(length(genesint), " intersecting with DE genes")
    write(paste(path, paste(genesint, collapse =","), sep="\t"), 
        file.path(results_folder, "PathwayswithgenesSiglpathway.tsv"), append=T)
    sigGenes[[path]] <- genesint
    totalgenes <- c(totalgenes, genesint)
}
totalgenes <- unique(totalgenes) 
genepaths <- list()
for (g in totalgenes) {
    paths <- c()
    for (p in names(sigGenes)) {
        if (g %in% sigGenes[[p]]) {
            paths <- c(paths, p)
        }
    }
    genepaths[[g]]<-unique(paths)
}

message(length(unique(totalgenes)), " total genes")
write("Gene\tPathways", file.path(results_folder, "TotalUniqueDEgenesinanypathofinterest.tsv"))
for (g in names(genepaths)) {
    write(paste(g, paste(genepaths[[g]], collapse =","), sep="\t"), 
        file.path(results_folder, "TotalUniqueDEgenesinanypathofinterest.tsv"), append=T)
}