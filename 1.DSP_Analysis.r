######### --Libraries--#########
set.seed(1)
options(mc.cores = 8)
library(GSEABase)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(openxlsx)
library(reshape2)
library(knitr)
library(cowplot)
library(edgeR)
library(pheatmap)
library(ClassDiscovery)
library(plotrix)
library(SetRank)
library(GeneSets.Homo.sapiens)
library(nlme)
library(lme4)
library(data.table)
library(ggrepel)
library(stringr)

######### --Functions --########
SetRankRun=FALSE
if (isTRUE(SetRankRun)) {
    # converters
    symbol2EntrezID <- createIDConverter(
        "Homo.sapiens", "SYMBOL",
        "ENTREZID"
    )
    IDConverter <- createIDConverter(
        "Homo.sapiens", "ENTREZID",
        "SYMBOL"
    )
}

scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

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

theme_Publication <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.6, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face = "italic"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

QC_histogram <- function(targetfile, qc_param="Nuclei",
                         color_by = NULL, results_path=getwd(),
                         thr = NULL,
                         scale_trans = NULL) {
    plt <- ggplot(targetfile,
                  aes_string(x = qc_param, color =color_by)) +
        geom_histogram(bins = 50,fill="white", alpha=0.5, position="identity") +
        # geom_vline(xintercept = thr, lty = "dashed", color = "black") +
        theme_bw() + guides(fill = "none") + scale_color_brewer(palette="Dark2")
        # facet_wrap(as.formula(paste("~", facet_by)), nrow = 2) +
        labs(x = qc_param, y = "Segments, #")
    ggsave(plt, file=file.path(results_path, paste0(qc_param, ".png")), 
        width=5, height=3.6, units="in", dpi=300)
}

pre_norm_figure <- function(Statdf, Statdf_m, results_path, figname="NormalizationPlot.png") {

    plt1 <- ggplot(
        Statdf_m,
        aes(x = Value, fill = Statistic)
    ) +
        geom_histogram(bins = 20) +
        theme_bw() +
        scale_x_continuous(trans = "log2") +
        facet_wrap(~Annotation, nrow = 1) +
        scale_fill_brewer(palette = 3, type = "qual") +
        labs(x = "Counts", y = "Segments, #")

    plt2 <- ggplot(
        Statdf,
        aes(x = NegProbe, y = Q3, color = Annotation)
    ) +
        geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
        geom_point() +
        guides(color = "none") +
        theme_bw() +
        scale_color_manual(values = c("red", "darkblue", "purple")) +
        scale_x_continuous(trans = "log2") +
        scale_y_continuous(trans = "log2") +
        theme(aspect.ratio = 1) +
        labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

    plt3 <- ggplot(
        Statdf,
        aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)
    ) +
        geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
        geom_point() +
        theme_bw() +
        scale_color_manual(values = c("red", "darkblue", "purple")) +
        scale_x_continuous(trans = "log2") +
        scale_y_continuous(trans = "log2") +
        theme(aspect.ratio = 1) +
        labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

    btm_row <- plot_grid(plt2, plt3,
        nrow = 1, labels = c("B", ""),
        rel_widths = c(0.43, 0.57)
    )
    plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
    ggsave(file.path(results_path, figname), height = 8, width = 8, dpi = 400)
}

Q3_norm <- function(Statsdf, counts, cpm=FALSE) {
    countsraw <- counts
    geomeanq3 <- ngeoMean(Statsdf$Q3)
    for (seg in colnames(counts)) {
        q3val <- Statsdf[seg, "Q3"]
        counts[, seg] <- counts[, seg] / q3val
        counts[, seg] <- counts[, seg] * geomeanq3
        if (isTRUE(cpm)) {
            countsum <- sum(countsraw[, seg])
            counts[, seg] <- (counts[, seg]/countsum) *1000000

        }
    }
    return(counts)
}

pca_fun <- function(exprs, labels, results_path,
                    base_file_name, target_columns,
                    figres = 100, size = 1, pca = FALSE) {
    library(factoextra)
    # Run PCA/SVD reduction
    if (isFALSE(pca)) {
        pca <- prcomp(t(exprs))
    }
    E <- get_eig(pca)
    cx <- sweep(t(exprs), 2, colMeans(t(exprs)), "-")
    sv <- svd(cx)


    vizualize_pca(
        file.path(results_path, paste0("svd_", base_file_name)),
        sv$u, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, E, size
    )
    vizualize_pca(
        file.path(results_path, paste0("pca_", base_file_name)),
        pca$x, labels[, target_columns[1]],
        labels[, target_columns[2]],
        figres, E, size
    )
    vizualize_scree_plot(
        file.path(
            results_path,
            paste0("scree_", base_file_name)
        ), pca, figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    is_pc1_0 <- loadingscores$PC1 > 0
    is_pc2_0 <- loadingscores$PC2 > 0

    loadingscores <- loadingscores[is_pc1_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC1)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc1", base_file_name, ".txt")),
        loadingscores["PC1"], figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc2_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC2)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc2", base_file_name, ".txt")),
        loadingscores["PC2"], figres
    )
    return(pca)
}

umap_fun <- function(exprs, labels, results_path,
                     base_file_name, target_columns,
                     figres = 100, size = 1, UMAP = FALSE) {
    library(umap)
    # Runs default paramaters of umap
    if (isFALSE(UMAP)) {
        UMAP <- umap(t(exprs))
    }
    vizualize_umap(
        file.path(results_path, paste0("umap_", base_file_name)),
        UMAP$layout, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, size
    )

    return(UMAP)
}

vizualize_umap <- function(plot_file, U, class1, class2, figres, size) {
    # Vizualize umap reduction
    library(Polychrome)
    minx <- min(U[, 1])
    maxx <- max(U[, 1])
    miny <- min(U[, 2])
    maxy <- max(U[, 2])
    if (length(levels(factor(class2))) <= 3) {
        if (length(levels(factor(class1))) <= 6) {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                scale_color_brewer(palette="Dark2") +
                scale_fill_brewer(palette="Dark2") +
                theme(legend.position = "right")
        } else {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                scale_color_brewer(palette="Dark2") +
                scale_fill_brewer(palette="Dark2") +
                theme(legend.position = "right") +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    } else {
        P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
        if (length(levels(factor(class1))) <= 6) {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                theme(legend.position = "right")
        } else if (length(levels(factor(class1))) > 6) {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                xlim(minx, maxx) + ylim(miny, maxy) +
                theme(legend.position = "right") +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    }
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = figres)
}

vizualize_pca <- function(plot_file, PCA, class1, class2, figres, E, size) {
    # Vizualize PCA  results
    library(Polychrome)
    minx <- min(PCA[, 1])
    maxx <- max(PCA[, 1])
    miny <- min(PCA[, 2])
    maxy <- max(PCA[, 2])
    if (length(levels(factor(class2))) <= 3) {
        if (length(levels(factor(class1))) <= 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                scale_color_brewer(palette="Dark2") +
                scale_fill_brewer(palette="Dark2") +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = "right")
        } else {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                scale_color_brewer(palette="Dark2") +
                scale_fill_brewer(palette="Dark2") +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = "right") + scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    } else {
        P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
        if (length(levels(factor(class1))) <= 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = "right") +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36))
        } else if (length(levels(factor(class1))) > 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = "right") +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    }
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = 300)
}

vizualize_scree_plot <- function(plot_file, PCA, figres) {
    # Vizualize principle component variation results
    scree.plot <- fviz_eig(PCA, addlabels = TRUE, hjust = -0.3)
    png(plot_file, width = 7, height = 6, units = "in", res = figres)
    print(scree.plot)
    dev.off()
}

save_loading_scores <- function(write_file, df, figres) {
    # Save list of genes that have a positive effect on variation of principle
    # component 1 and 2 sorted from most influential
    write.table(df, file = write_file)
}

gene_enrichment <- function(genes, results_folder, cluster) {
    inputGenes <- symbol2EntrezID(genes)
    network <- setRankAnalysis(inputGenes, collection,
        use.ranks = FALSE,
        setPCutoff = 0.01,
        fdrCutoff = 0.05
    )

    generate_folder(results_folder)
    #### IMPORTANT OUTPUT INFORMATION###
    # SetRank value -  value reflects the prominence of a gene set in a gene set network (based on PageRank algorithm developed by google)
    # pSetRank value - significance of the SetRank value
    exportSingleResult(network, inputGenes,
        collection, paste0(results_folder, "/de_unranked_", cluster),
        IDConverter = IDConverter
    )
    # png(file.path(results_folder, paste0("de_unranked_", cluster,".png")), res=100)
    # plot(network, layout = layout.spring)
    # dev.off()
    return(network)
}

run_gsea_and_ora <- function(finalrankedgenes, gmt.file, universe, region, results_folder) {
    library(fgsea)
    pathways <- gmtPathways(gmt.file)
    print(head(finalrankedgenes))
    fgseaRes <- fgsea(
        pathways = pathways,
        stats = finalrankedgenes,
        minSize = 15,
        maxSize = 1000
    )
    message("done gsea")
    print(fgseaRes)
    fgseaRes <- fgseaRes[order(pval), ]
    fgseaResSig <- fgseaRes[fgseaRes$padj <= 0.05, ]
    fgseaResdf <- as.data.frame(fgseaRes)
    fwrite(as.data.frame(fgseaResdf), file.path(results_folder, paste0(region, "_GSEA_allResults.csv")))
    fwrite(as.data.frame(fgseaResSig), file.path(results_folder, paste0(region, "_GSEA_allResultsSig.csv")))
    message("done writing gsea")

    foraRes <- fora(pathways, names(finalrankedgenes), universe, minSize = 5, maxSize = Inf)
    foraRes <- foraRes[order(pval), ]
    foraResSig <- foraRes[foraRes$padj <= 0.05, ]
    fwrite(as.data.frame(foraRes), file.path(results_folder, paste0(region, "_ORA_allResults.csv")))
    fwrite(as.data.frame(foraResSig), file.path(results_folder,paste0(region, "_ORA_allResultsSig.csv")))
    message("done ora")

    # for (path in hmpathways) {
    #     GENES <- foraRes[foraRes$pathway == path, "overlapGenes"][[1]][[1]]
    #     temp <- as.matrix(as.data.frame(deHGNCpvsnp[rownames(deHGNCpvsnp) %in% GENES, ]))
    #     # --heatmap
    #     if (length(GENES) > 1) {
    #         png(file.path(results_folder, paste0("heatmap_", path, ".png")), width = 10, height = 10, units = "in", res = 300)
    #         # par(mar = c(4, 4, -1, 2))
    #         global_modules <- heatmap.L.4(temp,
    #             figmargins = c(10, 10),
    #             cutoff = 1, distmethod = "euclidean", cexcol = 2,
    #             clustermethod = "ward.D2", clusterdim = "row", labRow = rownames(temp),
    #             colsep = c(7), labCol = c("D0", "D3", "D7", "D126", "D129", "D133", "PreChal")
    #         )
    #         dev.off()
    #     }
    # }
    return(list("ORA" = foraRes, "sigGSEA" =fgseaResSig, "GSEA" = fgseaRes))
}

diff_gene_exprs <- function(exprs, target, modelFormula, 
        subset=NULL, contrast.var=NULL, contrast.levels=NULL, 
        pairwise=TRUE, nCores=8) {
    library(nlme)
    library(lme4)
    deFunc <- function(i, groupVar, pDat, modelFormula, exprs, 
        pairwise = TRUE) {
        dat <- data.frame(expr = exprs[i,], pDat)
        if (length(colnames(pDat))>1) {
            lmOut <- suppressWarnings(lmerTest::lmer(modelFormula, 
                dat))
        } else {
            message("STATUS: linear model")
            lmOut <- suppressWarnings(lme4::lm(
                modelFormula,
                dat
            ))
        }
        if (pairwise == FALSE) {
            lsm <- lmerTest::ls_means(lmOut, which = groupVar, 
                pairwise = FALSE)
        }
        else {
            lsm <- lmerTest::ls_means(lmOut, which = groupVar, 
                pairwise = TRUE)
        }
        lmOut <- matrix(stats::anova(lmOut)[groupVar, "Pr(>F)"], 
            ncol = 1, dimnames = list(groupVar, "Pr(>F)"))
        lsmOut <- matrix(cbind(lsm[, "Estimate"], lsm[, "Pr(>|t|)"]), 
            ncol = 2, dimnames = list(gsub(groupVar, "", 
                rownames(lsm)), c("Estimate", "Pr(>|t|)")))
        return(list(anova = lmOut, lsmeans = lsmOut))
    }
    results <- c()
    mTerms <- all.vars(modelFormula)
    print (mTerms)
    if (!contrast.var %in% mTerms) {
        stop("Error 1: contrast.var needs to be defined as fixed effect in the model.\n")
    }
    if (is.null(contrast.var)) {
        stop("Error 2: contrast.var needs to be defined as fixed effect in the model.\n")
    }
    target <- as.data.frame(target)
    for (i in colnames(target)) {
        if (is.null(contrast.levels)) {
            pDat[, i] <-factor(target[, i])
        } else {
            if (i==contrast.var) {
                target[, i] <- factor(target[, i], levels=contrast.levels)
            } else {
                target[, i] <- as.factor(target[, i])
            }
        }
    }
    if (!is.null(subset)) {
        for(fac in unique(target[, subset])) {
            print(fac)
            ind <- target[, subset] == fac
            print (ind)
            if (length(mTerms)==1) {
                message("STATUS: Only one term ")
                pDat <- target[ind, mTerms ]
                names(pDat) <- rownames(target[ind, ])
                pDat <- as.data.frame(pDat)
                colnames(pDat) <- mTerms
            } else { 
                message("STATUS: Multiple terms ")
                pDat <- target[ind, mTerms ]
            }
            exprstmp <- exprs[, rownames(pDat)]
            print("expression")
            print(dim(exprstmp))
            mixedOut <- parallel::mclapply(rownames(exprstmp), 
                deFunc,  contrast.var, pDat, formula(paste("expr", 
                    as.character(modelFormula)[2], sep = " ~ ")),
                exprstmp, mc.cores = nCores)
            print(head(mixedOut))
            mixedOut <- rbind(array(lapply(mixedOut, function(x) x[["anova"]])), 
                array(lapply(mixedOut, function(x) x[["lsmeans"]])))
            colnames(mixedOut) <- rownames(exprstmp)
            rownames(mixedOut) <- c("anova", "lsmeans")
            r_test <- do.call(rbind, mixedOut["lsmeans", ])
            tests <- rownames(r_test)
            r_test <- as.data.frame(r_test)
            r_test$Contrast <- tests
            r_test$Gene <- 
                unlist(lapply(colnames(mixedOut),
                            rep, nrow(mixedOut["lsmeans", ][[1]])))
            r_test$Subset <- fac
            r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
            r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate",
                                "Pr(>|t|)", "FDR")]
            results <- rbind(results, r_test)
        }
    } else { 
        pDat <- target[ mTerms ]
        exprstmp <- exprs[, rownames(pDat)]
        mixedOut <- parallel::mclapply(rownames(exprstmp), 
            deFunc,  contrast.var, pDat, formula(paste("expr", 
                as.character(modelFormula)[2], sep = " ~ ")),
            exprstmp, mc.cores = nCores)
        mixedOut <- rbind(array(lapply(mixedOut, function(x) x[["anova"]])), 
            array(lapply(mixedOut, function(x) x[["lsmeans"]])))
        colnames(mixedOut) <- rownames(exprstmp)
        rownames(mixedOut) <- c("anova", "lsmeans")
        r_test <- do.call(rbind, mixedOut["lsmeans", ])
        tests <- rownames(r_test)
        r_test <- as.data.frame(r_test)
        r_test$Contrast <- tests
        r_test$Gene <- 
                unlist(lapply(colnames(mixedOut),
                            rep, nrow(mixedOut["lsmeans", ][[1]])))
        r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
        r_test <- r_test[, c("Gene", "Contrast", "Estimate",
                            "Pr(>|t|)", "FDR")]
        results <- rbind(results, r_test)    
        # use lapply in case you have multiple levels of your test factor to
        # correctly associate gene name with it's row in the results table
    }
    return(results)
}

fitPoisBG <- function(object, id, iterations = 10, tol = 1e-3, size_scale = c("sum", "first")) {

        size_scale <- match.arg(size_scale)

        uniid <- unique(as.character(id))
        n_feature <- NROW(object)
        n_sample <- NCOL(object)
        ind_na <- which(is.na(object), arr.ind = TRUE)
        # colnames(object) <- id
        sizefact <- apply(object, 2, mean, na.rm = TRUE)

        if (size_scale == "first") {
            scale_fac <- sizefact[1]
        } else if (size_scale == "sum") {
            scale_fac <- sum(sizefact)
        }


        sizefact <- sizefact / scale_fac

        sizefact0 <- sizefact
        sizefact_mat <- matrix(rep(sizefact, n_feature), n_feature, n_sample, byrow = TRUE)
        sizefact_mat[ind_na] <- NA
        featfact0 <- matrix(0, n_feature, length(uniid))
        for (iter in seq_len(iterations)) {
            featfact <- sapply(uniid, function(x) {
                apply(object[, x == id, drop = FALSE], 1, sum, na.rm = TRUE) /
                    apply(sizefact_mat[, x == id, drop = FALSE], 1, sum, na.rm = TRUE)
            })

            featfact_mat <- featfact[, id]
            featfact_mat[ind_na] <- NA

            sizefact <- apply(object, 2, sum, na.rm = TRUE) / apply(featfact_mat, 2, sum, na.rm = TRUE)

            if (size_scale == "first") {
                scale_fac <- sizefact[1]
            } else if (size_scale == "sum") {
                scale_fac <- sum(sizefact)
            }

            sizefact <- sizefact / scale_fac


            sizefact_mat <- matrix(rep(sizefact, n_feature), n_feature, n_sample, byrow = TRUE)
            sizefact_mat[ind_na] <- NA

            message(sprintf(
                "Iteration = %s, squared error = %e",
                iter,
                sum((sizefact - sizefact0)^2) + sum((featfact - featfact0)^2)
            ))

            if (sum((sizefact - sizefact0)^2) + sum((featfact - featfact0)^2) < tol) {
                break
            }

            sizefact0 <- sizefact
            featfact0 <- featfact
        }
        message("Model converged.")

        return(list(
            sizefact = sizefact,
            featfact = featfact,
            countmat = object,
            id = id
        ))
}

diagPoisBG <- function(object, padj = FALSE, padj_method = "BH",
 cutoff = 1e-6, generate_ppplot = TRUE, results_folder=getwd()) {
    countmat <- object$countmat
    if (NCOL(object$featfact) == 1) {
        countmat_expected <- (object$featfact %*% t(object$sizefact))
    } else {
        countmat_expected <- sweep(object$featfact[, object$id], 2, object$sizefact, FUN = "*")
    }

    countmat <- as.matrix(countmat)

    lowtail_prob1 <- ppois(q = countmat, lambda = countmat_expected)
    lowtail_prob2 <- ppois(q = countmat - 1, lambda = countmat_expected)

    lowtail_prob <- (lowtail_prob1 + lowtail_prob2) / 2

    uptail_prob <- 1 - lowtail_prob

    # simualte data (do it only once)
    countmat_simu <- t(apply(countmat_expected, 1, function(x) rpois(rep(1, ncol(countmat)), lambda = x)))
    lowtail_prob1_simu <- ppois(q = countmat_simu, lambda = countmat_expected)
    lowtail_prob2_simu <- ppois(q = countmat_simu - 1, lambda = countmat_expected)
    lowtail_prob_simu <- (lowtail_prob1_simu + lowtail_prob2_simu) / 2

    if (generate_ppplot) {
        y <- sort(lowtail_prob, na.last = TRUE)
        y_simu <- sort(lowtail_prob_simu, na.last = TRUE)
        png(file.path(results_folder, "dispersion_plot.png"))
        plot(y_simu, y, ylim = c(0, 1), xlab = "Empircal CDF from simulated data", ylab = "Empircal CDF", main = "Poisson model")
        graphics::abline(a = 0, b = 1)
        dev.off()
    }


    disper <- mean((countmat - countmat_expected)^2 / countmat_expected, na.rm = TRUE)

    if (padj) {
        lowtail_prob <- matrix(p.adjust(lowtail_prob, method = padj_method), nrow = nrow(lowtail_prob), ncol = ncol(lowtail_prob))
        uptail_prob <- matrix(p.adjust(uptail_prob, method = padj_method), nrow = nrow(uptail_prob), ncol = ncol(uptail_prob))
    }



    low_outlier <- which(lowtail_prob < cutoff, arr.ind = TRUE)
    up_outlier <- which(uptail_prob < cutoff, arr.ind = TRUE)

    return(list(
        lowtail_prob = lowtail_prob,
        uptail_prob = uptail_prob,
        disper = disper,
        outlier = list(
            low_outlier = low_outlier,
            up_outlier = up_outlier
        )
    ))
}

BGScoreTest_sp <- function(object, BGmod, adj = 1, probenum, removeoutlier = FALSE, useprior = FALSE) {
    id <- BGmod$id
    uniid <- unique(as.character(id))

    if (removeoutlier == TRUE) {
        #   boxobj <- apply(BGmod$featfact, 2, function(x) boxplot(x, plot = FALSE))

        featfact <- apply(BGmod$featfact, 2, function(x) {
            boxobj <- graphics::boxplot(x, plot = FALSE)
            message(sprintf("%s negative probes are removed prior to the score test.", length(boxobj$out)))
            x[which(x %in% boxobj$out)] <- NA
            x
        })
    } else {
        featfact <- BGmod$featfact
    }

    sizefact <- BGmod$sizefact


    if (useprior == FALSE) {
        if (missing(probenum)) {
            prodfact <- lapply(uniid, function(x) sizefact[x == id] * mean(adj * featfact[, x], na.rm = TRUE))
            names(prodfact) <- uniid
            #   scores <- apply(countmat, 2, function(x) sum(x - ab)/sqrt(sum(ab)))
            scores_sp <- sapply(uniid, function(x) apply(object[, x == id, drop = FALSE], 1, function(y) sum(y - prodfact[[x]]) / sqrt(sum(prodfact[[x]]))))
        } else {
            if (is.null(names(probenum))) names(probenum) <- rownames(object)
            scores_sp <- sapply(names(probenum), function(feat) {
                prodfact <- lapply(uniid, function(x) sizefact[x == id, drop = FALSE] * mean(probenum[feat] * featfact[, x], na.rm = TRUE))
                names(prodfact) <- uniid
                #   scores <- apply(countmat, 2, function(x) sum(x - ab)/sqrt(sum(ab)))
                sapply(uniid, function(x) sum(object[feat, x == id, drop = FALSE] - prodfact[[x]]) / sqrt(sum(prodfact[[x]])))
            })
            scores_sp <- t(scores_sp)
        }
    } else {
        if (missing(probenum)) {
            featfact0 <- colMeans(adj * featfact, na.rm = TRUE)
            sigma <- apply(adj * featfact, 2, var, na.rm = TRUE) / featfact0^2
            deno <- lapply(uniid, function(x) (sizefact[x == id] * sigma[x] * featfact0[x] + 1) * featfact0[x])
            names(deno) <- uniid

            scores_sp <- sapply(uniid, function(x) apply(object[, x == id, drop = FALSE], 1, function(y) sum((y - sizefact[x == id] * featfact0[x]) / deno[[x]]) / sqrt(sum(sizefact[x == id] / deno[[x]]))))
            #  scores <- apply(scores2, 1, mean)
        } else {
            if (is.null(names(probenum))) names(probenum) <- rownames(object)
            scores_sp <- sapply(names(probenum), function(feat) {
                featfact0 <- colMeans(probenum[feat] * featfact, na.rm = TRUE)
                sigma <- apply(probenum[feat] * featfact, 2, var, na.rm = TRUE) / featfact0^2
                deno <- lapply(uniid, function(x) (sizefact[x == id] * sigma[x] * featfact0[x] + 1) * featfact0[x])
                names(deno) <- uniid

                sapply(uniid, function(x) sum((object[feat, x == id, drop = FALSE] - sizefact[x == id] * featfact0[x]) / deno[[x]]) / sqrt(sum(sizefact[x == id] / deno[[x]])))
            })

            scores_sp <- t(scores_sp)
        }
    }

    pvalues <- pnorm(scores_sp, lower.tail = FALSE)

    return(list(
        pvalues = pvalues,
        scores_sp = scores_sp
    ))
}

fitNBth <-  function(object, features_high, probenum, sizefact_BG, sizefact_start = sizefact_BG, size_scale = c("sum", "first"), threshold_start, threshold_fix = FALSE, tol = 1e-7, iterations = 8,
start_para = c(threshold_start, 1), lower_sizefact = 0, lower_threshold = threshold_start / 5) {
    size_scale <- match.arg(size_scale)
    sizefact0 <- sizefact <- sizefact_start
    threshold <- threshold_start
    # mat <- matrix(1, nrow(object), 1)
    if (is.null(names(probenum))) names(probenum) <- rownames(object)
    for (iter in seq_len(iterations)) {
        para <- NBth_paraopt(object[features_high, ], probenum[features_high], sizefact, sizefact_BG, threshold, start = start_para)
        # result <- mleprobeNBall(object[,features_high], mat, sizefact0, sizefact,
        #                         matrix(0,1,1), threshold, 0,
        #                         c(rep(0,ncol(mat)), 1, threshold), 0)
        features_NA <- features_high[unique(which(is.na(para), arr.ind = TRUE)[, 2])]
        features_remain <- setdiff(features_high, features_NA)


        for (i in seq_len(length(sizefact))) {
            fun <- NBth_scalenll(object[features_remain, i], probenum[features_remain], t(para[1, features_remain]), t(para[2, features_remain]), sizefact_BG[i], threshold)
            sizefact[i] <- optim(c(sizefact_start[i]), fun, lower = c(lower_sizefact), method = "L-BFGS-B")$par
        }

        if (size_scale == "first") {
            scale_fac <- sizefact[1]
        } else if (size_scale == "sum") {
            scale_fac <- sum(sizefact)
        }

        sizefact <- sizefact / scale_fac


        if (!threshold_fix) {
            fun1 <- NBth_thnll(object[features_remain, ], probenum[features_remain], sizefact, sizefact_BG, scale_fac * t(para[1, features_remain]), t(para[2, features_remain]))

            threshold <- optim(c(threshold_start), fun1, lower = c(lower_threshold), method = "L-BFGS-B")$par
        }

        message(sprintf("Iteration = %s, squared error = %s", iter, sum((sizefact - sizefact0)^2)))


        if (sum((sizefact - sizefact0)^2) < tol) break

        sizefact0 <- sizefact
        
    }
    message("Model converged.")

    rownames(para) <- c("signal", "r")
    
    return(list(
        para0 = NA,
        para = para,
        sizefact = sizefact,
        preci1 = NA,
        conv0 = NA,
        conv = NA,
        Im = NA,
        features_high = features_high,
        features_all = NA,
        threshold = threshold
    ))
}

######### --Load DCC files --##########
count_data <- "1.count_data"
generate_folder(count_data)

#read in translation file.. dsp IDs to slide IDS
dsp2slidetrans <- read.xlsx("./fastqIDCompletetransfile.xlsx", sheet = 1)
dsp2slidetrans$slidename <- str_replace_all(dsp2slidetrans$slidename, " ", "_")

dsp2slidetrans$segment <- str_remove_all(dsp2slidetrans$segment, " ")
dsp2slidetrans$name <- paste(dsp2slidetrans$slidename, dsp2slidetrans$roi, dsp2slidetrans$segment, sep="_")
dsp2slidetrans$name <- str_replace_all(dsp2slidetrans$name, "slide_01_", "slide_01_scan_")
dsp2slidetrans$name <- str_replace_all(dsp2slidetrans$name, "slide_02_", "slide_02_scan_")
dsp2slidetrans$name <- str_remove_all(dsp2slidetrans$name, "\\+")

#ROI
countsROIfilter <- read.xlsx("./221110_CS_GeoMx_RawCountMatrix_JG_LWedit_JGedit2_LWedit2.xlsx", # New counts, raw unfiltered
    sheet = 4,
    # rowNames = T
)
orignames <- colnames(countsROIfilter)
colnames(countsROIfilter) <- gsub("Full.ROI", "FullROI", colnames(countsROIfilter), fixed = T)
colnames(countsROIfilter) <- gsub(".|.", "_", colnames(countsROIfilter), fixed = T)
colnames(countsROIfilter) <- gsub(".", "_", colnames(countsROIfilter), fixed = T)
colnames(countsROIfilter) <- gsub("+", "", colnames(countsROIfilter), fixed = T)

targetrun1 <- read.xlsx("./Copy of M-718 All Data WTA_no filter.xlsx", sheet=1)
targetrun1$Sample_ID <- gsub("Full.ROI", "FullROI", targetrun1$Sample_ID, fixed = T)
targetrun1$Sample_ID <- gsub(".|.", "_", targetrun1$Sample_ID, fixed = T)
targetrun1$Sample_ID <- gsub(".", "_", targetrun1$Sample_ID, fixed = T)
targetrun1$Sample_ID <- gsub("+", "", targetrun1$Sample_ID, fixed = T)
rownames(targetrun1) <- targetrun1$Sample_ID
targetrun1 <- targetrun1[, c("SequencingSaturation", "RawReads", "StitchedReads", 
                                    "AlignedReads", "DeduplicatedReads")]
targetrun1ROI <- targetrun1[rownames(targetrun1) %in% colnames(countsROIfilter),]

targetinforun2 <- read.xlsx("../221110_CS_GeoMx_SegmentSummary.xlsx", sheet = 1)
targetinforun2$Segment.name <- gsub(" ", "", targetinforun2$Segment.name)
rownames(targetinforun2) <- paste(targetinforun2$Scan.name, targetinforun2$ROI.name,
             targetinforun2$Segment.name, sep="_")
targetinforun2 <- targetinforun2[, c("Sequencing.saturation", "Raw.reads", "Stitched.reads", 
                                    "Aligned.reads", "Deduplicated.reads")]
colnames(targetinforun2) <-  c("SequencingSaturation", "RawReads", "StitchedReads", 
                                    "AlignedReads", "DeduplicatedReads")

targetinforun2ROI <- targetinforun2[rownames(targetinforun2) %in% colnames(countsROIfilter), ]
targetinfo <- rbind(targetrun1ROI, targetinforun2ROI)
target_roifilter <- countsROIfilter[1:5, ]
target_roifilter <- as.data.frame(t(target_roifilter))
target_roifilter <- target_roifilter[-1, ]
colnames(target_roifilter) <- target_roifilter[1, ]
target_roifilter <- target_roifilter[-1, ]
target_roifilter$NewSampleName <- paste("ROI", target_roifilter$ROI, target_roifilter$Region, target_roifilter$Animal, sep="_")
target_roifilter$run <- c(rep("Run1", 72), rep("Run2", 45 ))

target_roifilter <- merge(target_roifilter, targetinfo,by.x="row.names", by.y="row.names")
rownames(target_roifilter) <- target_roifilter$NewSampleName
colnames(target_roifilter)[1] <-"orig.id" 
target_roifilter <- as.data.frame(target_roifilter)

countsROIfilter <- countsROIfilter[c(-1, -2, -3, -4, -5), ]
countsROIfilter$X2 <- NULL
negprobesROI <- countsROIfilter[countsROIfilter$TargetName == "NegProbe-WTX",]
rownames(negprobesROI) <- paste(negprobesROI$TargetName, c(1:length(negprobesROI$TargetName)), sep = "_")
negprobesROI$TargetName <- NULL

countsROIfilter <- countsROIfilter[!(countsROIfilter$TargetName == "NegProbe-WTX"), ]
rownames(countsROIfilter) <- countsROIfilter$TargetName
countsROIfilter$TargetName <- NULL

if (all.equal(as.character(target_roifilter$orig.id), colnames(countsROIfilter))==TRUE) {
    message("Printing files...")
    colnames(countsROIfilter) <- target_roifilter$NewSampleName
    target_roifilter$orig_slide_ID <- rownames(target_roifilter)
    colnames(negprobesROI) <- target_roifilter$NewSampleName
    rownames(target_roifilter) <- target_roifilter$NewSampleName
    write.csv(target_roifilter, file.path(count_data, "targetROIfilter.csv"))
    write.csv(countsROIfilter, file.path(count_data, "countsROIfilter.csv"))
} else {
    countsROIfilter<- countsROIfilter[, target_roifilter$orig.id]
    if (all.equal(as.character(target_roifilter$orig.id), colnames(countsROIfilter))==TRUE) {
        message("Printing files attempt 2 ...")
        colnames(countsROIfilter) <- target_roifilter$NewSampleName
        target_roifilter$orig_slide_ID <- rownames(target_roifilter)
        colnames(negprobesROI) <- target_roifilter$NewSampleName
        rownames(target_roifilter) <- target_roifilter$NewSampleName
        write.csv(target_roifilter, file.path(count_data, "targetROIfilter.csv"))
        write.csv(countsROIfilter, file.path(count_data, "countsROIfilter.csv"))
    } else { 
        stop("TARGET FILE AND COUNT FILE FOR ROI DATA IS NOT IN SAME ORDER")
    }
}

dsp2slidetranssub <- dsp2slidetrans[dsp2slidetrans$name %in% target_roifilter$orig.id, ]
target_roifilterfull <- merge(target_roifilter, dsp2slidetranssub, by.x="orig.id", by.y="name")
write.csv(target_roifilterfull, file.path(count_data, "targetROIfilterfull.csv"))

#AOI
countsAOIfilter <- read.xlsx("../221110_CS_GeoMx_RawCountMatrix_JG_LWedit_JGedit2_LWedit2.xlsx", # New counts, raw unfiltered
    sheet = 5,
    # rowNames = T
)
colnames(countsAOIfilter) <- gsub(".|.", "_", colnames(countsAOIfilter), fixed = T)
colnames(countsAOIfilter) <- gsub(".", "_", colnames(countsAOIfilter), fixed = T)
colnames(countsAOIfilter) <- gsub("+", "", colnames(countsAOIfilter), fixed = T)

target_AOIfilter <- countsAOIfilter[1:5, ]
target_AOIfilter <- as.data.frame(t(target_AOIfilter))
target_AOIfilter <- target_AOIfilter[-1, ]
colnames(target_AOIfilter) <- target_AOIfilter[1, ]
target_AOIfilter <- target_AOIfilter[-1, ]
target_AOIfilter$NewSampleName <- paste("AOI", target_AOIfilter$ROI, target_AOIfilter$Region, target_AOIfilter$Animal, sep = "_")
target_AOIfilter$run <- c(rep("Run1", 23), rep("Run2", 1))

targetrun1AOI <- targetrun1[rownames(targetrun1) %in% colnames(countsAOIfilter),]
targetinforun2AOI <- targetinforun2[rownames(targetinforun2) %in% colnames(countsAOIfilter),]
targetinfoaoi <- rbind(targetrun1AOI,targetinforun2AOI )
target_AOIfilter <- merge(target_AOIfilter, targetinfoaoi,by.x="row.names", by.y="row.names")
rownames(target_AOIfilter) <- target_AOIfilter$NewSampleName
colnames(target_AOIfilter)[1] <-"orig.id" 
slidenames <- str_remove_all(target_AOIfilter$orig.id, "_\\d+_GFAP")
slidenames <- str_remove_all(slidenames, "_\\d+_NeuN")
slidenames <- str_remove_all(slidenames, "_\\d+_Olig2")
slidenames <- str_remove_all(slidenames, "_\\d+_Gfap")
target_AOIfilter$slidename <- slidenames
countsAOIfilter <- countsAOIfilter[c(-1, -2, -3, -4, -5), ]
countsAOIfilter$X2 <- NULL
negprobesAOI <- countsAOIfilter[countsAOIfilter$TargetName == "NegProbe-WTX", ]
rownames(negprobesAOI) <- paste(negprobesAOI$TargetName, c(1:length(negprobesAOI$TargetName)), sep = "_")
negprobesAOI$TargetName <- NULL
countsAOIfilter <- countsAOIfilter[!(countsAOIfilter$TargetName == "NegProbe-WTX"), ]
rownames(countsAOIfilter) <- countsAOIfilter$TargetName
countsAOIfilter$TargetName <- NULL

if (all.equal(as.character(target_AOIfilter$orig.id), colnames(countsAOIfilter)) == TRUE) {
    message("Printing attempting AOI")
    colnames(countsAOIfilter) <- target_AOIfilter$NewSampleName
    target_AOIfilter$orig_slide_ID <- rownames(target_AOIfilter)
    colnames(negprobesAOI) <- target_AOIfilter$NewSampleName
    rownames(target_AOIfilter) <- target_AOIfilter$NewSampleName
    write.csv(target_AOIfilter, file.path(count_data, "targetAOIfilter.csv"))
    write.csv(countsAOIfilter, file.path(count_data, "countsAOIfilter.csv"))
} else {
    countsAOIfilter<- countsAOIfilter[, target_AOIfilter$orig.id]
    if (all.equal(as.character(target_AOIfilter$orig.id), colnames(countsAOIfilter)) == TRUE) {
        message("Printing attempting 2 AOI")
        colnames(countsAOIfilter) <- target_AOIfilter$NewSampleName
        target_AOIfilter$orig_slide_ID <- rownames(target_AOIfilter)
        colnames(negprobesAOI) <- target_AOIfilter$NewSampleName
        rownames(target_AOIfilter) <- target_AOIfilter$NewSampleName
        write.csv(target_AOIfilter, file.path(count_data, "targetAOIfilter.csv"))
        write.csv(countsAOIfilter, file.path(count_data, "countsAOIfilter.csv"))
    } else {
       stop("TARGET FILE AND COUNT FILE FOR AOI DATA IS NOT IN SAME ORDER")
    }
}
target_roifilter$condition <- gsub("\\d+", "", target_roifilter$Animal)
target_AOIfilter$condition <- gsub("\\d+", "", target_AOIfilter$Animal)

dsp2slidetransaoisub <- dsp2slidetrans[dsp2slidetrans$name %in% target_AOIfilter$orig.id, ]
target_aoifilterfull <- merge(target_AOIfilter, dsp2slidetransaoisub, by.x = "orig.id", by.y = "name")
write.csv(target_aoifilterfull, file.path(count_data, "targetAOIfilterfull.csv"))

######### -- Boxplots of distribution of Data --########
countsROIfilter <- as.matrix(data.frame(countsROIfilter))
class(countsROIfilter) <- "numeric"

png(file.path(count_data, "ROIrawcountmatrix.png"), width = 9, height=6, units = "in", res=100)
boxplot(log2(countsROIfilter + 1), ylab = "log2 Expression",
         main = "Raw count matrix", names = target_roifilter$run, cex.axis = .6, las = 2)
dev.off()

countsAOIfilter <- as.matrix(data.frame(countsAOIfilter))
class(countsAOIfilter) <- "numeric"

png(file.path(count_data, "AOIrawcountmatrix.png"), width = 9, height = 6, units = "in", res = 100)
boxplot(log2(countsAOIfilter + 1),
    ylab = "log2 Expression",
    main = "Raw count matrix", names = target_AOIfilter$run, cex.axis = .6, las = 2
)
dev.off()

######### --QC--########

#ROI
qc_roi_results <- "1.qc_ROI"
generate_folder(qc_roi_results)
target_roifilter$Nuclei <- as.numeric(target_roifilter$Nuclei)
QC_histogram(target_roifilter, qc_param = "Nuclei", color_by = "run",
                  results_path = qc_roi_results )
target_roifilter$DeduplicatedReads <- as.numeric(target_roifilter$DeduplicatedReads)
QC_histogram(target_roifilter, qc_param = "DeduplicatedReads", color_by = "run",
                 results_path = qc_roi_results )
target_roifilter$AlignedReads <- as.numeric(target_roifilter$AlignedReads)
QC_histogram(target_roifilter, qc_param = "AlignedReads", color_by = "run",
                 results_path = qc_roi_results )
target_roifilter$AlignedReads <- as.numeric(target_roifilter$AlignedReads)

#% Sequencing saturation ([1-deduplicated reads/aligned reads]%): segments below ~50% require additional sequencing to capture full sample diversity and are not typically analyzed until improved.
QC_histogram(target_roifilter,
    qc_param = "SequencingSaturation", color_by = "run",
    results_path = qc_roi_results
)

# ROI
qc_aoi_results <- "1.qc_AOI"
generate_folder(qc_aoi_results)
target_AOIfilter$Nuclei <- as.numeric(target_AOIfilter$Nuclei)
QC_histogram(target_AOIfilter,
    qc_param = "Nuclei", color_by = "run",
    results_path = qc_aoi_results
)
target_AOIfilter$DeduplicatedReads <- as.numeric(target_AOIfilter$DeduplicatedReads)
QC_histogram(target_AOIfilter,
    qc_param = "DeduplicatedReads", color_by = "run",
    results_path = qc_aoi_results
)
target_AOIfilter$AlignedReads <- as.numeric(target_AOIfilter$AlignedReads)
QC_histogram(target_AOIfilter,
    qc_param = "AlignedReads", color_by = "run",
    results_path = qc_aoi_results
)
target_AOIfilter$AlignedReads <- as.numeric(target_AOIfilter$AlignedReads)

# % Sequencing saturation ([1-deduplicated reads/aligned reads]%): segments below ~50% require additional sequencing to capture full sample diversity and are not typically analyzed until improved.
QC_histogram(target_AOIfilter,
    qc_param = "SequencingSaturation", color_by = "run",
    results_path = qc_aoi_results
)

######## --Probe QC--########

# From this count matrix we do not know if there a multiple probes per gene therefore 
# I am focusing the probe QC on the negative controls 

# first I will calculate the geoMean for each of the neg probes and divdie it 
# by the geometric mean of all the probcounts
negprobesROI <- as.matrix(data.frame(negprobesROI))
class(negprobesROI) <- "numeric"
neggeomeanROIprob <- apply(negprobesROI, MARGIN = 1, FUN = ngeoMean) # apply(negprobesROI, 2, function(x) exp(mean(log(x))))
probecounts <- c()
for (sample in colnames(negprobesROI)) { 
    probecounts <- c(probecounts, negprobesROI[, sample])
}

totalngeMeancounts <- ngeoMean(probecounts)
neggeomeanROIprobfinalmetric <- neggeomeanROIprob/totalngeMeancounts
neggeomeanROIprobfinalmetric_melt <- melt(neggeomeanROIprobfinalmetric)
neggeomeanROIprobfinalmetric_melt$name <- rownames(neggeomeanROIprobfinalmetric_melt)
# this value represents the prob counts all segments is less than 0.1 (anything above this is 1 )
ggplot(neggeomeanROIprobfinalmetric_melt, aes(x=name, y = value)) +
    geom_bar(stat="identity", width=0.5)+ geom_hline(yintercept = 0.1, linetype="dashed", color="red") +
     coord_flip() + theme(axis.text.y = element_text(size=5)) + labs(y="geoMean neg probe / geoMean of count probes", x="")
ggsave(file.path(qc_roi_results, "NegProbeQC.png"))

#AOI
negprobesAOI <- as.matrix(data.frame(negprobesAOI))
class(negprobesAOI) <- "numeric"
negprobesAOIprob <- apply(negprobesAOI, MARGIN = 1, FUN = ngeoMean) # apply(negprobesROI, 2, function(x) exp(mean(log(x))))
probecounts <- c()
for (sample in colnames(countsAOIfilter)) {
    probecounts <- c(probecounts, countsAOIfilter[, sample])
}

totalngeMeancounts <- ngeoMean(probecounts)
neggeomeanAOIprobfinalmetric <- negprobesAOIprob / totalngeMeancounts
neggeomeanAOIprobfinalmetric_melt <- melt(neggeomeanAOIprobfinalmetric)
neggeomeanAOIprobfinalmetric_melt$name <- rownames(neggeomeanAOIprobfinalmetric_melt)
# this value represents the prob counts all segments is less than 0.1 (anything above this is 1 )
ggplot(neggeomeanAOIprobfinalmetric_melt, aes(x = name, y = value)) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
    coord_flip() +
    theme(axis.text.y = element_text(size = 5)) +
    labs(y = "geoMean neg probe / geoMean of count probes", x = "")
ggsave(file.path(qc_aoi_results, "NegProbeQC.png"))

######## --Calulating limit of quantification (LOQ) (background noise) -- ######## 

#The LOQ is calculated based on the distribution of negative control probes and is 
#intended to approximate the quantifiable limit of gene expression per segment. 
#Please note that this process is more stable in larger segments.
cutoff <- 2
minLOQ <- 2

##ROI 

negprobesROI <- as.matrix(data.frame(negprobesROI))
class(negprobesROI) <- "numeric"
neggeomeanROI <- apply(negprobesROI, MARGIN = 2, FUN=ngeoMean) #apply(negprobesROI, 2, function(x) exp(mean(log(x))))
negmeanROI <- apply(negprobesROI, MARGIN = 2, FUN=mean) 
neggeoSDROI <- apply(negprobesROI, MARGIN = 2, FUN = ngeoSD)
if (all.equal(names(neggeomeanROI), rownames(target_roifilter))==TRUE) {
    target_roifilter$negprobgeomean <- neggeomeanROI
    target_roifilter$negprobgeoSD <- neggeoSDROI
    target_roifilter$LOQ <- pmax(minLOQ, neggeomeanROI * neggeoSDROI^cutoff) # pmax is assigning a value of 2 to the LOQ if it is lessthan 2

} else {
    stop("1. Columns need to be in the same order ")
}
write.csv(target_roifilter, file.path(count_data, "targetROIfilterLOQ.csv"))
QC_histogram(target_roifilter,
    qc_param = "LOQ", color_by = "run",
    results_path = qc_roi_results
)

## AOI
negprobesAOI <- as.matrix(data.frame(negprobesAOI))
class(negprobesAOI) <- "numeric"
neggeomeanAOI <- apply(negprobesAOI, MARGIN = 2, FUN = ngeoMean) # apply(negprobesROI, 2, function(x) exp(mean(log(x))))
neggeoSDAOI <- apply(negprobesAOI, MARGIN = 2, FUN = ngeoSD)
if (all.equal(names(neggeomeanAOI), rownames(target_AOIfilter)) == TRUE) {
    target_AOIfilter$negprobgeomean <- neggeomeanAOI
    target_AOIfilter$negprobgeoSD <- neggeoSDAOI
    target_AOIfilter$LOQ <- pmax(minLOQ, neggeoSDAOI * neggeomeanAOI^cutoff) # pmax is assigning a value of 2 to the LOQ if it is lessthan 2
} else {
    stop("1. Columns need to be in the same order ")
}
write.csv(target_AOIfilter, file.path(count_data, "targetAOIfilterLOQ.csv"))
QC_histogram(target_AOIfilter,
    qc_param = "LOQ", color_by = "run",
    results_path = qc_aoi_results
)

### Number of genes above the LOQ
LOQ_Mat <- c()

LOQ_Mat <- t(apply(countsROIfilter,
    MARGIN = 1,
    FUN = function(x) {
        x > target_roifilter$LOQ
    }
))
LOQ_Mat <- LOQ_Mat[rownames(countsROIfilter), ]

target_roifilter$GenesDetected <- colSums(LOQ_Mat, na.rm = TRUE)
target_roifilter$GeneDetectionRate <- target_roifilter$GenesDetected / nrow(countsROIfilter)
target_roifilter$DetectionThreshold <-
    cut(target_roifilter$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%")
    )
ggplot(
    target_roifilter,
    aes(x = DetectionThreshold)
) +
    geom_bar(aes(fill = run)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() + scale_fill_brewer(palette = "Dark2") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
        x = "Gene Detection Rate",
        y = "Segments, #",
        fill = "run"
    )
ggsave(file.path(qc_roi_results, "DetectionThreshold.png"), width=4, height = 5, dpi=300)

ggplot(
    target_roifilter,
    aes(x = DetectionThreshold)
) +
    geom_bar(aes(fill = Region)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    scale_fill_manual(values=c('red','darkblue',"purple"))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
        x = "Gene Detection Rate",
        y = "Segments, #",
        fill = "Region"
    )
ggsave(file.path(qc_roi_results, "DetectionThresholdRegion.png"), width = 4, height = 5, dpi = 300)

ggplot(
    target_roifilter,
    aes(x = DetectionThreshold)
) +
    geom_bar(aes(fill = Animal)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    # scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
        x = "Gene Detection Rate",
        y = "Segments, #",
        fill = "Region"
    )
ggsave(file.path(qc_roi_results, "DetectionThresholdAnimal.png"), width = 4, height = 5, dpi = 300)

ggplot(
    target_roifilter,
    aes(x = Animal, y = LOQ, fill=Region)
) +
    geom_bar(stat = "identity",  position=position_dodge()) + scale_fill_manual(values=c('red','darkblue',"purple"))+
    # geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw()
    # scale_fill_brewer(palette = "Paired") +
    # scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    # labs(
    #     x = "Gene Detection Rate",
    #     y = "Segments, #",
    #     fill = "Region"
    # )
ggsave(file.path(qc_roi_results, "DetectionThresholdAnimalv2.png"), width = 5, height = 5, dpi = 300)

ggplot(
    target_roifilter,
    aes(x = DetectionThreshold)
) +
    geom_bar(aes(fill = condition)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    scale_fill_manual(values=c('#999999','black'))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
        x = "Gene Detection Rate",
        y = "Segments, #",
        fill = "Condition"
    )
ggsave(file.path(qc_roi_results, "DetectionThresholdCondition.png"), width = 4, height = 5, dpi = 300)

detecttable <- kable(table(
    target_roifilter$DetectionThreshold,
    target_roifilter$run
))

# filter based on 1% of the genes being above the LOQ
thresholdgd <- 0.01
target_roifilter10 <- target_roifilter[target_roifilter$GeneDetectionRate >= thresholdgd, ]
message("STATUS: Number of samples above the ",  thresholdgd, " threshold ", nrow(target_roifilter10)) 
#figure showing distribution zika animals and brain region
countsROIfilter10 <- countsROIfilter[, rownames(target_roifilter10)]
LOQ_Mat10 <- LOQ_Mat[, colnames(countsROIfilter10)]
Genes_LOQ <- rowSums(LOQ_Mat, na.rm = TRUE)
Genes_LOQ <- as.data.frame(Genes_LOQ)
colnames(Genes_LOQ) <- c("DetectedSegments")
Genes_LOQ$DetectionRate <- Genes_LOQ$DetectedSegments / nrow(target_roifilter10)
#Looking at genes detected across the datasets
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
    unlist(lapply(
        c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
        function(x) {
            sum(Genes_LOQ$DetectionRate >= x)
        }
    ))
plot_detect$Rate <- plot_detect$Number / nrow(countsROIfilter)
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
        vjust = 1.6, color = "black", size = 4
    ) +
    scale_fill_gradient2(
        low = "orange2", mid = "lightblue",
        high = "dodgerblue3", midpoint = 0.65,
        limits = c(0, 1),
        labels = scales::percent
    ) +
    theme_bw() +
    scale_y_continuous(
        labels = scales::percent, limits = c(0, 1),
        expand = expansion(mult = c(0, 0))
    ) +
    labs(
        x = "% of Segments",
        y = "Genes Detected, % of Panel > LOQ"
    )
ggsave(file.path(qc_roi_results, "GenesDetected.png"))

#Filter genes 
genesfiltered <- Genes_LOQ[Genes_LOQ$DetectionRate > 0.15, ]
countsROIfilter10 <- countsROIfilter10[rownames(genesfiltered),]
message("STATUS: Number of genes kept ", nrow(genesfiltered), " out of ", nrow(Genes_LOQ))

### AOI Number of genes above the LOQ
LOQ_Mat <- c()

LOQ_Mat <- t(apply(countsAOIfilter,
    MARGIN = 1,
    FUN = function(x) {
        x > target_AOIfilter$LOQ
    }
))
LOQ_Mat <- LOQ_Mat[rownames(countsAOIfilter), ]

target_AOIfilter$GenesDetected <- colSums(LOQ_Mat, na.rm = TRUE)
target_AOIfilter$GeneDetectionRate <- target_AOIfilter$GenesDetected / nrow(countsAOIfilter)
target_AOIfilter$DetectionThreshold <-
    cut(target_AOIfilter$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%")
    )
ggplot(
    target_AOIfilter,
    aes(x = DetectionThreshold)
) +
    geom_bar(aes(fill = run)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() + scale_fill_brewer(palette = "Dark2") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
        x = "Gene Detection Rate",
        y = "Segments, #",
        fill = "run"
    )
ggsave(file.path(qc_aoi_results, "DetectionThreshold.png"), width=4, height = 5, dpi=300)

ggplot(
    target_AOIfilter,
    aes(x = DetectionThreshold)
) +
    geom_bar(aes(fill = Region)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    scale_fill_manual(values=c('red','darkblue',"purple"))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
        x = "Gene Detection Rate",
        y = "Segments, #",
        fill = "Region"
    )
ggsave(file.path(qc_aoi_results, "DetectionThresholdRegion.png"), width = 4, height = 5, dpi = 300)

ggplot(
    target_AOIfilter,
    aes(x = DetectionThreshold)
) +
    geom_bar(aes(fill = Animal)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    # scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
        x = "Gene Detection Rate",
        y = "Segments, #",
        fill = "Region"
    )
ggsave(file.path(qc_aoi_results, "DetectionThresholdAnimal.png"), width = 4, height = 5, dpi = 300)

ggplot(
    target_AOIfilter,
    aes(x = Animal, y = LOQ, fill=Region)
) +
    geom_bar(stat = "identity",  position=position_dodge()) + scale_fill_manual(values=c('red','darkblue',"purple"))+
    # geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw()
    # scale_fill_brewer(palette = "Paired") +
    # scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    # labs(
    #     x = "Gene Detection Rate",
    #     y = "Segments, #",
    #     fill = "Region"
    # )
ggsave(file.path(qc_aoi_results, "DetectionThresholdAnimalv2.png"), width = 5, height = 5, dpi = 300)

ggplot(
    target_AOIfilter,
    aes(x = DetectionThreshold)
) +
    geom_bar(aes(fill = condition)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    scale_fill_manual(values=c('#999999','black'))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
        x = "Gene Detection Rate",
        y = "Segments, #",
        fill = "condition"
    )
ggsave(file.path(qc_aoi_results, "DetectionThresholdCondition.png"), width = 4, height = 5, dpi = 300)

detecttable <- kable(table(
    target_AOIfilter$DetectionThreshold,
    target_AOIfilter$run
))
write.csv(target_AOIfilter, file.path(count_data, "targetAOIgenedetection.csv"))
# AOI filter based on 1% of the genes being above the LOQ
thresholdgd <- 0.01
target_AOIfilter10 <- target_AOIfilter[target_AOIfilter$GeneDetectionRate >= thresholdgd, ]
message("STATUS: Number of samples above the AOI",  thresholdgd, " threshold ", nrow(target_AOIfilter10)) 
#figure showing distribution zika animals and brain region
countsAOIfilter10 <- countsAOIfilter[, rownames(target_AOIfilter10)]
LOQ_Mat10 <- LOQ_Mat[, colnames(countsAOIfilter10)]
Genes_LOQ <- rowSums(LOQ_Mat, na.rm = TRUE)
Genes_LOQ <- as.data.frame(Genes_LOQ)
colnames(Genes_LOQ) <- c("DetectedSegments")
Genes_LOQ$DetectionRate <- Genes_LOQ$DetectedSegments / nrow(target_AOIfilter10)
#Looking at genes detected across the datasets
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
    unlist(lapply(
        c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
        function(x) {
            sum(Genes_LOQ$DetectionRate >= x)
        }
    ))
plot_detect$Rate <- plot_detect$Number / nrow(countsAOIfilter)
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
        vjust = 1.6, color = "black", size = 4
    ) +
    scale_fill_gradient2(
        low = "orange2", mid = "lightblue",
        high = "dodgerblue3", midpoint = 0.65,
        limits = c(0, 1),
        labels = scales::percent
    ) +
    theme_bw() +
    scale_y_continuous(
        labels = scales::percent, limits = c(0, 1),
        expand = expansion(mult = c(0, 0))
    ) +
    labs(
        x = "% of Segments",
        y = "Genes Detected, % of Panel > LOQ"
    )
ggsave(file.path(qc_aoi_results, "GenesDetected.png"))

#Filter genes 
genesfiltered <- Genes_LOQ[Genes_LOQ$DetectionRate > 0.25, ]
genesfilteredAOI25 <- Genes_LOQ[Genes_LOQ$DetectionRate > 0.25, ]
genesfilteredAOI15 <- Genes_LOQ[Genes_LOQ$DetectionRate > 0.15, ]

countsAOIfilter10 <- countsAOIfilter10[rownames(genesfiltered),]
message("STATUS: AOI Number of genes kept ", nrow(genesfiltered), " out of ", nrow(Genes_LOQ))

######## --Normalization ROI --########
norm_results <- "1.norm_roi_results"
generate_folder(norm_results)

neggeomeanROI10 <- neggeomeanROI[colnames(countsROIfilter10)]
Stat_data <-
    data.frame(
        row.names = colnames(countsROIfilter10),
        Segment = colnames(countsROIfilter10),
        Annotation = target_roifilter10$Region,
        Q3 = unlist(apply(countsROIfilter10, 2,
            quantile, 0.75,
            na.rm = TRUE
        )), # identifies # genes above the the 0.75 th quantile  
        NegProbe = neggeomeanROI10
    )

Stat_data_m <- melt(Stat_data,
    measure.vars = c("Q3", "NegProbe"),
    variable.name = "Statistic", value.name = "Value"
)
pre_norm_figure(Stat_data, Stat_data_m, norm_results)

# based on this figure I am  going to remove the one GM 
# sample that has a higher negprobe geomeaan
Stat_data$q3vsNegprobe <- Stat_data$Q3 / Stat_data$NegProbe
removesamples <- rownames(Stat_data[Stat_data$q3vsNegprobe < 1, ])

message("STATUS: removing samples as they had a low q3 to negprobe geomean ratio ", length(removesamples))

Stat_data_removedoutliers <- Stat_data[!Stat_data$Segment %in% removesamples, ]
Stat_data_removedoutliers_m <- melt(Stat_data_removedoutliers,
    measure.vars = c("Q3", "NegProbe"),
    variable.name = "Statistic", value.name = "Value"
)
pre_norm_figure(Stat_data_removedoutliers, Stat_data_removedoutliers_m, norm_results, figname = "NormalizationPlot_removedoutliers.png")

target_roifilter10final <- target_roifilter10[! rownames(target_roifilter10) %in% removesamples, ]
countsROIfilter10final <- countsROIfilter10[, rownames(target_roifilter10final)]

neggeomeanROI10final <- negprobesROI[,  rownames(target_roifilter10final)]

message("STATUS: Number of samples remaining for normalization and analysis is ", nrow(target_roifilter10final))
write.csv(target_roifilter10final, file.path(count_data, "targetROIfilterLOQfinal.csv"))
write.csv(countsROIfilter10final, file.path(count_data, "countsROIfilterfinal.csv"))
if (all.equal(colnames(countsROIfilter10final), rownames(Stat_data_removedoutliers))) {
    norm_matrix <- Q3_norm(Stat_data_removedoutliers, countsROIfilter10final)
} else { 
    stop("WARNING: Stats_data and countmatrix are not in the same order so not normalizing ")
}
write.csv(norm_matrix, file.path(norm_results, "1.norm_matrix.csv"))

png(file.path(norm_results, "filtered_Q3norm_boxplot.png"))
boxplot(norm_matrix,
    col = "#2CA02C", , main = "Q3 Norm Counts", cex.axis = .6, las = 2,
     xlab = "Segment", names=target_roifilter10final$run,
    frame=FALSE, ylab = "Q3 Norm Counts"
)
dev.off()

norm_matrix_log <- log2(norm_matrix)
write.csv(norm_matrix_log, file.path(norm_results, "1.lognorm_matrix.csv"))
png(file.path(norm_results, "filtered_logQ3norm_boxplot.png"))
boxplot(norm_matrix_log,
    col = "#2CA02C", , main = "Q3 Norm (log2) Counts", cex.axis = .6, las = 2,
    log = "y", xlab = "Segment", names = target_roifilter10final$run,
    frame = FALSE, ylab = "Q3 Norm (log2) Counts"
)
dev.off()

# ######## --Normalization AOI --########
# norm_results <- "1.norm_AOI_results"
# generate_folder(norm_results)

# neggeomeanAOI10 <- neggeomeanAOI[colnames(countsAOIfilter10)]
# Stat_data <-
#     data.frame(
#         row.names = colnames(countsAOIfilter10),
#         Segment = colnames(countsAOIfilter10),
#         Annotation = target_AOIfilter10$Region,
#         Q3 = unlist(apply(countsAOIfilter10, 2,
#             quantile, 0.75,
#             na.rm = TRUE
#         )), # identifies # genes above the the 0.75 th quantile  
#         NegProbe = neggeomeanAOI10
#     )

# Stat_data_m <- melt(Stat_data,
#     measure.vars = c("Q3", "NegProbe"),
#     variable.name = "Statistic", value.name = "Value"
# )
# pre_norm_figure(Stat_data, Stat_data_m, norm_results)

# # based on this figure I am  going to remove the one GM 
# # sample that has a higher negprobe geomeaan
# Stat_data$q3vsNegprobe <- Stat_data$Q3 / Stat_data$NegProbe
# removesamples <- rownames(Stat_data[Stat_data$q3vsNegprobe < 1, ])

# message("STATUS: removing samples as they had a low q3 to negprobe geomean ratio ", length(removesamples))

# Stat_data_removedoutliers <- Stat_data[!Stat_data$Segment %in% removesamples, ]
# Stat_data_removedoutliers_m <- melt(Stat_data_removedoutliers,
#     measure.vars = c("Q3", "NegProbe"),
#     variable.name = "Statistic", value.name = "Value"
# )
# pre_norm_figure(Stat_data_removedoutliers, Stat_data_removedoutliers_m, norm_results, figname = "NormalizationPlot_removedoutliers.png")

# target_AOIfilter10final <- target_AOIfilter10[! rownames(target_AOIfilter10) %in% removesamples, ]
# countsAOIfilter10final <- countsAOIfilter10[, rownames(target_AOIfilter10final)]

# neggeomeanAOI10final <- negprobesAOI[,  rownames(target_AOIfilter10final)]

# message("STATUS: Number of samples remaining for normalization and analysis is ", nrow(target_AOIfilter10final))
# write.csv(target_AOIfilter10final, file.path(count_data, "targetAOIfilterLOQfinal.csv"))
# write.csv(countsAOIfilter10final, file.path(count_data, "countsAOIfilterfinal.csv"))
# if (all.equal(colnames(countsAOIfilter10final), rownames(Stat_data_removedoutliers))) {
#     norm_matrixAOI <- Q3_norm(Stat_data_removedoutliers, countsAOIfilter10final)
# } else { 
#     stop("WARNING: Stats_data and countmatrix are not in the same order so not normalizing ")
# }
# write.csv(norm_matrixAOI, file.path(norm_results, "1.norm_matrix.csv"))

# png(file.path(norm_results, "filtered_Q3norm_boxplot.png"))
# boxplot(norm_matrixAOI,
#     col = "#2CA02C", , main = "Q3 Norm Counts", cex.axis = .6, las = 2,
#     log = "y", xlab = "Segment", names=target_AOIfilter10final$run,
#      frame=FALSE, ylab = "Q3 Norm Counts"
# )
# dev.off()

# norm_matrixAOI_log <- log2(norm_matrixAOI)
# write.csv(norm_matrixAOI_log, file.path(norm_results, "1.lognorm_matrix.csv"))
# png(file.path(norm_results, "filtered_logQ3norm_boxplot.png"))
# boxplot(norm_matrixAOI_log,
#     col = "#2CA02C", , main = "Q3 Norm (log2) Counts", cex.axis = .6, las = 2,
#     xlab = "Segment", names = target_AOIfilter10final$run,
#     frame = FALSE, ylab = "Q3 Norm (log2) Counts"
# )
# dev.off()

# ######## --Normalization AOI without QC --########
# norm_results <- "1.norm_AOI_results_no_QC"
# generate_folder(norm_results)

# ## No filtering ##
# neggeomeanAOI <- neggeomeanAOI[colnames(countsAOIfilter)]

# Stat_data <-
#     data.frame(
#         row.names = colnames(countsAOIfilter),
#         Segment = colnames(countsAOIfilter),
#         Annotation = target_AOIfilter$Region,
#         Q3 = unlist(apply(countsAOIfilter, 2,
#             quantile, 0.75,
#             na.rm = TRUE
#         )), # identifies # genes above the the 0.75 th quantile  
#         NegProbe = neggeomeanAOI
#     )

# Stat_data_m <- melt(Stat_data,
#     measure.vars = c("Q3", "NegProbe"),
#     variable.name = "Statistic", value.name = "Value"
# )
# pre_norm_figure(Stat_data, Stat_data_m, norm_results)

# norm_matrixAOI_no_QC <- Q3_norm(Stat_data, countsAOIfilter)
# norm_matrixAOI_no_QC_log2 <- log2(norm_matrixAOI_no_QC)

# ## 25% filter##
# norm_results <- "1.norm_AOI_results_no_QC_25"
# generate_folder(norm_results)

# countsAOIfilter25 <- countsAOIfilter[rownames(genesfilteredAOI25), ]
# write.csv(countsAOIfilter25, file.path(norm_results, "countsAOIfilter25.csv"))
# neggeomeanAOI <- neggeomeanAOI[colnames(countsAOIfilter25)]

# Stat_data <-
#     data.frame(
#         row.names = colnames(countsAOIfilter25),
#         Segment = colnames(countsAOIfilter25),
#         Annotation = target_AOIfilter$Region,
#         Q3 = unlist(apply(countsAOIfilter25, 2,
#             quantile, 0.75,
#             na.rm = TRUE
#         )), # identifies # genes above the the 0.75 th quantile  
#         NegProbe = neggeomeanAOI
#     )

# Stat_data_m <- melt(Stat_data,
#     measure.vars = c("Q3", "NegProbe"),
#     variable.name = "Statistic", value.name = "Value"
# )
# pre_norm_figure(Stat_data, Stat_data_m, norm_results)

# norm_matrixAOI_no_QC_25 <- Q3_norm(Stat_data, countsAOIfilter25)
# norm_matrixAOI_no_QC_25_log2 <- log2(norm_matrixAOI_no_QC_25)

# ## 15% filter##
# norm_results <- "1.norm_AOI_results_no_QC_15"
# generate_folder(norm_results)
# countsAOIfilter15 <- countsAOIfilter[rownames(genesfilteredAOI15), ]
# neggeomeanAOI <- neggeomeanAOI[colnames(genesfilteredAOI15)]

# Stat_data <-
#     data.frame(
#         row.names = colnames(countsAOIfilter15),
#         Segment = colnames(countsAOIfilter15),
#         Annotation = target_AOIfilter$Region,
#         Q3 = unlist(apply(countsAOIfilter15, 2,
#             quantile, 0.75,
#             na.rm = TRUE
#         )), # identifies # genes above the the 0.75 th quantile  
#         NegProbe = neggeomeanAOI
#     )

# Stat_data_m <- melt(Stat_data,
#     measure.vars = c("Q3", "NegProbe"),
#     variable.name = "Statistic", value.name = "Value"
# )
# pre_norm_figure(Stat_data, Stat_data_m, norm_results)

# norm_matrixAOI_no_QC_15 <- Q3_norm(Stat_data, countsAOIfilter15)
# norm_matrixAOI_no_QC_15_log2 <- log2(norm_matrixAOI_no_QC_15)

######## --Feature Reduction --########
feature_results <- "1.feature_red"
generate_folder(feature_results)

pca <- pca_fun(
    norm_matrix, target_roifilter10final,
    feature_results, "RegionRun.png",
    c("Region", "run"), 300, 3
)

pca <- pca_fun(
    norm_matrix, target_roifilter10final,
    feature_results, "Regioncondition.png",
    c("Region", "condition"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrix, target_roifilter10final,
    feature_results, "RegionAnimal.png",
    c("Region", "Animal"), 300, 3, pca=pca
)

targettmp <-  target_roifilter10final[target_roifilter10final$condition=="Z", ]
norm_matrixtmp <-  norm_matrix[, rownames(targettmp)]
pca <- pca_fun(
    norm_matrixtmp, targettmp,
    feature_results, "RegionAnimalZika.png",
    c("Region", "Animal"), 300, 3
)
pca <- pca_fun(
    norm_matrixtmp, targettmp,
    feature_results, "RegionAnimalZika.pdf",
    c("Region", "Animal"), 300, 3, pca=pca
)
pca <- pca_fun(
    log2(norm_matrixtmp), targettmp,
    feature_results, "LogRegionAnimalZika.png",
    c("Region", "Animal"), 300, 3
)
pca <- pca_fun(
    log2(norm_matrixtmp), targettmp,
    feature_results, "LogRegionAnimalZika.pdf",
    c("Region", "Animal"), 300, 3, pca=pca
)

DEgenes <- read.csv("./1.de_results_full/sigLFCs.csv", row.names=1)
targettmp <-  target_roifilter10final[target_roifilter10final$condition=="Z", ]
norm_matrixtmp<-  norm_matrix[rownames(DEgenes), rownames(targettmp)]
pca <- pca_fun(
    norm_matrixtmp, targettmp,
    feature_results, "RegionAnimalZikaDE.png",
    c("Region", "Animal"), 300, 3
)
pca <- pca_fun(
    norm_matrixtmp, targettmp,
    feature_results, "RegionAnimalZikaDE.pdf",
    c("Region", "Animal"), 300, 3, pca=pca
)
pca <- pca_fun(
    log2(norm_matrixtmp), targettmp,
    feature_results, "LogRegionAnimalZikaDE.png",
    c("Region", "Animal"), 300, 3
)
pca <- pca_fun(
    log2(norm_matrixtmp), targettmp,
    feature_results, "LogRegionAnimalZikaDE.pdf",
    c("Region", "Animal"), 300, 3,  pca=pca
)

umap <- umap_fun(
   norm_matrix, target_roifilter10final,
    feature_results, "RegionRun.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
   norm_matrix, target_roifilter10final,
    feature_results, "Regioncondition.png",
    c("Region", "condition"), 300, 3, UMAP=umap
)

umap <- umap_fun(
   norm_matrix, target_roifilter10final,
    feature_results, "RegionAnimal.png",
    c("Region", "Animal"), 300, 3, UMAP=umap
)

pca <- pca_fun(
    norm_matrix_log, target_roifilter10final,
    feature_results, "logRegionRun.png",
    c("Region", "run"), 300, 3
)

pca <- pca_fun(
    norm_matrix_log, target_roifilter10final,
    feature_results, "logRegioncondition.png",
    c("Region", "condition"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrix_log, target_roifilter10final,
    feature_results, "logRegionAnimal.png",
    c("Region", "Animal"), 300, 3, pca=pca
)

umap <- umap_fun(
   norm_matrix_log, target_roifilter10final,
    feature_results, "logRegionRun.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
   norm_matrix_log, target_roifilter10final,
    feature_results, "logRegioncondition.png",
    c("Region", "condition"), 300, 3, UMAP=umap
)

umap <- umap_fun(
   norm_matrix_log, target_roifilter10final,
    feature_results, "logRegionAnimal.png",
    c("Region", "Animal"), 300, 3, UMAP=umap
)

######## --Feature Reduction AOI --########
feature_results <- "1.feature_red_AOI"
generate_folder(feature_results)

pca <- pca_fun(
    norm_matrixAOI, target_AOIfilter10final,
    feature_results, "RegionRun.png",
    c("Region", "run"), 300, 3
)

pca <- pca_fun(
    norm_matrixAOI, target_AOIfilter10final,
    feature_results, "Regioncondition.png",
    c("Region", "condition"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrixAOI, target_AOIfilter10final,
    feature_results, "RegionAnimal.png",
    c("Region", "Animal"), 300, 3, pca=pca
)

umap <- umap_fun(
   norm_matrixAOI, target_AOIfilter10final,
    feature_results, "RegionRun.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
   norm_matrixAOI, target_AOIfilter10final,
    feature_results, "Regioncondition.png",
    c("Region", "condition"), 300, 3, UMAP=umap
)

umap <- umap_fun(
   norm_matrixAOI, target_AOIfilter10final,
    feature_results, "RegionAnimal.png",
    c("Region", "Animal"), 300, 3, UMAP=umap
)

pca <- pca_fun(
    norm_matrixAOI_log, target_AOIfilter10final,
    feature_results, "logRegionRun.png",
    c("Region", "run"), 300, 3
)

pca <- pca_fun(
    norm_matrixAOI_log, target_AOIfilter10final,
    feature_results, "logRegioncondition.png",
    c("Region", "condition"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrixAOI_log, target_AOIfilter10final,
    feature_results, "logRegionAnimal.png",
    c("Region", "Animal"), 300, 3, pca=pca
)

umap <- umap_fun(
   norm_matrixAOI_log, target_AOIfilter10final,
    feature_results, "logRegionRun.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
   norm_matrixAOI_log, target_AOIfilter10final,
    feature_results, "logRegioncondition.png",
    c("Region", "condition"), 300, 3, UMAP=umap
)

umap <- umap_fun(
   norm_matrixAOI_log, target_AOIfilter10final,
    feature_results, "logRegionAnimal.png",
    c("Region", "Animal"), 300, 3, UMAP=umap
)

######## --Feature Reduction AOI no QC--########
feature_results <- "1.feature_red_AOI_no_QC"
generate_folder(feature_results)

pca <- pca_fun(
    norm_matrixAOI_no_QC, target_AOIfilter,
    feature_results, "RegionRun.png",
    c("Region", "run"), 300, 3
)

pca <- pca_fun(
    norm_matrixAOI_no_QC, target_AOIfilter,
    feature_results, "Regioncondition.png",
    c("Region", "condition"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrixAOI_no_QC, target_AOIfilter,
    feature_results, "RegionAnimal.png",
    c("Region", "Animal"), 300, 3, pca=pca
)

umap <- umap_fun(
   norm_matrixAOI_no_QC, target_AOIfilter,
    feature_results, "RegionRun.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
   norm_matrixAOI_no_QC, target_AOIfilter,
    feature_results, "Regioncondition.png",
    c("Region", "condition"), 300, 3, UMAP=umap
)

umap <- umap_fun(
   norm_matrixAOI_no_QC, target_AOIfilter,
    feature_results, "RegionAnimal.png",
    c("Region", "Animal"), 300, 3, UMAP=umap
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_log2, target_AOIfilter,
    feature_results, "logRegionRun.png",
    c("Region", "run"), 300, 3
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_log2, target_AOIfilter,
    feature_results, "logRegioncondition.png",
    c("Region", "condition"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_log2, target_AOIfilter,
    feature_results, "logRegionAnimal.png",
    c("Region", "Animal"), 300, 3, pca=pca
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_log2, target_AOIfilter,
    feature_results, "logRegionRun.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_log2, target_AOIfilter,
    feature_results, "logRegioncondition.png",
    c("Region", "condition"), 300, 3, UMAP=umap
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_log2, target_AOIfilter,
    feature_results, "logRegionAnimal.png",
    c("Region", "Animal"), 300, 3, UMAP=umap
)

######## --Feature Reduction AOI no QC 25--########
outliers <- setdiff(rownames(target_AOIfilter), rownames(target_AOIfilter10final))
outlierstatus <- c()
for (i in rownames(target_AOIfilter)) {
    if (i %in% outliers) {
        outlierstatus <- c(outlierstatus, "Outlier")
    } else { 
        outlierstatus <- c(outlierstatus, "No Outlier")
    }
}
target_AOIfilter$Outliers <- outlierstatus
feature_results <- "1.feature_red_AOI_no_QC_25"
generate_folder(feature_results)

pca <- pca_fun(
    norm_matrixAOI_no_QC_25, target_AOIfilter,
    feature_results, "RegionRun.png",
    c("Region", "run"), 300, 3
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_25, target_AOIfilter,
    feature_results, "Regioncondition.png",
    c("Region", "condition"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_25, target_AOIfilter,
    feature_results, "RegionAnimal.png",
    c("Region", "Animal"), 300, 3, pca=pca
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_25, target_AOIfilter,
    feature_results, "RegionRun.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_25, target_AOIfilter,
    feature_results, "Regioncondition.png",
    c("Region", "condition"), 300, 3, UMAP=umap
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_25, target_AOIfilter,
    feature_results, "RegionAnimal.png",
    c("Region", "Animal"), 300, 3, UMAP=umap
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_25_log2, target_AOIfilter,
    feature_results, "logRegionRun.png",
    c("Region", "run"), 300, 3
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_25_log2, target_AOIfilter,
    feature_results, "logRegionOutliers.png",
    c("Region", "Outliers"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_25_log2, target_AOIfilter,
    feature_results, "logRegioncondition.png",
    c("Region", "condition"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_25_log2, target_AOIfilter,
    feature_results, "logRegionAnimal.png",
    c("Region", "Animal"), 300, 3, pca=pca
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_25_log2, target_AOIfilter,
    feature_results, "logRegionRun.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
    norm_matrixAOI_no_QC_25_log2, target_AOIfilter,
    feature_results, "logRegionOutliers.png",
    c("Region", "Outliers"), 300, 3,
    UMAP = umap
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_25_log2, target_AOIfilter,
    feature_results, "logRegioncondition.png",
    c("Region", "condition"), 300, 3, UMAP=umap
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_25_log2, target_AOIfilter,
    feature_results, "logRegionAnimal.png",
    c("Region", "Animal"), 300, 3, UMAP=umap
)

######## --Feature Reduction AOI no QC 15--########
feature_results <- "1.feature_red_AOI_no_QC_15"
generate_folder(feature_results)

pca <- pca_fun(
    norm_matrixAOI_no_QC_15, target_AOIfilter,
    feature_results, "RegionRun.png",
    c("Region", "run"), 300, 3
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_15, target_AOIfilter,
    feature_results, "Regioncondition.png",
    c("Region", "condition"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_15, target_AOIfilter,
    feature_results, "RegionAnimal.png",
    c("Region", "Animal"), 300, 3, pca=pca
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_15, target_AOIfilter,
    feature_results, "RegionRun.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_15, target_AOIfilter,
    feature_results, "Regioncondition.png",
    c("Region", "condition"), 300, 3, UMAP=umap
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_15, target_AOIfilter,
    feature_results, "RegionAnimal.png",
    c("Region", "Animal"), 300, 3, UMAP=umap
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_15_log2, target_AOIfilter,
    feature_results, "logRegionRun.png",
    c("Region", "run"), 300, 3
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_15_log2, target_AOIfilter,
    feature_results, "logRegionOutliers.png",
    c("Region", "Outliers"), 300, 3,
    pca = pca
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_15_log2, target_AOIfilter,
    feature_results, "logRegioncondition.png",
    c("Region", "condition"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrixAOI_no_QC_15_log2, target_AOIfilter,
    feature_results, "logRegionAnimal.png",
    c("Region", "Animal"), 300, 3, pca=pca
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_15_log2, target_AOIfilter,
    feature_results, "logRegionRun.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
    norm_matrixAOI_no_QC_15_log2, target_AOIfilter,
    feature_results, "logRegionOutliers.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_15_log2, target_AOIfilter,
    feature_results, "logRegioncondition.png",
    c("Region", "condition"), 300, 3, UMAP=umap
)

umap <- umap_fun(
   norm_matrixAOI_no_QC_15_log2, target_AOIfilter,
    feature_results, "logRegionAnimal.png",
    c("Region", "Animal"), 300, 3, UMAP=umap
)

######## --Clustering with variable genes --########

#calculate the coefficient of variation (CV) for each gene (g) using the formula CVg=SDg/meang. 
variable_gene_results <- "1.variable_gene_results"
generate_folder(variable_gene_results)
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- apply(norm_matrix_log, MARGIN = 1, calc_CV)
CV_dat_sort <- sort(CV_dat, decreasing =TRUE)
write.csv(CV_dat_sort, file.path(variable_gene_results, "variablegenes.csv"))
GOI <- names(CV_dat_sort)[CV_dat_sort > quantile(CV_dat_sort, 0.9)]
write.csv(GOI, file.path(variable_gene_results, "variablegenes90thquartile.csv"))

message("STATUS: Number of genes above the 90th quartile ", length(GOI))
ds <- distanceMatrix(t(norm_matrix_log[GOI, ]), metric = "pearson")
hc <- hclust(ds, method = "ward.D2")
x <- cutree(hc, k =5)
geneorder <- hc$labels[hc$order]
clustergenes <- x[geneorder]
P36 <- createPalette(length(unique(clustergenes)), c("#ff0000", "#00ff00", "#0000ff"))
names(P36) <- unique(clustergenes)
clust <- data.frame(clustergenes)
colnames(clust) <- c("cluster")
P36animal <- createPalette(length(unique(target_roifilter10final$Animal)), c("#ff0000", "#00ff00", "#0000ff"))
names(P36animal) <- unique(target_roifilter10final$Animal)

colorslist <- c()
colornames <- c()
for (i in clustergenes) {
    colorslist <- c(colorslist, P36[[i]])
    colornames <- c(colornames, color.id(P36[[i]])[1])
}
annot_col <- list(
    run=c("Run1"="#0b830b", "Run2"="orange"),
    condition = c("Z" = "black", "C" = "gray"),
    Region=c("DWM"='red', "GM"='darkblue',"SWM"="purple"),
    Animal=P36animal,
    cluster = P36
)
pl <- pheatmap(norm_matrix_log[GOI, ], 
    scale = "row", show_rownames = FALSE, show_colnames = FALSE,
    border_color = NA,  clustering_method = "ward.D2",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    breaks = seq(-3, 3, 0.05), annotation_row = clust, annotation_colors  =annot_col,
    color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
    annotation_col = target_roifilter10final[, c("Animal", "condition", "Region", "run")],
    filename = file.path(variable_gene_results, "variablegenehm_0.9.png"))

clust$cluster <- colornames

write.csv(clust, file.path(variable_gene_results, "all_clusters.csv"))

if (isTRUE(SetRankRun)) {

    # allGenesHGNC <- unique(rownames(norm_matrix))
    # referenceSet <- symbol2EntrezID(allGenesHGNC)

    # collection <- buildSetCollection(allDBs,
    #     referenceSet = referenceSet,
    #     maxSetSize = 500
    # )
    # saveRDS(collection, "collectionallDBs.rds")
    collection <- readRDS("../collectionallDBs.rds")
}
clust$genes <- rownames(clust)
if (isTRUE(SetRankRun)) {
    message("STATUS: Finding gene enrichments for each cluster in DE")
    for (cluster in unique(clust$cluster)) {
        print(paste0("STATUS: gene enrichments for module ", cluster))
        genes <- clust[clust$cluster==cluster,"genes"]
        if (length(genes) > 0) {
            network <- gene_enrichment(genes, file.path(variable_gene_results, "SetRank_results/"), cluster)
        }
    }
}

######## --Based on PCA and looking of DSP with Caleb I am taking the remaining C15 samples out of the analysis --########
message("STATUS: removing animal C15 for technical reasons and it appears as on outlier the PCA plot ")
target_roifilter10final<- target_roifilter10final[target_roifilter10final$Animal!="C15",]
norm_matrix <- norm_matrix[, rownames(target_roifilter10final)]
write.csv(norm_matrix, file.path(norm_results, "1.norm_matrix_removedoutliers.csv"))
write.csv(target_roifilter10final, file.path(norm_results, "target_roifilter10final_removedoutliers.csv"))

norm_matrix_log <- norm_matrix_log[, rownames(target_roifilter10final)]
write.csv(norm_matrix_log, file.path(norm_results, "1.norm_matrix_log_removedoutliers.csv"))

######## --Feature Reduction with outliers removed--########
feature_results <- "1.feature_red_removedoutliers"
generate_folder(feature_results)

pca <- pca_fun(
    norm_matrix, target_roifilter10final,
    feature_results, "RegionRun.png",
    c("Region", "run"), 300, 3
)
pca <- pca_fun(
    norm_matrix, target_roifilter10final,
    feature_results, "RegionRun.pdf",
    c("Region", "run"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrix, target_roifilter10final,
    feature_results, "Regioncondition.png",
    c("Region", "condition"), 300, 3, pca=pca
)
pca <- pca_fun(
    norm_matrix, target_roifilter10final,
    feature_results, "Regioncondition.pdf",
    c("Region", "condition"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrix, target_roifilter10final,
    feature_results, "RegionAnimal.png",
    c("Region", "Animal"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrix, target_roifilter10final,
    feature_results, "RegionAnimal.pdf",
    c("Region", "Animal"), 300, 3, pca=pca
)

umap <- umap_fun(
   norm_matrix, target_roifilter10final,
    feature_results, "RegionRun.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
   norm_matrix, target_roifilter10final,
    feature_results, "Regioncondition.png",
    c("Region", "condition"), 300, 3, UMAP=umap
)

umap <- umap_fun(
   norm_matrix, target_roifilter10final,
    feature_results, "RegionAnimal.png",
    c("Region", "Animal"), 300, 3, UMAP=umap
)

pca <- pca_fun(
    norm_matrix_log, target_roifilter10final,
    feature_results, "logRegionRun.png",
    c("Region", "run"), 300, 3
)
pca <- pca_fun(
    norm_matrix_log, target_roifilter10final,
    feature_results, "logRegionRun.pdf",
    c("Region", "run"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrix_log, target_roifilter10final,
    feature_results, "logRegioncondition.png",
    c("Region", "condition"), 300, 3, pca=pca
)
pca <- pca_fun(
    norm_matrix_log, target_roifilter10final,
    feature_results, "logRegioncondition.pdf",
    c("Region", "condition"), 300, 3, pca=pca
)

pca <- pca_fun(
    norm_matrix_log, target_roifilter10final,
    feature_results, "logRegionAnimal.png",
    c("Region", "Animal"), 300, 3, pca=pca
)
pca <- pca_fun(
    norm_matrix_log, target_roifilter10final,
    feature_results, "logRegionAnimal.pdf",
    c("Region", "Animal"), 300, 3, pca=pca
)

umap <- umap_fun(
   norm_matrix_log, target_roifilter10final,
    feature_results, "logRegionRun.png",
    c("Region", "run"), 300, 3
)

umap <- umap_fun(
   norm_matrix_log, target_roifilter10final,
    feature_results, "logRegioncondition.png",
    c("Region", "condition"), 300, 3, UMAP=umap
)

umap <- umap_fun(
   norm_matrix_log, target_roifilter10final,
    feature_results, "logRegionAnimal.png",
    c("Region", "Animal"), 300, 3, UMAP=umap
)
######## --Clustering with variable genes with outliers removed--########

#calculate the coefficient of variation (CV) for each gene (g) using the formula CVg=SDg/meang. 
variable_gene_results <- "1.variable_gene_results_removedoutliers"
generate_folder(variable_gene_results)
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- apply(norm_matrix_log, MARGIN = 1, calc_CV)
CV_dat_sort <- sort(CV_dat, decreasing =TRUE)
write.csv(CV_dat_sort, file.path(variable_gene_results, "variablegenes.csv"))
GOI <- names(CV_dat_sort)[CV_dat_sort > quantile(CV_dat_sort, 0.9)]
write.csv(GOI, file.path(variable_gene_results, "variablegenes90thquartile.csv"))

message("STATUS: Number of genes above the 90th quartile ", length(GOI))
ds <- distanceMatrix(t(norm_matrix_log[GOI, ]), metric = "pearson")
hc <- hclust(ds, method = "ward.D2")
x <- cutree(hc, k =5)
geneorder <- hc$labels[hc$order]
clustergenes <- x[geneorder]
P36 <- createPalette(length(unique(clustergenes)), c("#ff0000", "#00ff00", "#0000ff"))

names(P36) <- unique(clustergenes)
clust <- data.frame(clustergenes)
colnames(clust) <- c("cluster")
colorslist <- c()
colornames <- c()
for (i in clustergenes) {
    colorslist <- c(colorslist, P36[[i]])
    colornames <- c(colornames, color.id(P36[[i]])[1])
}
annot_col <- list(
    run=c("Run1"="#0b830b", "Run2"="orange"),
    condition = c("Z" = "black", "C" = "gray"),
    Region=c("DWM"='red', "GM"='darkblue',"SWM"="purple"),
    Animal=P36animal,
    cluster = P36
)
pl <- pheatmap(norm_matrix_log[GOI, ], 
    scale = "row", show_rownames = FALSE, show_colnames = FALSE,
    border_color = NA,  clustering_method = "ward.D2",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    breaks = seq(-3, 3, 0.05), annotation_row = clust, annotation_colors  =annot_col,
    color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
    annotation_col = target_roifilter10final[, c("Animal", "condition", "Region", "run")],
    filename = file.path(variable_gene_results, "variablegenehm_0.9.png"))

clust$cluster <- colornames

write.csv(clust, file.path(variable_gene_results, "all_clusters.csv"))

if (isTRUE(SetRankRun)) {

    # allGenesHGNC <- unique(rownames(norm_matrix))
    # referenceSet <- symbol2EntrezID(allGenesHGNC)

    # collection <- buildSetCollection(allDBs,
    #     referenceSet = referenceSet,
    #     maxSetSize = 500
    # )
    # saveRDS(collection, "collectionallDBs.rds")
    collection <- readRDS("../collectionallDBs.rds")
}
clust$genes <- rownames(clust)
if (isTRUE(SetRankRun)) {
    message("STATUS: Finding gene enrichments for each cluster in DE")
    for (cluster in unique(clust$cluster)) {
        print(paste0("STATUS: gene enrichments for module ", cluster))
        genes <- clust[clust$cluster==cluster,"genes"]
        if (length(genes) > 0) {
            network <- gene_enrichment(genes, file.path(variable_gene_results, "SetRank_results/"), cluster)
        }
    }
}
######## --Clustering with variable genes AOI --########

#calculate the coefficient of variation (CV) for each gene (g) using the formula CVg=SDg/meang. 
variable_gene_results <- "1.variable_gene_AOI_results"
generate_folder(variable_gene_results)
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- apply(norm_matrixAOI_log, MARGIN = 1, calc_CV)
CV_dat_sort <- sort(CV_dat, decreasing =TRUE)
write.csv(CV_dat_sort, file.path(variable_gene_results, "variablegenes.csv"))
GOI <- names(CV_dat_sort)[CV_dat_sort > quantile(CV_dat_sort, 0.9)]
write.csv(GOI, file.path(variable_gene_results, "variablegenes90thquartile.csv"))

message("STATUS: Number of genes above the 90th quartile ", length(GOI))
ds <- distanceMatrix(t(norm_matrixAOI_log[GOI, ]), metric = "pearson")
hc <- hclust(ds, method = "ward.D2")
x <- cutree(hc, k =5)
geneorder <- hc$labels[hc$order]
clustergenes <- x[geneorder]
P36 <- createPalette(length(unique(clustergenes)), c("#ff0000", "#00ff00", "#0000ff"))
names(P36) <- unique(clustergenes)
clust <- data.frame(clustergenes)
colnames(clust) <- c("cluster")
P36animalAOI <- createPalette(length(unique(target_AOIfilter10final$Animal)), c("#ff0000", "#00ff00", "#0000ff"))
names(P36animalAOI) <- unique(target_AOIfilter10final$Animal)

colorslist <- c()
colornames <- c()
for (i in clustergenes) {
    colorslist <- c(colorslist, P36[[i]])
    colornames <- c(colornames, color.id(P36[[i]])[1])
}
annot_col <- list(
    run=c("Run1"="#0b830b", "Run2"="orange"),
    condition = c("Z" = "black", "C" = "gray"),
    Region=c("OLIG2"='#00d9ff', "NEUN"="pink","GFAP"="#20f09d"),
    Animal=P36animalAOI,
    cluster = P36
)
pl <- pheatmap(norm_matrixAOI_log[GOI, ], 
    scale = "row", show_rownames = FALSE, show_colnames = FALSE,
    border_color = NA,  clustering_method = "ward.D2",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    breaks = seq(-3, 3, 0.05), annotation_row = clust, annotation_colors  =annot_col,
    color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
    annotation_col = target_AOIfilter10final[, c("Animal", "condition", "Region", "run")],
    filename = file.path(variable_gene_results, "variablegenehm_0.9.png"))

clust$cluster <- colornames

write.csv(clust, file.path(variable_gene_results, "all_clusters.csv"))
######## --Differential Gene Expression AOI --########
de_results <- "1.de_results_full_AOI"
generate_folder(de_results)
modelFormula <- ~ condition + (1 | slidename)
if (isTRUE(all.equal(colnames(norm_matrixAOI_no_QC_25_log2), rownames(target_AOIfilter)))) {
    results <- diff_gene_exprs(norm_matrixAOI_no_QC_25_log2, target_AOIfilter,
        modelFormula = modelFormula, subset = "Region", contrast.var = "condition",
        contrast.levels = c("Z", "C"))
} else { 
    stop('ERROR: Columns and rownames for AOIs are not in the right order')
}
write.csv(results, file.path(de_results, "DEresults.csv"))

# visualization de volcano plots
results$Color <- "NS or LFC < 0.26"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.26] <- "NS or LFC < 0.26"
results$Color <- factor(results$Color,
    levels = c(
        "NS or FC < 0.5", "P < 0.05",
        "FDR < 0.05", "FDR < 0.001"
    )
)
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()

for (region in unique(target_AOIfilter[, "Region"])) {
    ind <- results$Subset == region
    top_g <- c(
        top_g,
        results[ind, "Gene"][
            order(results[ind, "invert_P"], decreasing = TRUE)[1:15]
        ],
        results[ind, "Gene"][
            order(results[ind, "invert_P"], decreasing = FALSE)[1:15]
        ]
    )
}

top_g <- unique(top_g)
results <- results[, -1 * ncol(results)] # remove invert_P from matrix
# Graph results
p <- ggplot(
    results,
    aes(
        x = Estimate, y = -log10(`Pr(>|t|)`),
        color = Color, label = Gene
    )
) +
    geom_vline(xintercept = c(0.26, -0.26), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point(size = 0.6) +
    labs(
        x = "Zika vs Control log2(FC)",
        y = "Significance, -log10(P)",
        color = "Significance"
    ) +
    scale_color_manual(
        values = c(
            `FDR < 0.001` = "dodgerblue",
            `FDR < 0.05` = "lightblue",
            `P < 0.05` = "orange2",
            `NS or LFC < 0.26` = "gray"
        ),
        guide = guide_legend(override.aes = list(size = 4))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    geom_text_repel(
        data = subset(results, Gene %in% top_g & FDR < 0.05),
        size = 2.3, point.padding = 0.15, color = "black",
        min.segment.length = .1, box.padding = .2, lwd = 2,
        max.overlaps = 10
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom") +
    facet_wrap(~Subset, scales = "free_y") +
    theme(strip.background = element_rect(fill = "white")) # c("DWM"='red', "GM"='darkblue',"SWM"="purple")
ggsave(file.path(de_results, "volcano.png"), width = 8, height = 4, units = "in", dpi = 300, bg = "white")
ggsave(file.path(de_results, "volcano.svg"), width = 8, height = 4, units = "in", dpi = 300, bg = "white")

# Subset by just FDR
GOIlfcup <- unique(subset(results, `Estimate` > 0.26)$Gene)
GOIlfcdwn <- unique(subset(results, `Estimate` < -0.26)$Gene)
GOIlfc <- c(GOIlfcup, GOIlfcdwn)
GOIpval <- unique(subset(results, `FDR` < 0.05)$Gene)
GOI <- intersect(GOIlfc, GOIpval)
message("STATUS: Number of significant genes FDR < 0.05 & abs(LFC) > 0.26, ", length(GOI))
resultssig <- results[results$Gene %in% GOI, ]
siglfcs <- dcast(resultssig, Gene ~ Subset, value.var = "Estimate")
sigpvals <- dcast(resultssig, Gene ~ Subset, value.var = "FDR")
sigpvalues <- dcast(resultssig, Gene ~ Subset, value.var = "Pr(>|t|)")

rownames(siglfcs) <- siglfcs$Gene
siglfcs$Gene <- NULL

rownames(sigpvals) <- sigpvals$Gene
sigpvals$Gene <- NULL

rownames(sigpvalues) <- sigpvalues$Gene
sigpvalues$Gene <- NULL

write.csv(siglfcs, file.path(de_results, "sigLFCsFDRcutoff.csv"))
write.csv(sigpvals, file.path(de_results, "sigFDRsFDRcutoff.csv"))
write.csv(sigpvalues, file.path(de_results, "sigpvalueFDRcutoff.csv"))

# Subset by just pval
GOIlfcup <- unique(subset(results, `Estimate` > 0.26)$Gene)
GOIlfcdwn <- unique(subset(results, `Estimate` < -0.26)$Gene)
GOIlfc <- c(GOIlfcup, GOIlfcdwn)
GOIpval <- unique(subset(results, `Pr(>|t|)` < 0.01)$Gene)
GOI <- intersect(GOIlfc, GOIpval)
message("STATUS: Number of significant genes FDR < 0.05 & abs(LFC) > 0.26, ", length(GOI))
resultssig <- results[results$Gene %in% GOI, ]
siglfcs <- dcast(resultssig, Gene ~ Subset, value.var = "Estimate")
sigpvals <- dcast(resultssig, Gene ~ Subset, value.var = "FDR")
sigpvalues <- dcast(resultssig, Gene ~ Subset, value.var = "Pr(>|t|)")

rownames(siglfcs) <- siglfcs$Gene
siglfcs$Gene <- NULL

rownames(sigpvals) <- sigpvals$Gene
sigpvals$Gene <- NULL

rownames(sigpvalues) <- sigpvalues$Gene
sigpvalues$Gene <- NULL

write.csv(siglfcs, file.path(de_results, "sigLFCspvalcutoff.csv"))
write.csv(sigpvals, file.path(de_results, "sigFDRspvalcutoff.csv"))
write.csv(sigpvalues, file.path(de_results, "sigpvaluepvalcutoff.csv"))

######## -- GSEA on each region AOI --########
GESA_de_results <- file.path("1.de_results_full_AOI", "GSEAresults")
generate_folder(GESA_de_results)
GMTFILE <- "/share/lwhitmo/GeneEnrichmentDBs/msigdb_v7.0_files_to_download_locally/msigdb_v7.0_GMTs/c5.bp.v7.0.symbols.gmt"

### GFAP
siglfcsGFAP <- data.frame(resultssig[(resultssig$Subset == "GFAP" & resultssig$`Pr(>|t|)` < 0.05), c("Gene", "Estimate")])
write.csv(siglfcsGFAP, file.path(de_results, "GFAP_sigLFCs.csv"))

siglfcsGFAPordered <- siglfcsGFAP[order(-siglfcsGFAP$Estimate), ]
rownames(siglfcsGFAPordered) <- siglfcsGFAPordered$Gene
siglfcsGFAPordered$Gene <- NULL
siglfcsGFAPorderedfinal <- as.numeric(unlist(siglfcsGFAPordered))
names(siglfcsGFAPorderedfinal) <- rownames(siglfcsGFAPordered)

ResultsGFAP <- run_gsea_and_ora(siglfcsGFAPorderedfinal, GMTFILE,
     rownames(norm_matrixAOI_no_QC_25_log2), "GFAP", GESA_de_results)
siggsea <- ResultsGFAP$sigGSEA
siggseanes <- siggsea[, c("pathway", "NES"), ]
siggseanes$pathway <- str_remove_all(siggseanes$pathway, "GO_")
siggseanes <- as.data.frame(siggseanes)
siggseanes <- siggseanes[order(siggseanes$NES), ]

rownames(siggseanes) <- siggseanes$pathway
siggseanes$pathway <- NULL

contrasting <- colorRampPalette(rev(c("darkgreen", "green", "white", "purple", "#380c53")))(100)
myBreaks <- c(
    seq(min(-3), 0, length.out = ceiling(100 / 2) + 1),
    seq(max(3) / 100, max(3), length.out = floor(100 / 2))
)

pl <- pheatmap(t(siggseanes),
    # scale = "row",
    show_rownames = TRUE, show_colnames = TRUE,
    border_color = T, cluster_rows = FALSE, cluster_cols = FALSE,
    breaks = myBreaks,
    color = contrasting,
    filename = file.path(GESA_de_results, "sigGSEA_DE_GFAP.png"),
    width = 6, height = 4
)

### OLIG2
siglfcsOLIG2 <- data.frame(resultssig[(resultssig$Subset == "OLIG2" & resultssig$`Pr(>|t|)` < 0.05), c("Gene", "Estimate")])
write.csv(siglfcsOLIG2, file.path(de_results, "OLIG2_sigLFCs.csv"))
siglfcsOLIG2ordered <- siglfcsOLIG2[order(-siglfcsOLIG2$Estimate), ]
rownames(siglfcsOLIG2ordered) <- siglfcsOLIG2ordered$Gene
siglfcsOLIG2ordered$Gene <- NULL
siglfcsOLIG2orderedfinal <- as.numeric(unlist(siglfcsOLIG2ordered))
names(siglfcsOLIG2orderedfinal) <- rownames(siglfcsOLIG2ordered)

ResultsOLIG2 <- run_gsea_and_ora(siglfcsOLIG2orderedfinal, GMTFILE, rownames(norm_matrix_log), "OLIG2", GESA_de_results)

siggsea <- ResultsOLIG2$sigGSEA
siggseanes <- siggsea[, c("pathway", "NES"), ]
siggseanes$pathway <- str_remove_all(siggseanes$pathway, "GO_")
siggseanes <- as.data.frame(siggseanes)
siggseanes <- siggseanes[order(siggseanes$NES), ]
rownames(siggseanes) <- siggseanes$pathway
siggseanes$pathway <- NULL

# contrasting <- colorRampPalette(rev(c("darkgreen", "green", "white", "purple", "#380c53")))(100)

# pl <- pheatmap(t(siggseanes),
#     # scale = "row",
#     show_rownames = TRUE, show_colnames = TRUE,
#     border_color = T, cluster_rows = FALSE, cluster_cols = FALSE,
#     breaks = myBreaks,
#     color = contrasting, fontsize_col = 6,
#     filename = file.path(GESA_de_results, "sigGSEA_DE_OLIG2.png"),
#     width = 8, height = 4
# )

### NeuN
siglfcsNeuN <- data.frame(resultssig[(resultssig$Subset == "NeuN" & resultssig$`Pr(>|t|)` < 0.05), c("Gene", "Estimate")])
write.csv(siglfcsNeuN, file.path(de_results, "NeuN_sigLFCs.csv"))
siglfcsNeuNordered <- siglfcsNeuN[order(-siglfcsNeuN$Estimate), ]
rownames(siglfcsNeuNordered) <- siglfcsNeuNordered$Gene
siglfcsNeuNordered$Gene <- NULL
siglfcsNeuNorderedfinal <- as.numeric(unlist(siglfcsNeuNordered))
names(siglfcsNeuNorderedfinal) <- rownames(siglfcsNeuNordered)

ResultsNeuN <- run_gsea_and_ora(siglfcsNeuNorderedfinal, GMTFILE, rownames(norm_matrix_log), "NeuN", GESA_de_results)

siggsea <- ResultsNeuN$sigGSEA
siggseanes <- siggsea[, c("pathway", "NES"), ]
siggseanes$pathway <- str_remove_all(siggseanes$pathway, "GO_")
siggseanes <- as.data.frame(siggseanes)
siggseanes <- siggseanes[order(siggseanes$NES), ]
rownames(siggseanes) <- siggseanes$pathway
siggseanes$pathway <- NULL

# contrasting <- colorRampPalette(rev(c("darkgreen", "green", "white", "purple", "#380c53")))(100)

# pl <- pheatmap(t(siggseanes),
#     # scale = "row",
#     show_rownames = TRUE, show_colnames = TRUE,
#     border_color = T, cluster_rows = FALSE, cluster_cols = FALSE,
#     breaks = myBreaks,
#     color = contrasting, fontsize_col = 6,
#     filename = file.path(GESA_de_results, "sigGSEA_DE_NeuN.png"),
#     width = 8, height = 4
# )

######## --Differential Gene Expression --########
de_results <- "1.de_results_full"
generate_folder(de_results)
modelFormula = ~ condition + (1 | run)
results <- diff_gene_exprs(norm_matrix_log, target_roifilter10final,
        modelFormula = modelFormula, subset="Region", contrast.var = "condition",
        contrast.levels = c( "Z", "C"))
write.csv(results, file.path(de_results, "DEresults.csv"))

#visualization de volcano plots 
results$Color <- "NS or LFC < 0.26"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.26] <- "NS or LFC < 0.26"
results$Color <- factor(results$Color,
    levels = c(
        "NS or FC < 0.5", "P < 0.05",
        "FDR < 0.05", "FDR < 0.001"
    )
)
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()

for(region in unique(target_roifilter10final[, "Region"])) {
    ind <- results$Subset == region
    top_g <- c(top_g,
               results[ind, 'Gene'][
                   order(results[ind, 'invert_P'], decreasing = TRUE)[1:15]],
               results[ind, 'Gene'][
                   order(results[ind, 'invert_P'], decreasing = FALSE)[1:15]])
}

top_g <- unique(top_g)
results <- results[, -1 * ncol(results)] # remove invert_P from matrix
# Graph results
p<- ggplot(
    results,
    aes(
        x = Estimate, y = -log10(`Pr(>|t|)`),
        color = Color, label = Gene
    )
) +
    geom_vline(xintercept = c(0.26, -0.26), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point(size=0.6) +
    labs(
        x = "Zika vs Control log2(FC)",
        y = "Significance, -log10(P)",
        color = "Significance"
    ) +
    scale_color_manual(
        values = c(
            `FDR < 0.001` = "dodgerblue",
            `FDR < 0.05` = "lightblue",
            `P < 0.05` = "orange2",
            `NS or LFC < 0.26` = "gray"
        ),
        guide = guide_legend(override.aes = list(size = 4))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    geom_text_repel(
        data = subset(results, Gene %in% top_g & FDR < 0.05),
        size = 2.3, point.padding = 0.15, color = "black",
        min.segment.length = .1, box.padding = .2, lwd = 2,
        max.overlaps = 10
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom") +
    facet_wrap(~Subset, scales = "free_y") +
    theme(strip.background = element_rect(fill ="white")) #c("DWM"='red', "GM"='darkblue',"SWM"="purple")
ggsave(file.path(de_results, "volcano.png"), width=8, height=4, units="in", dpi=300, bg="white")
ggsave(file.path(de_results, "volcano.svg"), width = 8, height = 4, units = "in", dpi = 300, bg="white")

#visualization de volcano plots 
results$Color <- "NS or LFC < 0.26"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.26] <- "NS or LFC < 0.26"
results$Color <- factor(results$Color,
    levels = c(
        "NS or FC < 0.5", "P < 0.05",
        "FDR < 0.05", "FDR < 0.001"
    )
)
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()

for(region in unique(target_roifilter10final[, "Region"])) {
    ind <- results$Subset == region
    top_g <- c(top_g,
               results[ind, 'Gene'][
                   order(results[ind, 'invert_P'], decreasing = TRUE)[1:15]],
               results[ind, 'Gene'][
                   order(results[ind, 'invert_P'], decreasing = FALSE)[1:15]])
}

top_g <- unique(top_g)
results <- results[, -1 * ncol(results)] # remove invert_P from matrix
# Graph results
p<- ggplot(
    results,
    aes(
        x = Estimate, y = -log10(`Pr(>|t|)`),
        color = Color, label = Gene
    )
) +
    geom_vline(xintercept = c(0.26, -0.26), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    geom_point(size=0.6) +
    labs(
        x = "Zika vs Control log2(FC)",
        y = "Significance, -log10(P)",
        color = "Significance"
    ) +
    scale_color_manual(
        values = c(
            `FDR < 0.001` = "dodgerblue",
            `FDR < 0.05` = "lightblue",
            `P < 0.05` = "orange2",
            `NS or LFC < 0.26` = "gray"
        ),
        guide = guide_legend(override.aes = list(size = 4))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    geom_text_repel(
        data = subset(results, Gene %in% top_g & FDR < 0.05),
        size = 2.3, point.padding = 0.15, color = "black",
        min.segment.length = .1, box.padding = .2, lwd = 2,
        max.overlaps = 10
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom") +
    facet_wrap(~Subset, scales = "free_y") +
    theme(strip.background = element_rect(fill ="white")) #c("DWM"='red', "GM"='darkblue',"SWM"="purple")
ggsave(file.path(de_results, "volcano.png"), width=8, height=4, units="in", dpi=300, bg="white")
ggsave(file.path(de_results, "volcano.svg"), width = 8, height = 4, units = "in", dpi = 300, bg="white")

# Heatmap of significant genes norm values
GOIlfcup <- unique(subset(results, `Estimate` > 0.26)$Gene)
GOIlfcdwn <- unique(subset(results, `Estimate` < -0.26)$Gene)
GOIlfc <- c(GOIlfcup, GOIlfcdwn)
GOIpval <- unique(subset(results, `FDR` < 0.05)$Gene)
GOI <- intersect(GOIlfc, GOIpval)
message("STATUS: Number of significant genes FDR < 0.05 & abs(LFC) > 0.26, ", length(GOI))
ds <- distanceMatrix(t(norm_matrix_log[GOI, ]), metric = "pearson")
hc <- hclust(ds, method = "ward.D2")
x <- cutree(hc, k = 5)
geneorder <- hc$labels[hc$order]
clustergenes <- x[geneorder]
P36 <- createPalette(length(unique(clustergenes)), c("#ff0000", "#00ff00", "#0000ff"))

names(P36) <- unique(clustergenes)
clust <- data.frame(clustergenes)
colnames(clust) <- c("cluster")
colorslist <- c()
colornames <- c()
for (i in clustergenes) {
    colorslist <- c(colorslist, P36[[i]])
    colornames <- c(colornames, color.id(P36[[i]])[1])
}
annot_col <- list(
    run = c("Run1" = "#0b830b", "Run2" = "orange"),
    condition = c("Z" = "black", "C" = "gray"),
    Region = c("DWM" = "red", "GM" = "darkblue", "SWM" = "purple"),
    Animal = P36animal,
    cluster = P36
)
pl <- pheatmap(norm_matrix_log[GOI, ],
    scale = "row", show_rownames = FALSE, show_colnames = FALSE,
    border_color = NA, clustering_method = "ward.D2",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    breaks = seq(-3, 3, 0.05), annotation_row = clust, annotation_colors = annot_col,
    color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
    annotation_col = target_roifilter10final[, c("Animal", "condition", "Region", "run")],
    filename = file.path(de_results, "DEgenesnormvalues.png")
)
annot_col <- list(
    run = c("Run1" = "#0b830b", "Run2" = "orange"),
    condition = c("Z" = "black", "C" = "gray"),
    Region = c("DWM" = "red", "GM" = "darkblue", "SWM" = "purple"),
    Animal = P36animal
)
target_roifilter10finalorder <- target_roifilter10final[order(target_roifilter10final$Region, 
    target_roifilter10final$condition, target_roifilter10final$Animal),]
pl <- pheatmap(norm_matrix_log[GOI, rownames(target_roifilter10finalorder) ],
    scale = "row", show_rownames = FALSE, show_colnames = FALSE,
    border_color = T, clustering_method = "ward.D2",
    cluster_cols = F,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    gaps_col=c(28, 66),
    breaks = seq(-3, 3, 0.05), annotation_colors = annot_col,
    color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
    annotation_col = target_roifilter10final[, c("Animal", "condition", "Region", "run")],
    filename = file.path(de_results, "DEgenesnormvaluesOrdered.png")
)

# Heatmap of significant LFCs at least in one brain region 
GOIlfcup <- unique(subset(results, `Estimate` > 0.26)$Gene)
GOIlfcdwn <- unique(subset(results, `Estimate` < -0.26)$Gene)
GOIlfc <- c(GOIlfcup, GOIlfcdwn)
GOIpval <- unique(subset(results, `FDR` < 0.05)$Gene)
GOI <- intersect(GOIlfc, GOIpval)
message("STATUS: Number of significant genes FDR < 0.05 & abs(LFC) > 0.26, ", length(GOI))
resultssig <- results[results$Gene %in% GOI, ]
siglfcs <- dcast(resultssig, Gene ~ Subset,  value.var = "Estimate")
sigpvals <- dcast(resultssig, Gene ~ Subset, value.var = "FDR")
sigpvalues <- dcast(resultssig, Gene ~ Subset, value.var = "Pr(>|t|)")

rownames(siglfcs) <- siglfcs$Gene
siglfcs$Gene <- NULL

rownames(sigpvals) <- sigpvals$Gene
sigpvals$Gene <- NULL

rownames(sigpvalues) <- sigpvalues$Gene
sigpvalues$Gene <- NULL

write.csv(siglfcs, file.path(de_results, "sigLFCs.csv"))
write.csv(sigpvals, file.path(de_results, "sigFDRs.csv"))
write.csv(sigpvalues, file.path(de_results, "sigpvalue.csv"))

siglfcs <- as.matrix(data.frame(siglfcs))
class(siglfcs) <- "numeric"
ds <- distanceMatrix(t(siglfcs), metric = "pearson")
hc <- hclust(ds, method = "ward.D2")
x <- cutree(hc, k = 4)
geneorder <- hc$labels[hc$order]
clustergenes <- x[geneorder]
P36 <- createPalette(length(unique(clustergenes)), c("#ff0000", "#00ff00", "#0000ff"))
names(P36) <- unique(clustergenes)
clust <- data.frame(clustergenes)
colnames(clust) <- c("cluster")
colorslist <- c()
colornames <- c()
for (i in clustergenes) {
    colorslist <- c(colorslist, P36[[i]])
    colornames <- c(colornames, color.id(P36[[i]])[1])
}
annot_col <- list(
    cluster = P36
)
contrasting <- colorRampPalette(rev(c("red", "orange", "white", "lightblue", "blue")))(100)
myBreaks <- c(
    seq(min(-1), 0, length.out = ceiling(100 / 2) + 1),
    seq(max(1) / 100, max(1), length.out = floor(100 / 2))
)
pl <- pheatmap(siglfcs,
    # scale = "row", 
    show_rownames = FALSE, show_colnames = TRUE,
    border_color = NA, clustering_method = "ward.D2",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    treeheight_col = 20,
    breaks = myBreaks, 
    color = contrasting, 
    annotation_colors = annot_col,  
    annotation_row = clust,
    filename = file.path(de_results, "DEgenesLFC.pdf"),
    width = 3.5, height = 8
)

clust$cluster <- colornames
write.csv(clust, file.path(de_results, "all_clusters.csv"))
clust$genes <- rownames(clust)
if (isTRUE(SetRankRun)) {
    message("STATUS: Finding gene enrichments for each cluster in DE")
    for (cluster in unique(clust$cluster)) {
        print(paste0("STATUS: gene enrichments for module ", cluster))
        genes <- clust[clust$cluster == cluster, "genes"]
        if (length(genes) > 0) {
            network <- gene_enrichment(genes, file.path(de_results, "SetRank_results/"), cluster)
        }
    }
}
siglfcsclusterorder <- siglfcs[clust$genes, ]
sigpvalsclusterorder <- sigpvalues[clust$genes, ]
sigpvalsFDRclusterorder <- sigpvals[clust$genes, ]
if ((isTRUE(all.equal(rownames(siglfcsclusterorder), rownames(sigpvalsclusterorder)))) && (isTRUE(all.equal(rownames(sigpvalsFDRclusterorder), clust$genes)))) {
    tot <- cbind(clust, siglfcsclusterorder, sigpvalsclusterorder, sigpvalsFDRclusterorder)
    # tot$genes <- NULL
    colnames(tot) <- c("cluster", "LFC.DWM", "LFC.GM", "LFC.SWM", "pval.DWM", "pval.GM",
         "pval.SWM", "FDR.DWM", "FDR.GM", "FDR.SWM", "cluster")
    write.csv(tot, file.path(de_results, "suppde.csv"))
} else {
    stop("ERROR: gene names not in the same order")
}

######## -- Generate files for cytoscape  -- ########
IDConverter <- createIDConverter(
    "Homo.sapiens", "ENTREZID",
    "SYMBOL"
)
cytoscape_files_results <- "Cytoscapefiles"
generate_folder(cytoscape_files_results)

# get pathways of interest 
axonguidancepath <- allDBs[allDBs$termID=="GO:0007411", "geneID"]
axonguidancepath <- IDConverter(as.character(axonguidancepath))
axonpath <- allDBs[allDBs$termID=="GO:0030424", "geneID"]
axonpath <- IDConverter(as.character(axonpath))
myelinsheathpath <-  allDBs[allDBs$termID=="GO:0043209", "geneID"]
myelinsheathpath <- IDConverter(as.character(myelinsheathpath))
myelinpath <-  allDBs[allDBs$termID=="GO:0042552", "geneID"]
myelinpath <- IDConverter(as.character(myelinpath))

pinkmem <- read.table(
    "/share/lwhitmo/projects/Gale_GeoMax_JennyGo/WhitmoreAnalysis/DSP_Q3_Analysis_removedC15/1.de_results_full/SetRank_results/de_unranked_mistyrose2_membership.txt", 
     header=TRUE)
axonguid = pinkmem[pinkmem$GO.0007411=="X", c("gene", "GO.0007411")]$gene
axon = pinkmem[pinkmem$GO.0030424=="X", c("gene", "GO.0030424")]$gene
write(paste0("Gene\tInteraction\tPathway"), file.path(cytoscape_files_results, "cytoscapepathways.sif"))
for (g in axonguid) {
    print(g[1])
    if (g %in% axonguidancepath) {
        write(paste0(g,"\tinteracts with\t", "Axon Guidance"), 
        file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }  else { 
        message("WARNING: ", g, "not in in axon guidance path which it should ")
    }
    if (g %in% axonpath) {
        write(paste0(g,"\tinteracts with\t", "Axon"),file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }
    if (g %in% myelinpath) {
        write(paste0(g,"\tinteracts with\t", "Myelin"), file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }
    if (g %in% myelinsheathpath) {
        write(paste0(g,"\tinteracts with\t", "Myelin Sheath"), file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }
}
for (g in axon) {
    if (g %in% axonpath) {
        write(paste0(g,"\tinteracts with\t", "Axon"), 
        file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }  else { 
        message("WARNING: ", g, "not in in axon path which it should ")
    }
    if (g %in% axonguidancepath) {
        write(paste0(g,"\tinteracts with\t", "Axon Guidance"),file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }
    if (g %in% myelinpath) {
        write(paste0(g,"\tinteracts with\t", "Myelin"), file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }
    if (g %in% myelinsheathpath) {
        write(paste0(g,"\tinteracts with\t", "Myelin Sheath"), file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }
}
bluemem <- read.table(
    "/share/lwhitmo/projects/Gale_GeoMax_JennyGo/WhitmoreAnalysis/DSP_Q3_Analysis_removedC15/1.de_results_full/SetRank_results/de_unranked_blue_membership.txt",
     header=T)
myelinsheath = bluemem[bluemem$GO.0043209=="X", c("gene", "GO.0043209")]$gene
myelin = bluemem[bluemem$GO.0042552=="X", c("gene", "GO.0042552")]$gene
for (g in myelinsheath) {
    if (g %in% myelinsheathpath) {
        write(paste0(g,"\tinteracts with\t", "Myelin Sheath"), 
        file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }  else { 
        message("WARNING: ", g, "not in myelin sheath path which it should ")
    }
    if (g %in% axonguidancepath) {
        write(paste0(g,"\tinteracts with\t", "Axon Guidance"),
        file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }
    if (g %in% myelinpath) {
        write(paste0(g,"\tinteracts with\t", "Myelin"), 
        file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }
    if (g %in% axonpath) {
        write(paste0(g,"\tinteracts with\t", "Axon"), 
        file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }
}
for (g in myelin) {
    if (g %in% myelinpath) {
        write(paste0(g,"\tinteracts with\t", "Myelin"), 
        file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }  else { 
        message("WARNING: ", g, "not in myelin path which it should ")
    }
    if (g %in% axonguidancepath) {
        write(paste0(g,"\tinteracts with\t", "Axon Guidance"),
        file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }
    if (g %in% myelinsheathpath) {
        write(paste0(g,"\tinteracts with\t", "Myelin Sheath"), 
        file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }
    if (g %in% axonpath) {
        write(paste0(g,"\tinteracts with\t", "Axon"), 
        file.path(cytoscape_files_results, "cytoscapepathways.sif"), append=TRUE)
    }
}
lfcscyto <- siglfcs[unique(c(axon, axonguid, myelinsheath, myelin)), c("DWM", "GM")]
pvalscyto <-  sigpvalsclusterorder[unique(c(axon, axonguid, myelinsheath, myelin)), c("DWM", "GM")]
fdrpvalscyto <-  sigpvalsFDRclusterorder[unique(c(axon, axonguid, myelinsheath, myelin)), c("DWM", "GM")]
neglog2fdrpvalcyto <- -log2(fdrpvalscyto)
totalcytoscape <- cbind(lfcscyto, pvalscyto, fdrpvalscyto, neglog2fdrpvalcyto)

colnames(totalcytoscape) <- c("lfc.DWM", "lfc.GM", "pval.DWM", "pval.GM", 
    "FDR.DWM", "FDR.GM", 
    "-log2(FDR.DWM)", "-log2(FDR.GM)")
write.csv(totalcytoscape, file.path(cytoscape_files_results, "lfcs4network.csv"))

######## -- GSEA on each region--########
GESA_de_results <- file.path("1.de_results_full", "GSEAresults")
generate_folder(GESA_de_results)
GMTFILE <- "/share/lwhitmo/GeneEnrichmentDBs/msigdb_v7.0_files_to_download_locally/msigdb_v7.0_GMTs/c5.bp.v7.0.symbols.gmt"

###DWM
siglfcsDWM <- data.frame(resultssig[(resultssig$Subset=="DWM" & resultssig$FDR < 0.05), c("Gene", "Estimate")])
write.csv(siglfcsDWM, file.path(de_results, "DWM_sigLFCs.csv"))

siglfcsDWMordered <- siglfcsDWM[order(-siglfcsDWM$Estimate), ]
rownames(siglfcsDWMordered) <- siglfcsDWMordered$Gene
siglfcsDWMordered$Gene <- NULL
siglfcsDWMorderedfinal <- as.numeric(unlist(siglfcsDWMordered))
names(siglfcsDWMorderedfinal) <- rownames(siglfcsDWMordered)

ResultsDWM <- run_gsea_and_ora(siglfcsDWMorderedfinal, GMTFILE, rownames(norm_matrix_log), "DWM", GESA_de_results)
siggsea <- ResultsDWM$sigGSEA
siggseanes <- siggsea[, c("pathway", "NES"), ]
siggseanes$pathway <- str_remove_all(siggseanes$pathway, "GO_")
siggseanes <- as.data.frame(siggseanes)
siggseanes <- siggseanes[order(siggseanes$NES), ]

rownames(siggseanes) <- siggseanes$pathway
siggseanes$pathway <- NULL

contrasting <- colorRampPalette(rev(c("darkgreen", "green", "white", "purple", "#380c53")))(100)
myBreaks <- c(
    seq(min(-3), 0, length.out = ceiling(100 / 2) + 1),
    seq(max(3) / 100, max(3), length.out = floor(100 / 2))
)

pl <- pheatmap(t(siggseanes),
    # scale = "row",
    show_rownames = TRUE, show_colnames = TRUE,
    border_color = T, cluster_rows=FALSE, cluster_cols=FALSE,
    breaks = myBreaks,
    color = contrasting,
    filename = file.path(GESA_de_results, "sigGSEA_DE_DWM.png"),
    width = 6, height = 4
)

### SWM
siglfcsSWM <- data.frame(resultssig[(resultssig$Subset == "SWM" & resultssig$FDR < 0.05), c("Gene", "Estimate")])
write.csv(siglfcsSWM, file.path(de_results, "SWM_sigLFCs.csv"))

siglfcsSWMordered <- siglfcsSWM[order(-siglfcsSWM$Estimate), ]
rownames(siglfcsSWMordered) <- siglfcsSWMordered$Gene
siglfcsSWMordered$Gene <- NULL
siglfcsSWMorderedfinal <- as.numeric(unlist(siglfcsSWMordered))
names(siglfcsSWMorderedfinal) <- rownames(siglfcsSWMordered)

ResultsSWM <- run_gsea_and_ora(siglfcsSWMorderedfinal, GMTFILE, rownames(norm_matrix_log), "SWM", GESA_de_results)

### GM
siglfcsGM <- data.frame(resultssig[(resultssig$Subset == "GM" & resultssig$FDR < 0.05), c("Gene", "Estimate")])
write.csv(siglfcsGM, file.path(de_results, "GM_sigLFCs.csv"))
siglfcsGMordered <- siglfcsGM[order(-siglfcsGM$Estimate), ]
rownames(siglfcsGMordered) <- siglfcsGMordered$Gene
siglfcsGMordered$Gene <- NULL
siglfcsGMorderedfinal <- as.numeric(unlist(siglfcsGMordered))
names(siglfcsGMorderedfinal) <- rownames(siglfcsGMordered)

ResultsGM <- run_gsea_and_ora(siglfcsGMorderedfinal, GMTFILE, rownames(norm_matrix_log), "GM", GESA_de_results)

siggsea <- ResultsGM$sigGSEA
siggseanes <- siggsea[, c("pathway", "NES"), ]
siggseanes$pathway <- str_remove_all(siggseanes$pathway, "GO_")
siggseanes <- as.data.frame(siggseanes)
siggseanes <- siggseanes[order(siggseanes$NES),]
rownames(siggseanes) <- siggseanes$pathway
siggseanes$pathway <- NULL

contrasting <- colorRampPalette(rev(c("darkgreen", "green", "white", "purple", "#380c53")))(100)

pl <- pheatmap(t(siggseanes),
    # scale = "row",
    show_rownames = TRUE, show_colnames = TRUE,
    border_color = T, cluster_rows = FALSE, cluster_cols = FALSE,
    breaks = myBreaks,
    color = contrasting, fontsize_col=6,
    filename = file.path(GESA_de_results, "sigGSEA_DE_GM.png"),
    width = 8, height = 4
)

# total plot 
pathwaysuniq <- unique(c(ResultsDWM$sigGSEA$pathway, ResultsGM$sigGSEA$pathway))
dwmgsea <- ResultsDWM$GSEA[ResultsDWM$GSEA$pathway %in% pathwaysuniq, c(1:3, 6)]
dwmgsea$region <- rep("DWM", length(rownames(dwmgsea)))

gmgsea <- ResultsGM$GSEA[ResultsGM$GSEA$pathway %in% pathwaysuniq, c(1:3, 6)]
gmgsea$region <- rep("GM", length(rownames(gmgsea)))

totalgsea <- rbind(dwmgsea, gmgsea)
totalgsea$pathway <- str_remove_all(totalgsea$pathway, "GO_")
totalgsea$pathway <- str_replace_all(totalgsea$pathway, "_", " ")
totalgsea$logpval <- -log10(totalgsea$padj)
ggplot(totalgsea, aes(x=pathway, y=NES)) +
    geom_line(aes(group = pathway)) +
    geom_point(aes(color = logpval, shape=region), size=3) +
    geom_hline(yintercept=0, linetype="dashed") +
    # scale_colour_manual(values =c("DWM" = "red", "GM" = "darkblue")) +
    theme_Publication() + labs(x="") + ylim(-3, 4) +
    theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)
    )
ggsave(file.path(GESA_de_results, "GM_DWMGSEAplot.png"), width = 7.5, height = 5, units = "in", dpi = 500)
ggsave(file.path(GESA_de_results, "GM_DWMGSEAplot.svg"), width = 7.5, height = 5, units = "in", dpi = 500)
ggsave(file.path(GESA_de_results, "GM_DWMGSEAplot.pdf"), width = 7.5, height = 5, units = "in", dpi = 500)

ggplot(totalgsea, aes(x=pathway, y=NES)) +
    geom_line(aes(group = pathway)) +
    geom_point(aes(color = logpval, shape=region), size=3) +
    geom_hline(yintercept=0, linetype="dashed") +
    # scale_colour_manual(values =c("DWM" = "red", "GM" = "darkblue")) +
    theme_Publication() + labs(x="", color="-log10(adj pvalue)") + ylim(-3, 4) +
    theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)
    ) + coord_flip() + theme(legend.position="bottom", legend.direction = "horizontal", 
         legend.key.size = unit(0.4, "cm"), legend.title=element_text(size=8), legend.text=element_text(size=8))
ggsave(file.path(GESA_de_results, "GM_DWMGSEAplot_flipped.png"), width = 8.3, height = 8.5, units = "in", dpi = 500)
ggsave(file.path(GESA_de_results, "GM_DWMGSEAplot_flipped.svg"), width = 8.3, height = 8.5, units = "in", dpi = 500)
ggsave(file.path(GESA_de_results, "GM_DWMGSEAplot_flipped.pdf"), width = 8.3, height = 8.5, units = "in", dpi = 500)

write.csv(totalgsea, file.path(GESA_de_results, "totalgsearesults4plot.csv"))

######## -- Marker Gene Exploration --########
marker_results <- "1.marker_results"
generate_folder(marker_results)

modelFormula = ~ Region + (1 | run)
# DWM vs SWM
DWMvsSWMresults <- diff_gene_exprs(norm_matrix_log, target_roifilter10final,
        modelFormula = modelFormula, subset=NULL, contrast.var = "Region",
        contrast.levels = c( "DWM", "SWM"))
write.csv(DWMvsSWMresults, file.path(marker_results, "DEresultsDWMvsSWM.csv"))

DWMvsSWMresultsSig <- DWMvsSWMresults[DWMvsSWMresults$FDR < 0.05, ]
DWMvsSWMresultsSig <- DWMvsSWMresultsSig[order(DWMvsSWMresultsSig$Estimate),]

# DWM vs GM
DWMvsGMresults <- diff_gene_exprs(norm_matrix_log, target_roifilter10final,
        modelFormula = modelFormula, subset=NULL, contrast.var = "Region",
        contrast.levels = c( "DWM", "GM"))
write.csv(DWMvsGMresults, file.path(marker_results, "DEresultsDWMvsGM.csv"))

DWMvsGMresultsSig <- DWMvsGMresults[DWMvsGMresults$FDR < 0.05, ]
DWMvsGMresultsSig <- DWMvsGMresultsSig[order(DWMvsGMresultsSig$Estimate),]

# SWM vs GM
SWMvsGMresults <- diff_gene_exprs(norm_matrix_log, target_roifilter10final,
        modelFormula = modelFormula, subset=NULL, contrast.var = "Region",
        contrast.levels = c( "SWM", "GM"))
write.csv(SWMvsGMresults, file.path(marker_results, "DEresultsSWMvsGM.csv"))

SWMvsGMresultsSig <- SWMvsGMresults[SWMvsGMresults$FDR < 0.05, ]
SWMvsGMresultsSig <- SWMvsGMresultsSig[order(SWMvsGMresultsSig$Estimate),]

SWMmarkers <- head(DWMvsSWMresultsSig, 10)[,"Gene"]
DWMmarkers <- tail(DWMvsSWMresultsSig, 10)[,"Gene"]
GMmarkers <- head(DWMvsGMresultsSig, 10)[,"Gene"]
DWMmarkers1 <- tail(DWMvsGMresultsSig, 10)[,"Gene"]
GMmarkers1 <- head(SWMvsGMresultsSig, 10)[,"Gene"]
SWMmarkers1 <- tail(SWMvsGMresultsSig, 10)[,"Gene"]

GMmarkersfinal <- unique(c(GMmarkers, GMmarkers1))
message("STATUS: number of markers for GM is ", length(GMmarkersfinal))
DWMmarkersfinal <- unique(c(DWMmarkers, DWMmarkers1))
message("STATUS: number of markers for DWM is ", length(DWMmarkersfinal))

SWMmarkersfinal <- unique(c(SWMmarkers, SWMmarkers1))
message("STATUS: number of markers for DWM is ", length(SWMmarkersfinal))

markers <- unique(c(GMmarkersfinal, DWMmarkersfinal, SWMmarkersfinal))
message("STATUS: number of markers for all regions is ", length(markers))
write.csv(markers, file.path(marker_results, "1.location_markers.csv"))
contrasting3 <- colorRampPalette(rev(c("darkblue", "mediumblue", "dodgerblue", "white", "orange", "red", "darkred")))(100)
breaksList <- seq(-4, 4, by = .1)

annot_col <- list(
    run = c("Run1" = "#0b830b", "Run2" = "orange"),
    condition = c("Z" = "black", "C" = "gray"),
    Region = c("DWM" = "red", "GM" = "darkblue", "SWM" = "purple"),
    Animal = P36animal
)

# png(file.path(marker_results, "MarkerGene_heatmap_colnames.png"), width = 15, height = 10, units = "in", res = 600)
pheatmap(norm_matrix_log[markers, ],
    scale = "row",
    clustering_method = "ward.D2",
    color = colorRampPalette(rev(c(name = contrasting3)))(length(breaksList)),
    cluster_rows = T,
    breaks = breaksList,
    cluster_cols = T,
    show_rownames = T,
    show_colnames = F,
    border_color = "black",
    angle_col = "90",
    # cellwidth = 10,
    # cellheight = 15,
    fontsize_number = 12,
    fontsize = 12,
    cutree_cols = 3,
    cutree_rows=2,
    # gaps_col = c(6,12,18),
    # annotation = annot,
    annotation_col = target_roifilter10final[, c("Animal", "condition", "Region", "run")],
    filename = file.path(marker_results, "MarkerGene_heatmap_colnames.pdf"),
    width = 15, height = 10,
    annotation_colors  =annot_col,
    main = paste("")
)
# dev.off()

### AOI z-score heatmap of highlighted genes from Fig S3
results_folder <- "1.AOI_goi_results"
generate_folder(results_folder)
genes <- c("TUBA4A", "TUBB4A", "HACD2", "IL13",  "ATF2", "HDAC1", "RASD2", "IL36G", "CCL2", 'DCX', "CCL11", "MBP",
"FOS", "MOBP", "PLP1", "MAG", "TOMM20", "TOMM20", "LY6H", "HBA1", "TOMM70", "TRPM7", "TLR3", "FASLG")

genesinter <- intersect(rownames(norm_matrixAOI_log), genes)
notfoundgenes <- setdiff(genes, rownames(norm_matrixAOI_log))

P36animalAOI <- createPalette(length(unique(target_AOIfilter10final$Animal)), c("#ff0000", "#00ff00", "#0000ff"))
target_AOIfilter10finalsort <- target_AOIfilter10final[order(target_AOIfilter10final$Region, target_AOIfilter10final$condition), ]
names(P36animalAOI) <- unique(target_AOIfilter10final$Animal)

annot_col <- list(
    condition = c("Z" = "black", "C" = "gray"),
    Region = c("GFAP" = "green", "NEUN" ="#4e4ee8", "OLIG2" ="#7714b5"),
    Animal = P36animalAOI
)
contrasting3 <- colorRampPalette(rev(c("darkblue", "mediumblue", "dodgerblue", "white", "orange", "red", "darkred")))(100)
breaksList <- seq(-2, 2, by = .1)

pheatmap(norm_matrixAOI_log[genesinter, rownames(target_AOIfilter10finalsort)],
    scale = "row",
    clustering_method = "ward.D2",
    color = colorRampPalette(rev(c(name = contrasting3)))(length(breaksList)),
    cluster_rows = T,
    breaks = breaksList,
    cluster_cols = F,
    show_rownames = T,
    show_colnames = F,
    border_color = "black",
    angle_col = "90",
    # cellwidth = 10,
    # cellheight = 15,
    fontsize_number = 12,
    fontsize = 12,
    # cutree_cols = 3,
    gaps_col = c(5, 12),

    cutree_rows=2,
    # gaps_col = c(6,12,18),
    # annotation = annot,
    annotation_col = target_AOIfilter10final[, c("Animal", "condition", "Region")],
    filename = file.path(results_folder, "AOI_GOI.png"),
    width = 8, height =6,
    annotation_colors  =annot_col,
    main = paste("")
)

genesinter <- intersect(rownames(norm_matrixAOI_no_QC_25_log2), genes)
target_AOIfilterfinalsort <- target_AOIfilter[order(target_AOIfilter$Region, target_AOIfilter$condition), ]
P36animalAOI <- createPalette(length(unique(target_AOIfilterfinalsort$Animal)), c("#ff0000", "#00ff00", "#0000ff"))
names(P36animalAOI) <- unique(target_AOIfilterfinalsort$Animal)

annot_col <- list(
    Outliers=c("Outlier"="red", "No Outlier"="black"),
    condition = c("Z" = "black", "C" = "gray"),
    Region = c("GFAP" = "green", "NEUN" = "#4e4ee8", "OLIG2" = "#7714b5"),
    Animal = P36animalAOI
)
pheatmap(norm_matrixAOI_no_QC_25_log2[genesinter, rownames(target_AOIfilterfinalsort)],
    scale = "row",
    clustering_method = "ward.D2",
    color = colorRampPalette(rev(c(name = contrasting3)))(length(breaksList)),
    cluster_rows = T,
    breaks = breaksList,
    cluster_cols = F,
    show_rownames = T,
    show_colnames = F,
    border_color = "black",
    angle_col = "90",
    # cellwidth = 10,
    # cellheight = 15,
    fontsize_number = 12,
    fontsize = 12,
    # cutree_cols = 3,
    gaps_col = c(9, 16),
    cutree_rows = 2,
    # gaps_col = c(6,12,18),
    # annotation = annot,
    annotation_col = target_AOIfilterfinalsort[, c("Animal", "condition", "Region", "Outliers")],
    filename = file.path(results_folder, "AOI_GOI_no_QC_25_log2.png"),
    width = 8, height = 6,
    annotation_colors = annot_col,
    main = paste("")
)
# dev.off()

######## --Network genes-- ########
genesinter <- intersect(rownames(totalcytoscape), rownames(norm_matrixAOI_no_QC_25_log2))
target_AOIfilterfinalsort <- target_AOIfilter[order(target_AOIfilter$Region, target_AOIfilter$condition), ]
P36animalAOI <- createPalette(length(unique(target_AOIfilterfinalsort$Animal)), c("#ff0000", "#00ff00", "#0000ff"))
names(P36animalAOI) <- unique(target_AOIfilterfinalsort$Animal)

annot_col <- list(
    Outliers = c("Outlier" = "red", "No Outlier" = "black"),
    condition = c("Z" = "black", "C" = "gray"),
    Region = c("GFAP" = "green", "NEUN" = "#4e4ee8", "OLIG2" = "#7714b5"),
    Animal = P36animalAOI
)
pheatmap(norm_matrixAOI_no_QC_25_log2[genesinter, rownames(target_AOIfilterfinalsort)],
    scale = "row",
    clustering_method = "ward.D2",
    color = colorRampPalette(rev(c(name = contrasting3)))(length(breaksList)),
    cluster_rows = T,
    breaks = breaksList,
    cluster_cols = F,
    show_rownames = T,
    show_colnames = F,
    border_color = "black",
    angle_col = "90",
    # cellwidth = 10,
    # cellheight = 15,
    fontsize_number = 12,
    fontsize = 8,
    # cutree_cols = 3,
    gaps_col = c(9, 16),
    cutree_rows = 3,
    # gaps_col = c(6,12,18),
    # annotation = annot,
    annotation_col = target_AOIfilterfinalsort[, c("Animal", "condition", "Region", "Outliers")],
    filename = file.path(results_folder, "ROINetworkgenes.png"),
    width = 8, height = 8,
    annotation_colors = annot_col,
    main = paste("")
)
# ## Include AOI's for validation
# totaltaraget <- rbind(target_roifilter10final, target_AOIfilter10final)
# P36animalall <- createPalette(length(unique(totaltaraget$Animal)), c("#ff0000", "#00ff00", "#0000ff"))
# names(P36animalall) <- unique(totaltaraget$Animal)
# annot_col <- list(
#     run = c("Run1" = "#0b830b", "Run2" = "orange"),
#     condition = c("Z" = "black", "C" = "gray"),
#     Region = c("DWM" = "red", "GM" = "darkblue", "SWM" = "purple",
#         "OLIG2"='#00d9ff', "NEUN"="pink","GFAP"="#20f09d"),
#     Animal = P36animalall
# )

# # png(file.path(marker_results, "MarkerGene_heatmap_colnames.png"), width = 15, height = 10, units = "in", res = 600)
# genesNOTinAOI <- setdiff(markers, rownames(norm_matrixAOI_log))
# genesinAOI <- intersect(rownames(norm_matrixAOI_log),markers)

# pheatmap(cbind(scale_rows(norm_matrix_log[genesinAOI, ]), scale_rows(norm_matrixAOI_log[genesinAOI, ])),
#     scale = "none",
#     clustering_method = "ward.D2",
#     color = colorRampPalette(rev(c(name = contrasting3)))(length(breaksList)),
#     cluster_rows = T,
#     breaks = breaksList,
#     cluster_cols = T,
#     show_rownames = T,
#     show_colnames = F,
#     border_color = "black",
#     angle_col = "90",
#     # cellwidth = 10,
#     # cellheight = 15,
#     fontsize_number = 12,
#     fontsize = 12,
#     cutree_cols = 5,
#     cutree_rows=2,
#     # gaps_col = c(6,12,18),
#     # annotation = annot,
#     annotation_col = totaltaraget[, c("Animal", "condition", "Region", "run")],
#     filename = file.path(marker_results, "MarkerGene_heatmap_colnameswithAOIs.png"),
#     width = 15, height = 10,
#     annotation_colors  =annot_col,
#     main = paste("")
# )

######## --Myelination --###########
# bluemem <- read.table("/share/lwhitmo/projects/Gale_GeoMax_JennyGo/WhitmoreAnalysis/DSP_Q3_Analysis_removedC15/1.de_results_full/SetRank_results/de_unranked_blue_membership.txt",
#                 header=T)
# myelinsheath = bluemem[bluemem$GO.0043209=="X", c("gene", "GO.0043209")]$gene
# myelin = bluemem[bluemem$GO.0042552=="X", c("gene", "GO.0042552")]$gene

# uniqmyelin <- setdiff(myelin, myelinsheath)
# inmyelin <- intersect(myelin, myelinsheath)
# uniqmyelinsheath <- setdiff(myelinsheath, myelin)

# myelingenes <- c(uniqmyelin, inmyelin, uniqmyelinsheath)
# pheatmap(norm_matrix_log[myelingenes, ],
#     scale = "row",
#     clustering_method = "ward.D2",
#     color = colorRampPalette(rev(c(name = contrasting3)))(length(breaksList)),
#     cluster_rows = T,
#     breaks = breaksList,
#     cluster_cols = T,
#     show_rownames = T,
#     show_colnames = F,
#     border_color = "black",
#     angle_col = "90",
#     # cellwidth = 10,
#     # cellheight = 15,
#     fontsize_number = 12,
#     fontsize = 12,
#     cutree_cols = 3,
#     cutree_rows=2,
#     # gaps_col = c(6,12,18),
#     # annotation = annot,
#     annotation_col = target_roifilter10final[, c("Animal", "condition", "Region", "run")],
#     filename = file.path(de_results, "DE_myelingenes.png"),
#     width = 15, height = 10,
#     annotation_colors  =annot_col,
#     main = paste("")
# )
# # dev.off()

# annot_col <- list(
#     run = c("Run1" = "#0b830b", "Run2" = "orange"),
#     condition = c("Z" = "black", "C" = "gray"),
#     Region = c("DWM" = "red", "GM" = "darkblue", "SWM" = "purple",
#         "OLIG2"='#00d9ff', "NEUN"="pink","GFAP"="#20f09d"),
#     Animal = P36animalall
# )
# pheatmap(cbind(scale_rows(norm_matrix_log[myelingenes, ]), scale_rows(norm_matrixAOI_log[myelingenes, ])),
#     scale = "none",
#     clustering_method = "ward.D2",
#     color = colorRampPalette(rev(c(name = contrasting3)))(length(breaksList)),
#     cluster_rows = T,
#     breaks = breaksList,
#     cluster_cols = T,
#     show_rownames = T,
#     show_colnames = F,
#     border_color = "black",
#     angle_col = "90",
#     # cellwidth = 10,
#     # cellheight = 15,
#     fontsize_number = 12,
#     fontsize = 12,
#     cutree_cols = 6,
#     cutree_rows=2,
#     # gaps_col = c(6,12,18),
#     # annotation = annot,
#     annotation_col = totaltaraget[, c("Animal", "condition", "Region", "run")],
#     filename = file.path(de_results, "DE_myelingeneswithAOIs.png"),
#     width = 15, height = 10,
#     annotation_colors  =annot_col,
#     main = paste("")
# )
# # dev.off()

# pheatmap(norm_matrix_log[myelingenes, row.names(target_roifilter10final[target_roifilter10final$Region=="DWM",]) ],
#     scale = "row",
#     clustering_method = "ward.D2",
#     color = colorRampPalette(rev(c(name = contrasting3)))(length(breaksList)),
#     cluster_rows = T,
#     breaks = breaksList,
#     cluster_cols = T,
#     show_rownames = T,
#     show_colnames = F,
#     border_color = "black",
#     angle_col = "90",
#     # cellwidth = 10,
#     # cellheight = 15,
#     fontsize_number = 12,
#     fontsize = 12,
#     cutree_cols = 2,
#     cutree_rows=2,
#     # gaps_col = c(6,12,18),
#     # annotation = annot,
#     annotation_col = target_roifilter10final[, c("Animal", "condition", "Region", "run")],
#     filename = file.path(de_results, "DE_myelingenesDWMonly.png"),
#     width = 10, height = 7,
#     annotation_colors  =annot_col,
#     main = paste("")
# )

# annot_col <- list(
#     run = c("Run1" = "#0b830b", "Run2" = "orange"),
#     condition = c("Z" = "black", "C" = "gray"),
#     Region = c("DWM" = "red", "OLIG2"='#00d9ff', "NEUN"="pink","GFAP"="#20f09d"),
#     Animal = P36animalall
# )

# totaltaragetreduced <- totaltaraget[!(totaltaraget$Region %in% c("SWM", "GM")),] 

# normmatrixreducedscaled <- cbind(scale_rows(norm_matrix_log[myelingenes, row.names(target_roifilter10final[target_roifilter10final$Region=="DWM",]) ]),
#         scale_rows(norm_matrixAOI_log[myelingenes, ]))

# # Scale AOIs and ROIs independently 
# pheatmap(normmatrixreducedscaled,
#     scale = "none",
#     clustering_method = "ward.D2",
#     color = colorRampPalette(rev(c(name = contrasting3)))(length(breaksList)),
#     cluster_rows = T,
#     breaks = breaksList,
#     cluster_cols = T,
#     show_rownames = T,
#     show_colnames = F,
#     border_color = "black",
#     angle_col = "90",
#     # cellwidth = 10,
#     # cellheight = 15,
#     fontsize_number = 12,
#     fontsize = 12,
#     cutree_cols = 3,
#     cutree_rows=2,
#     # gaps_col = c(6,12,18),
#     # annotation = annot,
#     annotation_col = totaltaragetreduced[, c("Animal", "condition", "Region", "run")],
#     filename = file.path(de_results, "DE_myelingenesDWMonlywithAOIS.png"),
#     width = 10, height = 7,
#     annotation_colors  =annot_col,
#     main = paste("")
# )

