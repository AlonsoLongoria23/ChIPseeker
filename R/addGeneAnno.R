##' get gene annotation, symbol, gene name etc.
##'
##'
##' @title getGeneAnno
##' @param annoDb annotation package
##' @param geneID query geneID
##' @param type gene ID type
##' @param columns names of columns to be obtained from database
##' @return data.frame
##' @importFrom AnnotationDbi select
##' @author G Yu
getGeneAnno <- function(annoDb, geneID, type, columns){
    kk <- unlist(geneID)
    require(annoDb, character.only = TRUE)
    annoDb <- eval(parse(text = annoDb))

    ## Custom OrgDb packages may use organism-specific gene IDs (e.g. GID)
    ## that do not map cleanly onto ChIPseeker's inferred `type`.
    ## For these packages, use the package-native keytype directly.
    if (annoDb$packageName %in% c("org.Dpulex.eg.db",
                                  "org.Tthymallus.eg.db",
                                  "org.Tarcticus.eg.db")) {
        kt <- "GID"
    } else if (type == "Entrez Gene ID") {
        kt <- "ENTREZID"
    } else if (type == "Ensembl gene ID" || type == "Ensembl Gene ID") {
        kt <- "ENSEMBL"
    } else {
        message("geneID type is not supported...\tPlease report it to developer...\n")
        return(NA)
    }

    if (annoDb$packageName == "org.Dpulex.eg.db") {
        Dpulex_kk <- AnnotationDbi::mapIds(
            annoDb,
            keys = kk,
            keytype = "SYMBOL",
            column = "GID"
        )

        ann <- tryCatch(
            suppressWarnings(AnnotationDbi::select(
                annoDb,
                keys = Dpulex_kk,
                keytype = kt,
                columns = columns
            )),
            error = function(e) NULL
        )

        if (is.null(ann)) {
            warning("ID type not matched, gene annotation will not be added...")
            return(NA)
        }

        return(ann)

    } else if (annoDb$packageName %in% c("org.Tthymallus.eg.db", "org.Tarcticus.eg.db")) {
        ann <- tryCatch(
            suppressWarnings(AnnotationDbi::select(
                annoDb,
                keys = kk,
                keytype = kt,
                columns = columns
            )),
            error = function(e) NULL
        )

        if (is.null(ann)) {
            warning("ID type not matched, gene annotation will not be added...")
            return(NA)
        }

        return(ann)

    } else {
        i <- which(!is.na(kk))
        kk <- gsub("\\.\\d+$", "", kk)

        ann <- tryCatch(
            suppressWarnings(select(
                annoDb,
                keys = unique(kk[i]),
                keytype = kt,
                columns = columns
            )),
            error = function(e) NULL
        )

        if (is.null(ann)) {
            warning("ID type not matched, gene annotation will not be added...")
            return(NA)
        }

        idx <- getFirstHitIndex(ann[, kt])
        ann <- ann[idx, ]

        rownames(ann) <- ann[, kt]
        res <- ann[as.character(kk), ]
        return(res)
    }
}


addGeneAnno <- function(peak.gr, annoDb, type, columns) {
    geneAnno <- getGeneAnno(annoDb, peak.gr$geneId, type, columns)
    if (!all(is.na(geneAnno))) {
        for (cn in colnames(geneAnno)[-1]) {
            mcols(peak.gr)[[cn]] <- geneAnno[, cn]
        }
    }
    return(peak.gr)
}
