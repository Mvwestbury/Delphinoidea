
pathnode <- function (phylo, tipsonly = T)
{
    require(phangorn)
    di.tr <- dist.nodes(phylo)
    root.tr <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[,
        2])][1]
    tr.depth <- max(di.tr[as.numeric(colnames(di.tr)) == root.tr,
        ])
    if (tipsonly == TRUE) {
        roottotippath <- di.tr[as.numeric(rownames(di.tr)) ==
            root.tr, 1:length(phylo$tip.label)]
        nodesinpath <- sapply(1:length(phylo$tip.label), function(x) length(Ancestors(phylo,
            x)))
    }
    else {
        roottotippath <- di.tr[as.numeric(rownames(di.tr)) ==
            root.tr, ]
        nodesinpath <- sapply(1:(length(phylo$tip.label) + phylo$Nnode),
            function(x) length(Ancestors(phylo, x)))
    }
    #plot(roottotippath, nodesinpath, xlab = "Root-to-tip path length", ylab = "Number of parent nodes", pch = 20)
    return(list(roottotippath = roottotippath, nodesinpath = nodesinpath))
}

stemmy <- function(tre) sum(tre$edge.length[which(tre$edge[,2] > Ntip(tre))], na.rm = T) / sum(tre$edge.length, na.rm = T)

rtt.cov <- function(tre){
	if(!is.rooted(tre)) tre <- midpoint(tre)
        rtts <- pathnode(tre)[[1]]
        rttcov <- sd(rtts, na.rm = T) / mean(rtts, na.rm = T)
        return(rttcov)
}

pis <- function (x, what = "fraction", use.ambiguities = FALSE){
    if (!inherits(x, "DNAbin"))
        stop("'x' is not of class 'DNAbin'")
    what <- match.arg(what, c("absolute", "fraction", "index"))
    if (use.ambiguities) {
        warning("'use.ambiguities' is currently ignored ", "and IUPAC ambiguity symbols are treated as missing data")
        use.ambiguities <- FALSE
    }
    pars.inf <- function(x) {
        x <- table(x)
        x <- x[x > 1]
        n <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w",
            "y")
        if (length(x[!names(x) %in% n]) > 1)
            x <- TRUE
        else x <- FALSE
    }
    x <- as.character(x)
    out <- apply(x, 2, pars.inf)
    if (what %in% c("absolute", "fraction")) {
        out <- length(out[out])
        if (what == "fraction") {
            out <- round(out/ncol(x) * 100, digits = 2)
        }
    }
    else {
        out <- which(out)
    }
    out
}
