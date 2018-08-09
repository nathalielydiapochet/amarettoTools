 importReg = function(regfn) {
  read.csv(regfn, stringsAsFactors=FALSE)
 }
 
 importExprs = function(expfn, genesOnly=TRUE, colClasses=NA, ...) {
  if (!genesOnly) stop("currently only importing gene symbols")
  read.csv(expfn, stringsAsFactors=FALSE, colClasses=colClasses, ...)
 }

#' produce a graph that summarizes regulator-module relationships
#' @param dataFolder character(1) specifying location of regulators 
#' and gene_expression CSV files by module, with filenames of form
#' 'Module_[n]_...' by default; other filename formats can be
#' specified by modpatt parameter.  `_gene_expression.csv` and
#' `_regulators.csv` are the assumed suffixes.
#' @param modpatt character(1) regular expression used to pick out the CSV files
#' @return instance of graph::graphNEL, mode 'directed'
#' @note The returned graph has regulators, module names,
#' and module gene symbols as nodes.  A module has edges to
#' all member genes, and a regulator has edges to all modules
#' that it regulates.
#' @examples
#' data(regulGBM)
#' regulGBM
#' @export
makeRegulatorGraph = function(dataFolder=".", modpatt="^Module_") {
 csvnames = dir(dataFolder, pattern=modpatt, full=TRUE)
 ncsv = length(csvnames)
 if ((ncsv %% 2)!=0) warning("number of CSV files is not even")
 regnames = grep("_regulators.csv", csvnames, value=TRUE)
 exprnames = grep("_gene_expression.csv", csvnames, value=TRUE)
 allregs = lapply(regnames, importReg)
 regsyms = lapply(allregs, "[[", 1)
 modnums = gsub(".*Module_(.*)_gene_expression.csv", "\\1", exprnames)
 modnames = paste0("Module_", modnums)
 names(allregs) = modnames
 names(regsyms) = modnames
#
# probe an expression file to determine colClasses
#
 expr1 = read.csv(grep("gene_expression", csvnames, value=TRUE)[1],
    nrows=1)
 NSAMP = ncol(expr1)-1
 fl = c("character", rep("NULL", NSAMP))
#
# now just read gene names from module files
#
 allexprs = lapply(exprnames, importExprs, colClasses=fl)
 names(allexprs) = modnames
 allregsyms = unique(unlist(regsyms))
 nn = c(allregsyms, names(regsyms)) # node names
 regul = new("graphNEL", edgemode="directed", nodes=nn)
 for (i in modnames) regul <- addEdge(allregs[[i]][,1], i,
     regul, weights = allregs[[i]][,2])
 for (j in modnames) {
     regul <- addNode(setdiff(allexprs[[j]][[1]], 
          nodes(regul)), regul)
     regul <- addEdge(j, allexprs[[j]][[1]], regul)
     }
 regul
}
