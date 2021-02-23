#' @importFrom RPostgreSQL dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @import RPostgreSQL
#' @import org.Hs.eg.db
#' @importFrom methods new
#'
#' @title GeneHancerDB
#------------------------------------------------------------------------------------------------------------------------
#' @name GeneHancerDB-class
#' @rdname GeneHancerDB-class
#' @aliases GeneHancerDB
#'
#' @import methods

# GENEHANCER_VERSION <- "gh411"
GENEHANCER_VERSION <- "gh50"

.GeneHancerDB <- setClass("GeneHancerDB",
                          representation = representation(
                             db="character",
                             state="environment"
                             )
                          )
#------------------------------------------------------------------------------------------------------------------------
setGeneric('retrieveEnhancersFromDatabase', signature='obj', function(obj, targetGene, tissues)
              standardGeneric('retrieveEnhancersFromDatabase'))
setGeneric('listTissues', signature='obj', function(obj, targetGene) standardGeneric('listTissues'))
setGeneric('getEnhancerTissues', signature='obj', function(obj, targetGene) standardGeneric ('getEnhancerTissues'))
setGeneric('getEnhancers',  signature='obj', function(obj, targetGene, tissues="all", maxSize=10000) standardGeneric ('getEnhancers'))
setGeneric('queryByRegion', signature='obj', function(obj, chrom, start, end) standardGeneric('queryByRegion'))
#------------------------------------------------------------------------------------------------------------------------
#' Create a GeneHancerDB connection
#'
#' @rdname GeneHancerDB-class
#'
#' @return An object of the GeneHancerDB class
#'
#' @export
#'
GeneHancerDB <- function()
{
   db <- NA_character_;
   state <- new.env(parent=emptyenv())
   .GeneHancerDB(db=db, state=state)

} # ctor
#------------------------------------------------------------------------------------------------------------------------
#' Get all the enhancer regions for the gene
#'
#' @rdname retrieveEnhancersFromDatabase
#' @aliases retrieveEnhancersFromDatabase
#'
#' @param obj An object class GeneHancerDB
#' @param targetGene a HUGO gene symbol
#' @param tissues "all" or a vector of case-agnostic tissue names
#'
#' @seealso listTissues
#'
#' @export

setMethod('retrieveEnhancersFromDatabase',  'GeneHancerDB',

     function(obj, targetGene, tissues){

        if(length(tissues) == 1 & tissues[1] == "all")
           tissueClause <- ""
        else {
           tissueSet <- paste(tissues, collapse="','")
           tissueClause <- sprintf("AND t.tissue in ('%s') ", tissueSet)
           }

        query <- paste0("select e.chr as chrom, ",
                        "e.element_start as start, ",
                        "e.element_end as end, ",
                        "a.symbol as gene, ",
                        "a.eqtl_score as eqtl, ",
                        "a.chic_score as HiC, ",
                        "a.erna_score as erna, ",
                        "a.expression_score as coexpression, ",
                        "a.distance_score as distanceScore, ",
                        "a.tss_proximity as tssProximity, ",
                        "a.combined_score as combinedScore, ",
                        "a.is_elite as elite, ",
                        "t.source as source, ",
                        "t.tissue as tissue, ",
                        "e.type as type, ",
                        "a.ghid as ghid ",
                        "from associations AS a, ",
                        "tissues AS t, elements as e ",
                        "where a.symbol='%s' ",
                        "%s",
                        "AND a.ghid=t.ghid ",
                        "AND e.ghid=a.ghid")
        query <- sprintf(query, targetGene, tissueClause)

        db <- dbConnect(PostgreSQL(), user= "trena", password="trena",
                        dbname=GENEHANCER_VERSION, #"gh411",
                        host="khaleesi.systemsbiology.net", port="5432")
        tbl <- dbGetQuery(db, query)
        dbDisconnect(db)

        if(nrow(tbl) == 0){
           warning(sprintf("no GeneHancer regions for %s in tissues %s", targetGene, paste(tissues, collapse=",")))
           return(data.frame())
           }

        tbl$sig <- with(tbl, sprintf("%s:%d-%d", chrom, start, end))
        tbl.trimmed <- .eliminateDupsCollapseTissues(tbl)

          # our current best guess is that eQTL, Hi-C, and enhancer RNA are credible indicators
          # of enhancer/gene association.  so keep only the rows with a value in one or more
          # of these columns, or with a combined score > 5.
          # combinedscore is some unstated function of all the scores.  we include as a fallback
          # an alternative threshold, just in case.

        tbl.2 <- subset(tbl.trimmed, !(is.nan(eqtl) & is.nan(hic) & is.nan(erna)) | combinedscore >= 5)
        return(tbl.2)
        })

#------------------------------------------------------------------------------------------------------------------------
.eliminateDupsCollapseTissues <- function(tbl)
{
   tbl.2 <- tbl
   sig.uniq <- unique(tbl.2$sig)
   sig.census <- lapply(sig.uniq, function(sig) grep(sig, tbl.2$sig))
   names(sig.census) <- sig.uniq

   tissues.by.sig <- lapply(sig.uniq, function(sig) tbl.2[grep(sig, tbl.2$sig), "tissue"])
   tissues.collapsed.by.sig <- lapply(tissues.by.sig, function(tissues) paste(tissues, collapse=";"))
   names(tissues.collapsed.by.sig) <- sig.uniq

   dups <- which(duplicated(tbl.2$sig))

   tbl.3 <- tbl.2   # optimistic, remains true if no dups in tabple

   if(length(dups) > 0){
      tbl.3 <- tbl.2[-dups,]
      tissues.collapsed.by.sig
      length(tissues.collapsed.by.sig)
      indices <- unlist(lapply(names(tissues.collapsed.by.sig), function(sig) grep(sig, tbl.3$sig)))
      tbl.3$tissue[indices] <- as.character(tissues.collapsed.by.sig)
      }

   coi <- c("chrom","start","end","gene","eqtl","hic","erna","coexpression","distancescore","tssproximity","combinedscore","elite","source","type","ghid","tissue")
   tbl.4 <- tbl.3[, coi]
   invisible(tbl.4)

} # .eliminateDupsCollapseTissues
#------------------------------------------------------------------------------------------------------------------------
#' Return a character vector containing all of the tissues known to GeneHancer
#'
#' @rdname listTissues
#' @aliases listTissues
#'
#' @param obj An object of class GeneHancerDB
#' @param targetGene A character string, default NA in which case all tissues for all genes are returned
#'
#' @export

setMethod('listTissues', 'GeneHancerDB',

    function(obj, targetGene){
       query <- "select distinct tissue from tissues"
       if(!is.na(targetGene)){
          query.p1 <- "select distinct t.tissue from associations as a, tissues as t where "
          query.p2 <- sprintf("a.symbol='%s' AND a.ghid=t.ghid", targetGene)
          query <- paste0(query.p1, query.p2)
          }
       db <- dbConnect(PostgreSQL(), user= "trena", password="trena",
                       dbname=GENEHANCER_VERSION, #"gh411",
                       host="khaleesi")
       result <- dbGetQuery(db, query)$tissue
       dbDisconnect(db)
       return(result)
       })

#------------------------------------------------------------------------------------------------------------------------
#' Get all the enhancer tissues included in the current genehancer
#'
#' @rdname getEnhancerTissues
#' @aliases getEnhancerTissues
#'
#' @param obj An object of class GeneHancerDB
#' @param targetGene a character string
#'
#' @seealso getEnhancers
#'
#' @export
setMethod('getEnhancerTissues',  'GeneHancerDB',

     function(obj, targetGene){
        listTissues(obj, targetGene)
        })

#------------------------------------------------------------------------------------------------------------------------
#' Get all the enhancer regions for the gene
#'
#' @rdname getEnhancers
#' @aliases getEnhancers
#'
#' @param obj An object of class GeneHancerDB
#' @param targetGene default NA, in which case the current object's targetGene is used.
#'
#' @seealso setTargetGene
#'
#' @export

setMethod('getEnhancers',  'GeneHancerDB',

     function(obj, targetGene, tissues="all", maxSize=10000){
        if(is.null(targetGene)) return(data.frame())
        tbl <- retrieveEnhancersFromDatabase(obj, targetGene, tissues)
        if(nrow(tbl) == 0)
           return(data.frame())
        size <- with(tbl, 1 + end - start)
        deleters <- which(size > maxSize)
        if(length(deleters) > 0)
           tbl <- tbl[-deleters,]
        tbl
        })

#------------------------------------------------------------------------------------------------------------------------
#' Get all the enhancer regions for the gene
#'
#' @rdname queryByRegion
#' @aliases queryByRegion
#'
#' @param obj An object of class GeneHancerDB
#' @param chrom character, the chromosome name WITHOUT leading 'chr'
#' @param start numeric the starting base
#' @param end numeric the ending base
#'
#' @export

setMethod('queryByRegion',  'GeneHancerDB',

     function(obj, chrom, start, end){
        query <- sprintf("select e.chr as chrom,
                          e.element_start as start,
                          e.element_end as end,
                          e.ghid as ghid,
                          e.type as class,
                          a.combined_score as combinedScore,
                          a.is_elite as elite,
                          a.symbol as gene,
                          a.eqtl_score as eqtl,
                          a.chic_score as HiC,
                          a.erna_score as erna,
                          a.expression_score as coexpression
                          from elements AS e,
                          associations AS a
                          where e.chr='%s'
                          AND e.element_start >= %d
                          AND e.element_end <= %d
                          AND a.ghid=e.ghid", chrom, start, end)

        db <- dbConnect(PostgreSQL(), user= "trena", password="trena",
                        dbname=GENEHANCER_VERSION, #"gh411",
                        host="khaleesi.systemsbiology.net", port="5432")
        tbl <- dbGetQuery(db, query)
        dbDisconnect(db)
        tbl$width <- with(tbl, 1 + end - start)
        tbl
        })

#------------------------------------------------------------------------------------------------------------------------
