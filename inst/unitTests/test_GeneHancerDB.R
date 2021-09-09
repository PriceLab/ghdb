# test_GeneHancerDB.R
#------------------------------------------------------------------------------------------------------------------------
library(ghdb)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("ghdb")){
   ghdb <- GeneHancerDB()
   } # creating for several tests below

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_listTissues()
   test_.eliminateDupsCollapseTissues()
   test_foxo6()
   test_failureByTissue()
   test_queryUnknownGene()

   test_getEnhancers()

   test_queryByRegion()
   test_to.hg19()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))
   checkTrue("GeneHancerDB" %in% is(ghdb))

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
test_listTissues <- function()
{
   message(sprintf("--- test_listTissues"))
   tissues <- listTissues(ghdb, "MEF2C")
   checkTrue(length(tissues) > 100)
   checkTrue(all(c("placenta", "Placenta") %in% tissues))

   tissues.gata2 <- listTissues(ghdb, targetGene="GATA2")
     # 120-130 with gh 411, 152 with v50
   checkTrue(length(tissues.gata2) > 150)
   checkTrue(length(tissues.gata2) < 160)

} # test_listTissues
#------------------------------------------------------------------------------------------------------------------------
test_.eliminateDupsCollapseTissues <- function()
{
   message(sprintf("--- test_.eliminateDupsCollapseTissues"))

   load(system.file(package="TrenaProjectHG38", "extdata", "testing", "tbl.gh.results.toTestReductionPreservingTissues.RData"))
   checkEquals(dim(tbl), c(956, 17))
   tbl.trimmed <- ghdb:::.eliminateDupsCollapseTissues(tbl)
   checkEquals(dim(tbl.trimmed), c(50, 16))

      # do a rough assay of the collapsed tissue column
   spread <- fivenum(nchar(tbl.trimmed$tissue))  # [1]   16   62  156  338 2213
   checkTrue(spread[1] < 20)
   checkTrue(spread[5] > 2000)

} # test_.eliminateDupsCollapseTissues
#------------------------------------------------------------------------------------------------------------------------
# a simple test case, chosen for historical (and no longer important) reasons
test_foxo6 <- function()
{
   message(sprintf("--- test_foxo6"))

      # ENCODE and Ensembl use different capitalization
   tissues <-  c("placenta", "Placenta")
   tbl <- retrieveEnhancersFromDatabase(ghdb, "FOXO6", tissues)
   checkTrue(is.data.frame(tbl))
   checkEquals(nrow(tbl), 19)   # 19 for v54, 23 for v50, 13 for v411
   checkTrue(all(tissues %in% tbl$tissue))
   checkEquals(length(which(duplicated(tbl$sig))), 0)  # no duplicated regions
   checkTrue(length(grep(";", tbl$tissue)) > 0)        # 4 instances of "Placenta;placenta"

   tbl.all <- retrieveEnhancersFromDatabase(ghdb, "FOXO6", tissues="all")
   checkTrue(is.data.frame(tbl.all))
   checkEquals(nrow(tbl.all), 70)   # 43 for v411, 71 for v50, 70 for v54
   checkEquals(length(which(duplicated(tbl.all$sig))), 0)  # no duplicated regions
   checkTrue(max(nchar(tbl.all$tissue)) > 1200)            # 2151 for v54
   checkTrue(length(grep(";", tbl.all$tissue)) > 40)       # most tissue values are semicolo separated multiples

} # test_gata6
#------------------------------------------------------------------------------------------------------------------------
# a simple test case, chosen for historical (and no longer important) reasons
test_failureByTissue <- function()
{
   message(sprintf("--- test_failureByTissue"))

   trem2.tissues <- retrieveEnhancersFromDatabase(ghdb, "TREM2", tissues="all")$tissue
   tissues <-  c("brain", "Brain")
   checkTrue(length(intersect(tissues, trem2.tissues)) == 0)

   suppressWarnings(
      tbl <- retrieveEnhancersFromDatabase(ghdb, "TREM2", tissues)
      )
   checkEquals(nrow(tbl), 0)

      # double check now that at least one tissue reports for TREM2
   tbl <- retrieveEnhancersFromDatabase(ghdb, "TREM2", "placenta")
   checkTrue(nrow(tbl) > 0)

} # test_failureByTissue
#------------------------------------------------------------------------------------------------------------------------
test_queryUnknownGene <- function()
{
   message(sprintf("--- test_queryUnknownGene"))
   symbol <- "LINC00982"
   suppressWarnings(
      tbl <- retrieveEnhancersFromDatabase(ghdb, symbol, tissues="all")
      )
   checkEquals(dim(tbl), c(0,0))

} # test_queryUnknownGene
#------------------------------------------------------------------------------------------------------------------------
test_getEnhancers <- function()
{
   message(sprintf("--- test_getEnhancers"))

   all.tissues <- getEnhancerTissues(ghdb, "TREM2")
   tbl.trem2.all <- getEnhancers(ghdb, "TREM2")
   checkTrue(all(tbl.trem2.all$geneSymbol == "TREM2"))
   dim(tbl.trem2.all)
   checkTrue(nrow(tbl.trem2.all) > 8)   # gh50: 16   gh54: 16

   tbl.mef2c <- getEnhancers(ghdb, "MEF2C")
   checkTrue(all(tbl.mef2c$geneSymbol == "MEF2C"))
   checkTrue(nrow(tbl.mef2c) > 190)  # gh411: 12 , gh50: 201, gh54 206

   suppressWarnings(tbl.bogus <- getEnhancers(ghdb, "bogus99"))
   checkEquals(nrow(tbl.bogus), 0)

      #----------------------------------------------------------------------------------
      # hoxa5 has a nonsensical promoter, 86kb in size
      # simon fishilevich calls this an outlier, worth elminating
      # comes from ENCODE.  he suggests this filter:
      #    drop all elements >=10kb AND (overlap a gene TSS OR overlap >=3 gene exons)
      # i am implementing just the first clause for now
      #----------------------------------------------------------------------------------

   tbl.hoxa5 <- getEnhancers(ghdb, "HOXA5", maxSize=100000)
   checkEquals(nrow(tbl.hoxa5), 11)
   sizes <- with(tbl.hoxa5, 1 + end - start)
   checkEquals(length(which(sizes > 10000)), 2)  # gh411: 1   gh50: 2, gh54: 2

   tbl.hoxa5 <- getEnhancers(ghdb, "HOXA5", maxSize=10000)
   checkEquals(nrow(tbl.hoxa5), 9)   # 10 with gh411
   sizes <- with(tbl.hoxa5, 1 + end - start)
   checkEquals(length(which(sizes > 10000)), 0)

} # test_getEnhancers
#------------------------------------------------------------------------------------------------------------------------
test_queryByRegion <- function()
{
    message(sprintf("--- test_queryByRegion"))

    chrom <- 'chr1'
    start <- 173595097
    end   <- 173692966

    tbl <- queryByRegion(ghdb, chrom, start, end)
    checkTrue(nrow(tbl) > 50)   # gh54: 76
    checkTrue(all(tbl$chrom == chrom))
    checkTrue(all(tbl$start >= start))
    checkTrue(all(tbl$end <= end))

    tbl.best <- subset(tbl, combinedscore > 250)
    checkTrue(nrow(tbl.best) > 10)   # 15 with gh version 5.0 > 500, 19 >250 in gh54
    checkTrue(nrow(tbl.best) < 20)

} # test_queryByRegion
#------------------------------------------------------------------------------------------------------------------------
test_to.hg19 <- function()
{
   message(sprintf("--- test_to.hg19"))

   tbl.all <- retrieveEnhancersFromDatabase(ghdb, "FOXO6", tissues="all")
   checkEquals(nrow(tbl.all), 70)   # gh50=71, gh54=70

   tbl.hg19 <- to.hg19(ghdb, tbl.all)
   checkEquals(dim(tbl.all), dim(tbl.hg19))
   checkTrue(length(intersect(tbl.all$start, tbl.hg19$start)) == 0)
   checkTrue(all(tbl.all$chrom == tbl.hg19$chrom))

} # test_to.hg19
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
