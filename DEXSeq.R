library( "DEXSeq" )
#Samplesheet for DEXSeq analysis is present in the working directory
##Reading sample table and viewing. Including the rownames will helpe while viewing annotations
sampleTable <- read.csv("samplesheet.csv", header = T, sep = ",",row.names = c( "control1","test1","control2","test2","control3","test3"))
sampleTable
countFiles <- as.character(sampleTable$filePath)

#             name                                     filePath group condition
#control1 control1     ./S1_S10_L004.STAR.genome.Readcounts.txt  Rep1   control
#test1       test1     ./S2_S11_L004.STAR.genome.Readcounts.txt  Rep1  silenced
#control2 control2 ./Contr1_S15_L004.STAR.genome.Readcounts.txt  Rep2   control
#test2       test2  ./TRF22_S13_L004.STAR.genome.Readcounts.txt  Rep2  silenced
#control3 control3  ./Contr_S14_L004.STAR.genome.Readcounts.txt  Rep3   control
#test3       test3  ./TRF24_S12_L004.STAR.genome.Readcounts.txt  Rep3  silenced

#basename(countFiles)
#Storing my gff file
flattenedFile = list.files(pattern="gff$", full.names=TRUE)
#basename(flattenedFile)

##Structuring data and meta data
dxd = DEXSeqDataSetFromHTSeq(countFiles,sampleData=sampleTable,design= ~ sample + exon + condition:exon,flattenedfile=flattenedFile )

#Subsetting the analysis for TERF2. My csv only contains TERF2 ENSEMBLE ID because that's my gene of interest 
genesForSubset = read.table( 
  file.path("geneID.csv"), 
  stringsAsFactors=FALSE)[[1]]
dxd = dxd[geneIDs( dxd ) %in% genesForSubset,]

##Exploring the dataframe
#information regarding each column of the DEXSeqDataSet, is specified: colData(dxd)
#access the first 5 rows from the count data by doing: head( counts(dxd), 5 ). Note that the number of columns is 12, the first six (we have six samples) corresponding to the number of reads mapping to out exonic regions and the last six correspond to the sum of the counts mapping to the rest of the exons from the same gene on each sample.

#                         [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
#ENSG00000132604.11:E001  191  149  168  121  149  151 2240 1564 1790  1242  2323  1715
#ENSG00000132604.11:E002   50   45   38   36   47   44 2381 1668 1920  1327  2425  1822
#ENSG00000132604.11:E003  207  175  157  134  198  151 2224 1538 1801  1229  2274  1715
#ENSG00000132604.11:E004   95   60   79   47  108   69 2336 1653 1879  1316  2364  1797
#ENSG00000132604.11:E005    1    3    5    3    4    3 2430 1710 1953  1360  2468  1863

#access only the first five rows from the count belonging to the exonic regions (‘this’) (without showing the sum of counts from the rest of the exons from the same gene
#head( featureCounts(dxd), 5 )

#                         control1 test1 control2 test2 control3 test3
#ENSG00000132604.11:E001      191   149      168   121      149   151
#ENSG00000132604.11:E002       50    45       38    36       47    44
#ENSG00000132604.11:E003      207   175      157   134      198   151
#ENSG00000132604.11:E004       95    60       79    47      108    69
#ENSG00000132604.11:E005        1     3        5     3        4     3

#sampleAnnotation( dxd )
#head( featureCounts(dxd), 5 )

#Normalisation
dxd = estimateSizeFactors( dxd )

#Dispersion estimation: to estimate the variability of the data. This is necessary to be able to distinguish technical and biological variation (noise) from real effects on exon usage due to the different conditions.
#since we have single gene Dispersion  trend will not be captured buy the the estimateDispersion function 
dxd = estimateDispersions( dxd )
plotDispEsts( dxd )
#saving the dataset
save(dxd, file = file.path(".", "TERF2DEXSeqDataSet.RData"))

#Testing for differential exon usage
dxd = testForDEU( dxd )

#relative exon usage fold changes.
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition") #applying log2(x+ 1) instead of vst because Dispersion function was not parameric
dxr1 = DEXSeqResults( dxd )
#dxr1
#mcols(dxr1)$description

#Check how many exonic regions are significant with a false discovery rate of 10%:
table ( dxr1$padj < 0.1 )
#FALSE 
#34 
#At 20% i.e 0.2 value
#FALSE  TRUE 
#26     8 

# To check how many genes are affected. In our case since we have one gene this isn't required
#table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )
plotMA( dxr1, cex=0.8 )
plotMA( dxr1, cex=0.8, alpha = 0.2 )# adjusted p-values

##Estimating Batch Effect on Exon Expression
formulaFullModel    =  ~ sample + exon + group:exon + condition:exon
formulaReducedModel =  ~ sample + exon + group:exon

##
dxd = estimateDispersions( dxd, formula = formulaFullModel )
#you had estimated dispersions, replacing these
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.
#Warning message:
#  In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
#              Estimated rdf < 1.0; not estimating variance
#Seems, the software cannot find parameters a and b to fit the dispersion over mean trend using that formula. The better the prediction, the less error there is in that prediction. A PRE statistic takes values between 0 and 1. 0 means no reduction in error, 1 means that there is perfect prediction—the error is completely eliminated


dxd = testForDEU( dxd, reducedModel = formulaReducedModel, fullModel = formulaFullModel )
dxr2 = DEXSeqResults( dxd )
table( dxr2$padj < 0.1 )
table( before = dxr1$padj < 0.1, now = dxr2$padj < 0.1 )
table( before = dxr1$padj < 0.2, now = dxr2$padj < 0.2 )


#Visualization
plotDEXSeq( dxr2, "ENSG00000132604.11", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
## View Transcripts. DEXSeq is designed to find changes in relative exon usage, i.e., changes in the expression of individual exons that are not simply the consequence of overall up- or down-regulation of the gene.
library(grid, lib.loc = "C:/Program Files/R/R-3.6.2/library")
png("test.png", width = 10, height = 9, units = 'in', res = 700)
#After correction for Passages
plotDEXSeq( dxr2, "ENSG00000132604.11", displayTranscripts=TRUE, expression=FALSE, splicing=TRUE,legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

#DEXSeqHTML( dxr2, FDR=0.1, path=".",file="testForDEU.html", color=c("#FF000080", "#0000FF80") )

#Before correction for passages
plotDEXSeq( dxr1, "ENSG00000132604.11", expression=FALSE, splicing=TRUE,legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2,FDR=0.2 )



