source("Module_A.R")
source("Module_B.R")

sCancer <- "LUAD"
sPath1 <- "./LUAD/"

#' Download biospecimen clinical data
path_clinical <-
  DownloadBiospecimenClinicalData(cancerType = sCancer,
                                  saveFolderName = paste(sPath1, "clinical", sep = ""))


#' Download somatic mutation data
path_somaticMutation <-
  DownloadSomaticMutationData(cancerType = sCancer,
                              assayPlatform = "somaticMutation_DNAseq",
                              saveFolderName = paste(sPath1, "somatic", sep = ""))

#' Download copy number alternation data
path_copyNumber <-
  DownloadCNAData(cancerType = sCancer,
                  assayPlatform = "cna_cnv.hg19",
                  saveFolderName = paste(sPath1, "copy_number", sep = ""))

#' Download DNA methylation 450 data
path_methylation_450 <-
  DownloadMethylationData(cancerType = sCancer,
                          assayPlatform = "methylation_450",
                          saveFolderName = paste(sPath1, "methylation", sep = ""))

#' Download miRNA	expression data
path_miRNAExp <-
  DownloadmiRNASeqData(cancerType = sCancer,
                       assayPlatform = "mir_HiSeq.hg19.mirbase20",
                       saveFolderName = paste(sPath1, "miRNA", sep = ""))

#' Download gene expression data
path_geneExp <-
  DownloadRNASeqData(cancerType = sCancer,
                     assayPlatform = "gene.normalized_RNAseq",
                     # assayPlatform = "gene_RNAseq",
                     saveFolderName = paste(sPath1, "gene_exp", sep = ""))

#' Download RPPA protein expression data
path_protein_RPPA <-
  DownloadRPPAData(cancerType = sCancer,
                   # assayPlatform = "protein_RPPA",
                   saveFolderName = paste(sPath1, "protein_RPPA", sep = ""))

#' Download iTRAQ protein expression data
# path_protein_iTRAQ <-
#   DownloadCPTACData(cancerType = sCancer,
#                     # assayPlatform = "proteome_iTRAQ",
#                     saveFolderName = paste(sPath1, "protein_iTRAQ/", sep = ""))


#' =============================================================================
#' Part 2: Perform basic processing of downloaded data using Module B functions
#' =============================================================================

#' Process somatic mutation data
list_somaticMutation <-
  ProcessSomaticMutationData(inputFilePath = path_somaticMutation[1],
                             outputFileName = paste(sCancer,
                                                    "somaticMutation",
                                                    sep = "__"),
                             outputFileFolder = paste(sPath1, "somatic/", "processed", sep = ""))

#' Process copy number alternation data
#' calculate an average copy number for each gene in each sample
list_copyNumber <-
  ProcessCNAData(inputFilePath = path_copyNumber[1],
                 outputFileName = paste(sCancer,
                                        "copyNumber",
                                        sep = "__"),
                 refGenomeFile = "./SupportingFiles/Hg19GenePosition.txt",
                 outputFileFolder = paste(sPath1, "copy_number/", "processed", sep = ""))

#' Process DNA methylation 450 data
list_methylation_450 <-
  ProcessMethylation450Data(inputFilePath = path_methylation_450[1],
                            outputFileName = paste(sCancer,
                                                   "methylation_450",
                                                   sep = "__"),
                            outputFileFolder = paste(sPath1, "methylation/", "processed", sep = ""))

#' Process miRNA expression data
list_miRNAExp <-
  ProcessmiRNASeqData(inputFilePath = path_miRNAExp[1],
                      outputFileName = paste(sCancer,
                                             "miRNAExp",
                                             sep = "__"),
                      outputFileFolder = paste(sPath1, "miRNA/", "processed", sep = ""))

#' Process gene expression data
list_geneExp <-
  ProcessRNASeqData(inputFilePath = path_geneExp[1],
                    outputFileName = paste(sCancer,
                                           "geneExp",
                                           sep = "__"),
                    dataType = "geneExp",
                    outputFileFolder = paste(sPath1, "gene_exp/", "processed", sep = ""))

#' Process RPPA protein expression data
list_protein_RPPA <-
  ProcessRPPADataWithGeneAnnotation(inputFilePath = path_protein_RPPA[1],
                                    outputFileName = paste(sCancer,
                                                           "protein_RPPA",
                                                           sep = "__"),
                                    outputFileFolder = paste(sPath1, "protein_RPPA/", "processed", sep = ""))

#' Process iTRAQ protein expression data
# list_protein_iTRAQ <-
#   ProcessCPTACData(inputFilePath = path_protein_iTRAQ[1],
#                    outputFileName = paste(sCancer,
#                                           "protein_iTRAQ",
#                                           sep = "__"),
#                    outputFileFolder = paste(sPath1, "protein_iTRAQ/", "processed", sep = ""))