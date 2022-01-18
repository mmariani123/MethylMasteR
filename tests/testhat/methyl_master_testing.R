#!/usr/bin/env Rscript

##debug(sesame::openSesame)
##debug(methyl_master_hm450)

##debug(methyl_master)
##debug(methyl_master_epicopy)
##debug(epicopy.mm)
##debug(Epicopy::export_gistic)
##debug(getLRR.mm)
##debug(LRRtoCNA.mm)

############################ debug sesame ##################################
##debug(methyl_master)
##debug(methyl_master_sesame)
##undebug(AutoCorrectPeak.mm)
##undebug(methyl_master_formatting_sesame)
##debug(methyl_master_olaps_and_visualize)
##debug(methyl_master_plot_individual)
##debug(visualizeSegments.mm)
##debug(methyl_master_population_ranges)
##debug(methyl_master_plot_individual)

############################## debug hm450 #################################
##debug(methyl_master)
##debug(methyl_master_hm450)
##debug(findSegments2) ##Need to replace start
##and stop with bp instead of probeID
##debug(methyl_master_formatting_hm450)
##debug(methyl_master_olaps_and_visualize)

############################## debug champ #################################
##debug(methyl_master)
##debug(methyl_master_champ)
##debug(ChAMP::champ.process)
##debug(methyl_master_formatting_champ)
##debug(methyl_master_olaps_and_visualize)

############################## debug epicopy ###############################
##debug(methyl_master)
##debug(methyl_master_epicopy)
##debug(Epicopy::epicopy)
##debug(readSheet.mm)
##debug(read.metharray.exp)
##debug(methyl_master_formatting_epicopy)
##debug(methyl_master_epicopy_helper_functions)
##debug(minfi::read.metharray.sheet)
##debug(minfi::read.metharray.exp)
##debug(minfi::read.metharray)
##debug(.funnorm.mm)
##debug(getLRR.mm)
##
##undebug(methyl_master_epicopy)
##undebug(Epicopy::epicopy)
##undebug(minfi::read.metharray.sheet)
##undebug(minfi::read.metharray.exp)
##undebug(minfi::read.metharray)
##undebug(.funnorm.mm)
##undebug(getLRR.mm)

############################## debug compare ###############################
##debug(methyl_master_compare)

############################## BLCA ########################################

##input.dir <- paste0("G:\\My Drive\\dartmouth\\salas_lab",
##                    "\\working\\cnv\\cnv_testing",
##                    "\\sesame_combined_results")

##input.dir <- paste0("G:\\My Drive\\dartmouth",
##                         "\\salas_lab\\working\\cnv",
##                         "\\cnv_testing",
##                         "\\data\\blca_idat_files_pooled")

##output.dir <- paste0("G:\\My Drive\\dartmouth",
##                     "\\salas_lab\\working\\cnv",
##                     "\\cnv_testing",
##                     "\\sesame_compare_results")

##sample.sheet.path <- paste0("G:\\My Drive\\dartmouth",
##                            "\\salas_lab\\working\\cnv",
##                            "\\cnv_testing",
##                            "\\Sample_Sheet.csv")

##hm450.anno.file.path <- paste0("G:\\My Drive\\dartmouth",
##                               "\\salas_lab\\working\\cnv",
##                               "\\cnv_testing",
##                               "\\hm450.manifest.hg38.rda")

########################### KIRC 3p ##################################

input.dir <- "data"

output.dir <- "test"

sample.sheet.path <- paste0("data",
                            .Platform$file.sep,
                            Sample_Sheet_Test.csv)

hm450.anno.file.path <- paste0("data",
                               .Platform$file.sep,
                               "hm450.manifest.hg38.rda")

methyl_master(input.dir            = input.dir,
              output.dir           = output.dir,
              sample.sheet.path    = sample.sheet.path,
              r.lib.path           = .libPaths()[1],
              file.sep             = "\\\\",
              n.cores              = 1,
              os.type              = "windows",
              proj                 = "TCGA-BLCA",
              visualize            = TRUE,
              visualize.individual = FALSE, ##Only works for routine sesame
              routine              = "sesame",
              reference            = "internal", ##"comparison" or 'internal"
              reference.name       = NA, ##NA for median, "all" For epicopy use,
                                         ##see below notes
                                     ##For splitting analysis across two factors
                                     ##in a field, note MethylMaster only
                                     ##supports splitting on a signal metadata
                                     ##containing ONLY two factors currently
                                     ##e.g. 'male' and 'female' in the
                                     ##'gender_reported' field
              comparison           = c("tumor","normal"),
              overlap.density      = 0.1,
              sesame.data.cache    = "EPIC",
              sesame.data.normal   = 'EPIC.5.normal',
              sesame.ref.version   = "hg38",
              sesame.hm450.mean.correct = FALSE,
              sesame.form.thresholds = NULL, ##c(-0.3,0.3),
              sesame.form.add.meta = NULL, ##c("Tumor"), ##NULL
              hm450.workflow       = "B",
              hm450.anno.file.path = hm450.anno.file.path, ##Needed for hm450
              champ.padj           = 0.05,
              champ.control        = FALSE,
              champ.run.combat     = FALSE,
              champ.run.dmp        = FALSE, ##If only one pheno var must = FALSE
              champ.run.dmr        = FALSE, ##If only one pheno var must = FALSE
              champ.run.block      = FALSE, ##If only one pheno var must = FALSE
              champ.run.gsea       = FALSE, ##Requires dmp and dmr results
              champ.run.epimod     = FALSE, ##If only one pheno var must = FALSE
              epi.run.gistic       = TRUE,
              save.seg             = TRUE,
              olaps.split.field    = "Sample_ID",
              estimate.recurrence  = TRUE,
              ov.less.stringent.ra.setting = TRUE,
              ov.pvalue            = 0.05,
              ov.keep.extra.columns = TRUE,
              simplify.reduce      = weightedmean,
              create.dir           = TRUE,
              compare.list.files   = FALSE,
              compare.files.in     = compare.files.in,
              compare.names        = compare.names,
              compare.olaps.size   = 1 ##Overlap of one or more base pairs
                                       ##considered
              )

########################### Create a sample sheet #############################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

####First get primary ids from sample names from cbioportal
##
##clin.tcga.file.path <- paste0("C:\\Users\\User\\Desktop",
##                      "\\select_kidney_samples_cbioportal_vhl_01042021.xlsx")
##
##clin.tcga.file.path <- paste0("C:\\Users\\User",
##                              "\\Desktop\\clinical_info.csv")
##
##clin.sub.dir <- paste0("G:\\My Drive\\dartmouth\\salas_lab",
##                       "\\working\\cnv\\kirc_cbio_3p_select")
##
##methyl_master_tcga_clin_data()
##
####Then create a sample sheet with the primary ids
##
##output.path.name <- paste0("C:\\Users\\User\\Desktop",
##                           "\\Sample_Sheet_kirc_cbio_3p.csv")
##
##file.sep <- "\\"
##
##idat.dir <- paste0("G:\\My Drive\\dartmouth",
##                   "\\salas_lab\\working\\cnv",
##                   "\\kirc_idat_files_pooled")
##
##sample.sheet.path <-
##  "C:\\Users\\User\\Desktop\\Sample_Sheet_testing.csv"
##
##Sample_Group <- c("tumor")
##
##methyl_master_create_sample_sheet()
##
############################# Compare testing ###################################
##
##work.dir <- "G:\\My Drive\\dartmouth\\salas_lab\\working\\cnv\\cnv_testing"
##
##time.1  <- paste0(work.dir,
##                  "\\",
##          "test_sesame_internal\\time_mem.sesame.2021.12.31.12.57.13.txt")
##time.2  <- paste0(work.dir,
##                  "\\",
##          "test_sesame_internal\\sesame_tumor_overlaps_filt.csv")
##time.3  <- paste0(work.dir,
##                  "\\",
##          "test_hm450_internal\\time_mem.hm450.2021.12.29.18.20.28.txt")
##time.4  <- paste0(work.dir,
##                  "\\",
##          "test_hm450_internal\\hm450_tumor_overlaps_filt.csv")
##olaps.1 <- paste0(work.dir,
##                  "\\",
##          "test_champ_internal\\time_mem.champ.2021.12.29.12.08.41.txt")
##olaps.2 <- paste0(work.dir,
##                  "\\",
##                  "test_champ_internal\\champ_tumor_overlaps_filt.csv")
##olaps.3 <- paste0(work.dir,
##                  "\\",
##          "test_epicopy_internal\\time_mem.epicopy.2021.12.29.14.22.02.txt")
##olaps.4 <- paste0(work.dir,
##                  "\\",
##          "test_epicopy_internal\\epicopy_tumor_overlaps_filt.csv")
##
##compare.files.in <- c(time.1,
##                      time.2,
##                      time.3,
##                      time.4,
##                      olaps.1,
##                      olaps.2,
##                      olaps.3,
##                      olaps.4)
##
##compare.names <- c("sesame_internal",
##                   "hm450_internal",
##                   "champ_internal",
##                   "epicopy_internal"
##                   )

###############################################################################

##Epicopy specifying normals smaples:
##From https://github.com/sean-cho/Epicopy/blob/master/vignettes/Epicopy.Rmd
## Specifying normal samples

### For epicopy function

#### User's input
##If reference normal samples are included in the raw data files,
##a column specifying the normal status of the samples should be included.
##Normals have to be tagged using the character string normal
##(case insensitive).

#### EpicopyData/TCGA-derived normals
##Otherwise, users can specify the one of the three type of normals
##included in the `EpicopyData` package, derived from normal
##solid tissue arrayed by the Cancer Genome Atlas (TCGA).
##To use those, users may input one of four arguments
##'thyroid', 'breast', 'lung', or 'all', the last of which
##'##uses all available normals. The default uses all normal samples.

### For getLRR function

#### Default
##Defaults to `NULL` which uses all the normal samples
##included with the `EpicopyData` package.

#### EpicopyData/TCGA-derived normals
##To use `EpicopyData` included normals, as before,
##normals can be specified using one of four arguments.

#### User's input
##If the user has their own normal samples,
##they can specify either a numeric/integer
##index that identifies the positions of the
##normal samples in the `RGChannelSet` or a
##logical vector that flags normal samples as `TRUE`.

#### Use all samples
##The argument `Normals = NA` uses the mode/median (as specified by the user)
##of all the samples, regardless of status, as reference normals.
##The idea behind this is that the median copy number of a given genomic region
##of all the samples should center around zero.
##Recommended only when there are many samples in the array.
