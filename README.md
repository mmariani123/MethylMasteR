# MethylMasteR
Michael Mariani PhD, Salas Lab, Dartmouth College 2022

#Install

require(devtools)\
devtools::install(mmariani123/MethylMasteR)

## Run MethylMasteR

## Set your paths:

input.dir <- "path.to.idat.dir"\
output.dir <- "path.to.output.dir"\
sample.sheet.path <- "path.to.sample.sheet"\
hm450.anno.file.path <- "path.to.hm450.manifest.hg38.rda"

## Select your routine :

1.) "sesame" for SeSAMEe\
2.) "hm450" for HM450\
3.) "champ" for ChAMP\
4.) "epicopy" for Epicopy\
5.) "compare" for comparison routine

## Select the other parameters and run!

methyl_master(routine              = "sesame",\
              input.dir            = input.dir,\
              output.dir           = output.dir,\
              sample.sheet.path    = sample.sheet.path,\
              r.lib.path           = .libPaths()[1],\
              file.sep             = "\\\\",\
              n.cores              = 1,
              os.type              = "windows",\
              proj                 = "TCGA-BLCA",\
              visualize            = TRUE,
              visualize.individual = FALSE, ##Only works for routine sesame\
              reference            = "internal", ##"comparison" or 'internal"\
              reference.name       = NA, ##NA for median, "all" For epicopy use,\
                                         ##see below notes\
              split.by             = NULL, ##NULL ##gender_reported\
                                     ##For splitting analysis across two factors\
                                     ##in a field, note MethylMaster only\
                                     ##supports splitting on a signal metadata\
                                     ##containing ONLY two factors currently\
                                     ##e.g. 'male' and 'female' in the\
                                     ##'gender_reported' field\
              comparison           = c("tumor","cord"),\
              overlap.density      = 0.1,\
              sesame.data.cache    = "EPIC",\
              sesame.data.normal   = 'EPIC.5.normal',\
              sesame.ref.version   = "hg38",\
              sesame.hm450.mean.correct = FALSE,
              sesame.form.thresholds = NULL, ##c(-0.3,0.3),\
              sesame.form.add.meta = c("Tumor"), ##NULL\
              hm450.workflow       = "B",
              hm450.anno.file.path = hm450.anno.file.path, ##Needed for hm450\
              champ.padj           = 0.05,\
              champ.control        = FALSE,\
              champ.run.combat     = FALSE,\
              champ.run.dmp        = FALSE, ##If only one pheno var must = FALSE\
              champ.run.dmr        = FALSE, ##If only one pheno var must = FALSE\
              champ.run.block      = FALSE, ##If only one pheno var must = FALSE\
              champ.run.gsea       = FALSE, ##Requires dmp and dmr results\
              champ.run.epimod     = FALSE, ##If only one pheno var must = FALSE\
              epi.run.gistic       = TRUE,\
              save.seg             = TRUE,\
              olaps.split.field    = "Sample_ID",\
              estimate.recurrence  = TRUE,\
              ov.less.stringent.ra.setting = TRUE,\
              ov.pvalue            = 0.05,\
              ov.keep.extra.columns = TRUE,\
              simplify.reduce      = weightedmean,\
              create.dir           = TRUE,\
              compare.list.files   = FALSE,\
              compare.files.in     = compare.files.in,\
              compare.names        = compare.names,\
              compare.olaps.size   = 1 ##Overlap of one or more base pairs\
                                       ##considered\
              )

