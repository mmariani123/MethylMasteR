##test_that("multiplication works", {
##  expect_equal(2 * 2, 4)
##})

print(getwd())
input.dir <- "/data"

output.dir <- "/test"

sample.sheet.path <- paste0("/data",
                            .Platform$file.sep,
                            "Sample_Sheet_Test.csv")

hm450.anno.file.path <- paste0("/data",
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
