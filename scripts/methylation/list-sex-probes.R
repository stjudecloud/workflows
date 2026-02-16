library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
sex_probes <- rownames(ann)[ann$chr == "chrX" | ann$chr == "chrY"]

write.table(
  sex_probes,
  file = "sex_probes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
