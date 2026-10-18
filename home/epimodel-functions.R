# Writes home/epimodel-functions.tsv: every function EpiModel exports and the
# reference page that documents it. home/build-data.lua uses the list to find
# EpiModel calls in each example, since the site is rendered without R.
# Rerun from the project root after an EpiModel release:
#   Rscript home/epimodel-functions.R

db <- tools::Rd_db("EpiModel")
aliases <- do.call(rbind, lapply(names(db), function(f) {
  a <- tools:::.Rd_get_metadata(db[[f]], "alias")
  if (length(a)) data.frame(alias = a, topic = sub(".Rd", "", f, fixed = TRUE))
}))
fx <- sort(getNamespaceExports("EpiModel"))
topic <- aliases$topic[match(fx, aliases$alias)]
writeLines(
  c(paste0("# EpiModel ", packageVersion("EpiModel"), " exports and reference topics"),
    paste(fx, ifelse(is.na(topic), "", topic), sep = "\t")),
  "home/epimodel-functions.tsv"
)
