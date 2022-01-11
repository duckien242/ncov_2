# Drops "County" from Location field in GISAID metadata file so counties aren't represented 2x in the final build
# Example: some sequences list "New Haven" as the Location, and others "New Haven County"

library(tidyverse)
path <- getwd()
filename <- "metadata_nextstrain.tsv"
dat <- read_tsv(file = paste0(path, "/pre-analyses/", filename))
stopwords = "County"
dat$location <- gsub(paste0(stopwords,collapse = "|"),"", dat$location)
write_tsv(dat, file = paste0(path, "/pre-analyses/", filename))

