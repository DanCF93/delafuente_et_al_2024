#use this script to pool the data from multiple experiments into a single file

library(tidyverse)

#function to read all _Image.txt files and merge them together
read_txt_bind_rows <- function(path, pattern = "all_samples") {
  files = list.files(path, pattern, full.names = TRUE)
  lapply(files, read.delim) %>% bind_rows()
}


data <- read_txt_bind_rows(path = "../data/raw/")

write.table(data, "../data/tidy/cell_cycle_2d_48h_run_2_all_samples.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# make metadata from the data
metadata <- data %>%
  select(Metadata_Well:Metadata_WellRow) %>%
  distinct()

write.table(metadata, "../data/metadata/metadata.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
