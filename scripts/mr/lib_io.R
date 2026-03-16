read_tsv_any <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  utils::read.delim(
    con,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character"
  )
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

timestamp_id <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}
