# Load every package used by the scripts in code/ (excluding any misc folder)
# and write the resulting sessionInfo() to ./sessionInfo.txt.
# Run from the project root: Rscript code/sessionInfo.R

scripts <- list.files("code", pattern = "\\.(R|Rmd)$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
scripts <- scripts[!grepl("(^|/)misc/", scripts, ignore.case = TRUE)]
scripts <- setdiff(scripts, "code/sessionInfo.R")

code_lines <- unlist(lapply(scripts, readLines, warn = FALSE))
code_lines <- code_lines[!grepl("^\\s*#", code_lines)]

extract <- function(pattern) {
  hits <- regmatches(code_lines, gregexpr(pattern, code_lines, perl = TRUE))
  unlist(lapply(hits, function(h) sub(pattern, "\\1", h, perl = TRUE)))
}

pkgs <- c(
  extract("(?:library|require|requireNamespace)\\(\\s*[\"']?([A-Za-z][A-Za-z0-9.]*)"),
  extract("\\b([A-Za-z][A-Za-z0-9.]*):::?")
)
pkgs <- sort(unique(pkgs))
pkgs <- setdiff(pkgs, rownames(installed.packages(priority = "base")))

loaded <- vapply(pkgs, function(p) {
  suppressWarnings(suppressPackageStartupMessages(
    requireNamespace(p, quietly = TRUE) && require(p, character.only = TRUE, quietly = TRUE)
  ))
}, logical(1))

out <- c(
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Scripts scanned:", length(scripts)),
  "",
  capture.output(sessionInfo())
)
if (any(!loaded)) {
  out <- c(out, "", "Packages referenced in code/ but not installed / failed to load:",
           paste(" ", names(loaded)[!loaded]))
}

writeLines(out, "sessionInfo.txt")
message("Wrote sessionInfo.txt (", sum(loaded), " packages loaded, ",
        sum(!loaded), " failed)")
