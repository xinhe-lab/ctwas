
skip_if_no_LD_file <- function(LD_files) {
  if (any(!file.exists(LD_files))) {
    skip("LD files not available")
  }
}

skip_if_no_weight <- function(weight) {
  if (any(!file.exists(weight))) {
    skip("weight not available")
  }
}
