
skip_if_no_LD_file <- function(LD_matrix_files) {
  if (any(!file.exists(LD_matrix_files))) {
    skip("LD matrix files not available")
  }
}

skip_if_no_weight <- function(weight) {
  if (any(!file.exists(weight))) {
    skip("weight not available")
  }
}
