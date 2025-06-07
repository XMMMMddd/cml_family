# R script to source all Rcpp files in a specified directory

# 1. Specify the directory containing the Rcpp files
rcpp_files_directory <- "cml家庭函数/cml家庭函数overlap版本/overlap2的rcpp代码/"

# Check if the directory exists
if (!dir.exists(rcpp_files_directory)) {
  stop(paste("Directory not found:", rcpp_files_directory))
}

# 2. Get a list of all files in the directory
all_files <- list.files(path = rcpp_files_directory, full.names = TRUE)

# 3. Filter for files ending with .cpp
rcpp_file_paths <- all_files[grep("\\.cpp$", all_files, ignore.case = TRUE)]

# 4. Check if any .cpp files were found
if (length(rcpp_file_paths) == 0) {
  message(paste("No .cpp files found in directory:", rcpp_files_directory))
} else {
  message(paste("Found", length(rcpp_file_paths), ".cpp files to source:"))
  
  # Loop through the .cpp files and source them
  for (cpp_file in rcpp_file_paths) {
    message(paste("Sourcing:", cpp_file))
    tryCatch({
      Rcpp::sourceCpp(cpp_file)
      message(paste("Successfully sourced:", cpp_file))
    }, error = function(e) {
      warning(paste("Error sourcing", cpp_file, ":", e$message))
    })
  }
  message("Finished sourcing all .cpp files.")
}

# Clean up variables if desired
# rm(rcpp_files_directory, all_files, rcpp_file_paths, cpp_file)