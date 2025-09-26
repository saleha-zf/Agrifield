options(warn = -1)

library(plyr)
library(tidyverse)
library(raster)
library(celestial)
library(caret)
library(fastICA)
library(SOAR)
library(RStoolbox)
library(jsonlite)
library(data.table)
library(spdep)
library(FIELDimageR)

cat("=== STARTING TEST DATA GENERATION ===\n")

# Paths
output_path <- Sys.getenv('OUTPUT_DATA')
path.dir <- Sys.getenv('INPUT_DATA')
data.dir <- file.path(path.dir, "ref_agrifieldnet_competition_v1")
test_path <- "ref_agrifieldnet_competition_v1_labels_test"
image_path <- "ref_agrifieldnet_competition_v1_source"

# Helper function to extract IDs from URLs
map.func <- function(x, y = 7){
  map <- x %>%
    sapply(FUN = function(x){strsplit(x, '[,.:_/"]')[[1]][y]}) %>%
    sub("[[:punct:]]", '', .) %>%
    sub("'", '', .)
  return(map)
}

####### Load Test Collection JSON
test_coll_path <- file.path(data.dir, test_path, "collection.json")
cat("Reading test collection JSON from:\n", test_coll_path, "\n")
if(!file.exists(test_coll_path)) stop("ERROR: Test collection JSON not found!")

colls <- jsonlite::read_json(test_coll_path)
qq <- colls$links[3:(length(colls$links)-1)]
cat("Number of test links:", length(qq), "\n")

# Map each field to its actual numbered folder (_0000, _0001, ...)
ids <- c()
field_path <- c()
for(i in 1:length(qq)){
  id <- sprintf("%04d", i-1)  # zero-padded folder index
  folder_name <- paste0(test_path, "_", id)
  tif_path <- file.path(data.dir, test_path, folder_name, "field_ids.tif")
  ids <- c(ids, id)
  field_path <- c(field_path, tif_path)
}
test_fields <- data.frame(ids = ids, field_path = field_path, stringsAsFactors = FALSE)
cat("Number of test fields parsed:", nrow(test_fields), "\n")
if(nrow(test_fields) == 0) stop("No test fields found! Check JSON paths.")

##########################################
#### GENERATE RAW TEST DATA
##########################################

# test <- data.frame()

# # Updated bands matching actual filenames
# bands <- c(
#   B01 = "B01.B1.tif",
#   B02 = "B02.B2.tif",
#   B03 = "B03.B3.tif",
#   B04 = "B04.B4.tif",
#   B05 = "B05.B5.tif",
#   B06 = "B06.B6.tif",
#   B07 = "B07.B7.tif",
#   B08 = "B08.B8.tif",
#   B09 = "B09.B9.tif",
#   B11 = "B11.B11.tif",
#   B12 = "B12.B12.tif",
#   B8A = "B8A.B8A.tif"
# )


# for(i in 1:nrow(test_fields)){
#   cat("Processing test field:", i, "/", nrow(test_fields), "ID:", test_fields$ids[i], "\n")
  
#   # Load field raster
#   if(!file.exists(test_fields$field_path[i])){
#     cat("WARNING: Field raster does not exist:", test_fields$field_path[i], "\n")
#     next
#   }
  
#   f_mat <- raster(test_fields$field_path[i]) %>% as.matrix()
#   f_mat_rs <- which(rowSums(f_mat, na.rm = TRUE) > 0)
#   f_mat_cs <- which(colSums(f_mat, na.rm = TRUE) > 0)
  
#   if(length(f_mat_rs) == 0 | length(f_mat_cs) == 0){
#     cat("WARNING: Empty field matrix for field", test_fields$ids[i], "\n")
#     next
#   }
  
#   f_mat <- f_mat[f_mat_rs, f_mat_cs]
#   train_data <- data.frame()
# }
# for(b in names(bands)){
#   img_file <- file.path(data.dir, image_path, paste0(image_path, "_", test_fields$ids[i]), bands[b])
  
#   if(!file.exists(img_file)){
#     cat("WARNING: Band file missing:", img_file, "\n")
#     train_data[[b]] <- NA  # fill missing band with NA
#     next
#   }
  
#   mm <- raster(img_file) %>% as.matrix()
#   mm <- mm[f_mat_rs, f_mat_cs]
  
#   train_data[[b]] <- as.vector(mm)  # column name is exactly 'B01', 'B02', etc.
# }
#   cat("Rows in train_data for this field before filtering NAs:", nrow(train_data), "\n")
#   # Initialize final test data frame
test <- data.frame()
bands <- c(paste0("B0", 1:9), "B11", "B12")  # Band names exactly as expected

for(i in 1:nrow(test_fields)) {
  cat("Processing test field:", i, "/", nrow(test_fields), "ID:", test_fields$ids[i], "\n")
  
  # Load field raster
  if(!file.exists(test_fields$field_path[i])) {
    cat("WARNING: Field raster does not exist:", test_fields$field_path[i], "\n")
    next
  }
  
  f_mat <- raster(test_fields$field_path[i]) %>% as.matrix()
  f_mat_rs <- which(rowSums(f_mat, na.rm = TRUE) > 0)
  f_mat_cs <- which(colSums(f_mat, na.rm = TRUE) > 0)
  
  if(length(f_mat_rs) == 0 | length(f_mat_cs) == 0) {
    cat("WARNING: Empty field matrix for field", test_fields$ids[i], "\n")
    next
  }
  
  f_mat <- f_mat[f_mat_rs, f_mat_cs]
  n_pix <- length(f_mat)  # number of pixels in the field

  # Initialize train_data with first band to get correct number of rows
  first_band_file <- file.path(data.dir, image_path, paste0(image_path, "_", test_fields$ids[i]), "B01.tif")
  if(!file.exists(first_band_file)) {
    cat("WARNING: First band missing:", first_band_file, "\n")
    next
  }
  mm <- raster(first_band_file) %>% as.matrix()
  mm <- mm[f_mat_rs, f_mat_cs]
  train_data <- data.frame(B01 = as.vector(mm))

  # Read remaining bands
  for(b in bands[-1]) {
    img_file <- file.path(data.dir, image_path, paste0(image_path, "_", test_fields$ids[i]), paste0(b, ".tif"))
    if(!file.exists(img_file)) {
      cat("WARNING: Band file missing:", img_file, "\n")
      train_data[[b]] <- NA  # fill missing band with NA
      next
    }
    mm <- raster(img_file) %>% as.matrix()
    mm <- mm[f_mat_rs, f_mat_cs]
    train_data[[b]] <- as.vector(mm)
  }
  
  # Flatten field mask as fid
  fid <- as.vector(f_mat)

  # Combine into final data frame for this field
  dd <- data.frame(folder = test_fields$ids[i], fid = fid, train_data) %>% dplyr::filter(!is.na(fid))
  cat("Rows for this field:", nrow(dd), "\n")
  
  # Append to cumulative test data
  test <- dplyr::bind_rows(test, dd)
  cat("Cumulative test rows:", nrow(test), "\n")
  
  # Clean up
  rm(f_mat, fid, dd, train_data, mm, f_mat_rs, f_mat_cs); invisible(gc())
}


cat("=== RAW TEST DATA MATRIX CREATED ===\n")
cat("Total test rows:", nrow(test), "Columns:", ncol(test), "\n")
if(nrow(test) == 0) stop("ERROR: No test data generated! Cannot continue.")

##########################################
#### FEATURE ENGINEERING
##########################################
test <- test %>%
  mutate(
    ndvi =(B08 - B04)/ (B08 + B04),
    GLI = (2*B03-B04-B02)/(2*B03+B04+B02),
    CVI = (B08 / B03) * (B04 / B03),
    SIPI = (B08 - B02) / (B08 - B04),
    S2REP = 705 + 35 * ((((B07 + B04)/2) - B05)/(B06 - B05)),
    CCCI = ((B08 - B05) / (B08 + B05)) / ((B08 - B04) / (B08 + B04)),
    hue = atan(2*(B02-B03-B04)/30.5*(B03-B04)),
    RENDVI = (B06 - B05) / (B06 + B05), 
    RECI = (B08 / B04)-1,
    RECI2 = (B08 / B05)-1,
    evi = 2.5 * (B08 - B04) / ((B08 + 6.0 * B04 - 7.5 * B02) + 1.0),
    evi2 = 2.4 * (B08 - B04) / (B08 + B04 + 1.0),
    npcri = (B04 - B02) / (B04 + B02),
    ndwi = (B03 - B08) / (B03 + B08)
  )

cat("Feature engineering completed. Dimensions:", dim(test), "\n")

##########################################
#### AGGREGATE PER FIELD
##########################################
df_test <- test %>% filter(!is.na(fid)) %>%
  group_by(fid) %>%
  mutate(field_tile_count = n(), field_overlap_count = length(unique(folder))) %>%
  ungroup() %>%
  dplyr::select(-c(folder, id)) %>%
  group_by(fid) %>%
  summarise_all(list(median = median, max = max)) %>%
  ungroup()

cat("Aggregated df_test dimensions:", dim(df_test), "\n")

##########################################
#### JOIN FIELD DETAILS (TILE SIZE)
##########################################
test_details <- data.frame()
for(i in 1:nrow(test_fields)){
  folder_name <- paste0(test_path, "_", test_fields$ids[i])
  stac_file <- file.path(data.dir, test_path, folder_name, "stac.json")
  if(!file.exists(stac_file)){
    cat("WARNING: stac.json missing:", stac_file, "\n")
    next
  }
  stac <- jsonlite::read_json(stac_file)
  bbox <- unlist(stac$bbox)
  tile_width <- bbox[3] - bbox[1]
  tile_height <- bbox[4] - bbox[2]
  dd <- data.frame(fid = test_fields$ids[i], tile_width = tile_width, tile_height = tile_height)
  test_details <- rbind(test_details, dd)
}

cat("Test details loaded. Rows:", nrow(test_details), "\n")

df_test <- df_test %>%
  left_join(test_details, by = "fid") %>%
  mutate(field_tile_size = 20000*field_tile_count_median*tile_width*tile_height)

cat("Final df_test dimensions after joining details:", dim(df_test), "\n")

##########################################
#### WRITE TO CSV
##########################################
output_file <- file.path(output_path, "Final_Test.csv")
fwrite(df_test, file = output_file, row.names = FALSE)
cat("Final_Test.csv written successfully to:", output_file, "\n")
