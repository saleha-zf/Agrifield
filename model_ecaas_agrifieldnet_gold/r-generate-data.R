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

# Helper function
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
cat("[INFO] Number of test links:", length(qq), "\n")

# Map each tile to its numbered folder (_0000, _0001, ...)
ids <- c(); field_path <- c()
for(i in 1:length(qq)){
  id <- sprintf("%04d", i-1)
  folder_name <- paste0(test_path, "_", id)
  tif_path <- file.path(data.dir, test_path, folder_name, "field_ids.tif")
  ids <- c(ids, id)
  field_path <- c(field_path, tif_path)
}
test_fields <- data.frame(ids = ids, field_path = field_path, stringsAsFactors = FALSE)
cat("[INFO] Number of test fields parsed:", nrow(test_fields), "\n")
if(nrow(test_fields) == 0) stop("No test fields found! Check JSON paths.")

#### QUICK SKIP if Final_Test.csv exists
output_file <- file.path(output_path, "Final_Test.csv")
if(file.exists(output_file)){
  existing <- tryCatch(fread(output_file), error=function(e) NULL)
  if(!is.null(existing) && "fid" %in% colnames(existing)){
    if(length(unique(existing$fid)) == nrow(test_fields)){
      cat("[INFO] Final_Test.csv exists and covers all test fields - skipping generation.\n")
      quit(save = "no")
    }
  }
}

#### RAW TEST DATA MATRIX
test <- data.frame()
bands <- c(
  B01 = "B01.B1.tif", B02 = "B02.B2.tif", B03 = "B03.B3.tif",
  B04 = "B04.B4.tif", B05 = "B05.B5.tif", B06 = "B06.B6.tif",
  B07 = "B07.B7.tif", B08 = "B08.B8.tif", B09 = "B09.B9.tif",
  B11 = "B11.B11.tif", B12 = "B12.B12.tif", B8A = "B8A.B8A.tif"
)

for(i in 1:nrow(test_fields)) {
  cat("[INFO] Processing tile:", i, "/", nrow(test_fields), "ID:", test_fields$ids[i], "\n")
  
  mask_file <- test_fields$field_path[i]
  if(!file.exists(mask_file)){
    cat("[WARNING] field_ids.tif not found for", test_fields$ids[i], "\n"); next
  }
  
  f_mat_raster <- raster(mask_file)
  f_mat <- as.matrix(f_mat_raster)
  f_mat_rs <- which(rowSums(f_mat, na.rm = TRUE) > 0)
  f_mat_cs <- which(colSums(f_mat, na.rm = TRUE) > 0)
  if(length(f_mat_rs)==0 || length(f_mat_cs)==0) { cat("[WARNING] Empty mask:", test_fields$ids[i], "\n"); next }

  f_mat <- f_mat[f_mat_rs, f_mat_cs]
  fid_vec <- as.vector(f_mat)
  npix <- length(fid_vec)

  # Extract bands safely
  band_list <- list()
  for(b in names(bands)) {
    img_file <- file.path(data.dir, image_path, paste0(image_path, "_", test_fields$ids[i]), bands[b])
    if(!file.exists(img_file)){
      cat("[WARNING] Band missing:", img_file, "\n")
      band_list[[b]] <- rep(NA_real_, npix)
      next
    }
    rtmp <- raster(img_file)
    mm_tmp <- as.matrix(rtmp)[f_mat_rs, f_mat_cs]
    band_list[[b]] <- as.vector(mm_tmp)
    rm(rtmp, mm_tmp); invisible(gc())
  }
  
  train_data <- as.data.frame(band_list, stringsAsFactors = FALSE)
  dd <- cbind(
    folder = rep(test_fields$ids[i], npix),
    fid = fid_vec,
    train_data
  ) %>% as.data.frame(stringsAsFactors = FALSE)

  colnames(dd) <- make.names(colnames(dd), unique = TRUE)
  cat("[INFO] Tile", i, "columns:", paste(colnames(dd), collapse = ", "), "\n")
  cat("[INFO] Rows for this tile:", nrow(dd), "\n")

  test <- dplyr::bind_rows(test, dd)
  rm(f_mat, fid_vec, dd, train_data, f_mat_rs, f_mat_cs, f_mat_raster); invisible(gc())
}

cat("=== RAW TEST DATA MATRIX CREATED ===\n")
cat("[INFO] Total test rows:", nrow(test), "Columns:", ncol(test), "\n")
if(nrow(test) == 0) stop("ERROR: No test data generated! Cannot continue.")

#### FEATURE ENGINEERING
df_test <- test %>%
  mutate(
    ndvi = (B08 - B04) / (B08 + B04),
    GLI = (2*B03 - B04 - B02) / (2*B03 + B04 + B02),
    CVI = (B08 / B03) * (B05 / B04),
    SIPI = (B08 - B02) / (B08 - B04),
    S2REP = 705 + 35 * (((B07 + B04) / 2 - B05) / (B06 - B05)),
    CCCI = ((B08 / B05) - (B08 / B04)) / ((B08 / B05) + (B08 / B04)),
    hue = atan2(sqrt(3) * (B03 - B02), 2 * B04 - B02 - B03),
    RENDVI = (B08 - B06) / (B08 + B06),
    RECI = (B08 / B05) - 1,
    RECI2 = (B08 / B06) - 1,
    evi = 2.5 * (B08 - B04) / (B08 + 6*B04 - 7.5*B02 + 1),
    evi2 = 2.5 * (B08 - B04) / (B08 + 2.4*B04 + 1),
    npcri = (B04 - B02) / (B04 + B02),
    ndwi = (B03 - B08) / (B03 + B08)
  )
cat("[INFO] Feature engineering completed. Dimensions:", paste(dim(df_test), collapse=" "), "\n")

#### AGGREGATE PER FIELD
df_test <- df_test %>%
  filter(!is.na(fid)) %>%
  group_by(fid) %>%
  mutate(field_tile_count = n(), field_overlap_count = length(unique(folder))) %>%
  ungroup() %>%
  dplyr::select(-any_of(c("folder", "id"))) %>%
  group_by(fid) %>%
  summarise_all(list(median = ~median(.x, na.rm = TRUE),
                     max = ~max(.x, na.rm = TRUE))) %>%
  ungroup()
cat("[INFO] Aggregated df_test dimensions:", dim(df_test), "\n")

#### TILE DETAILS + AREA FEATURES
test_details <- data.frame()
test_area_data <- data.frame()
for(i in 1:nrow(test_fields)){
  folder_name <- paste0(test_path, "_", test_fields$ids[i])
  stac_file <- file.path(data.dir, test_path, folder_name, "stac.json")
  mask_file <- test_fields$field_path[i]
  img_dir <- file.path(data.dir, image_path, paste0(image_path, "_", test_fields$ids[i]))
  
  # --- Tile width/height ---
  bbox <- NULL
  if(file.exists(stac_file)){
    stac <- jsonlite::read_json(stac_file)
    if(!is.null(stac$bbox)) bbox <- unlist(stac$bbox)
  }
  if(is.null(bbox)) bbox <- c(NA,NA,NA,NA)
  tile_width <- bbox[3]-bbox[1]; tile_height <- bbox[4]-bbox[2]
  test_details <- rbind(test_details, data.frame(fid=test_fields$ids[i], tile_width, tile_height, stringsAsFactors=FALSE))
  
  # --- AREA features ---
  band_files <- c(
    B02 = file.path(img_dir, "B02.B2.tif"),
    B03 = file.path(img_dir, "B03.B3.tif"),
    B04 = file.path(img_dir, "B04.B4.tif"),
    B08 = file.path(img_dir, "B08.B8.tif")
  )
  if(all(file.exists(band_files))){
    mm2 <- raster::stack(band_files)
    veg <- c("EVI","SI","GLI","HUE","NDVI","GNDVI")
    area <- c()
    for(j in veg){
      mask_res <- FIELDimageR::fieldMask(mosaic=mm2, Red=3, Green=2, Blue=1, NIR=4, index=j, plot=F, cropAbove=T)
      area_res <- FIELDimageR::fieldArea(mosaic=mask_res$mask, n.core=4, plot=F)
      area <- c(area, area_res$areaPorcent$objArea)
    }
    test_area_data <- rbind(test_area_data, cbind(fid=test_fields$ids[i], t(area)))
  }
}
colnames(test_area_data) <- c("fid", paste0("Area_", c("EVI","SI","GLI","HUE","NDVI","GNDVI")))

#### JOIN TILE DETAILS & AREA FEATURES
df_test <- df_test %>%
  left_join(test_details, by="fid") %>%
  left_join(test_area_data, by="fid")

#### NEIGHBOUR DATA / coordinates + PCA + scaling
neighbour_data <- data.frame()
for(i in 1:nrow(test_fields)){
  mask_file <- test_fields$field_path[i]
  if(!file.exists(mask_file)) next
  f_mat_raster <- raster(mask_file)
  pts <- tryCatch(rasterToPoints(f_mat_raster, spatial=TRUE), error=function(e) NULL)
  if(is.null(pts)) next
  llprj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  llpts <- tryCatch(spTransform(pts, CRS(llprj)), error=function(e) pts)
  f_df <- as.data.frame(llpts)
  idcol <- ifelse("field_ids" %in% colnames(f_df), "field_ids", "layer")
  f_df <- f_df %>% group_by(!!sym(idcol)) %>% summarise_all(list(mean=mean)) %>% ungroup()
  colnames(f_df)[1] <- "fid"
  neighbour_data <- rbind(neighbour_data, f_df)
}

# --- PCA & scale ---
dr.dat <- neighbour_data %>% dplyr::select(x_mean= x_mean.mean, y_mean= y_mean.mean)
trans_pca <- preProcess(dr.dat, method=c("center","scale","pca"))
pc <- predict(trans_pca, dr.dat)
neighbour_data$y_mean_pca <- pc$PC1
neighbour_data$x_mean_pca <- pc$PC2

trans_scale <- preProcess(dr.dat, method=c("center","scale"))
sc <- predict(trans_scale, dr.dat)
neighbour_data$y_mean_scale <- as.numeric(sc$y_mean)
neighbour_data$x_mean_scale <- as.numeric(sc$x_mean)

# --- Compute lat/rot coordinates ---
neighbour_data <- neighbour_data %>%
  mutate(
    lat1 = cos(y_mean_scale)*cos(x_mean_scale),
    lat2 = cos(y_mean_scale)*sin(x_mean_scale),
    lat3 = sin(y_mean_scale),
    rot45_x = 0.707*y_mean_scale + 0.707*x_mean_scale,
    rot45_y = 0.707*y_mean_scale - 0.707*x_mean_scale,
    rot30_x = 0.5*y_mean_scale + 0.866*x_mean_scale,
    rot30_y = 0.5*y_mean_scale - 0.866*x_mean_scale,
    rot60_x = 0.866*y_mean_scale + 0.5*x_mean_scale,
    rot60_y = 0.866*y_mean_scale - 0.5*x_mean_scale,
    x_mean2 = x_mean_scale*60,
    y_mean2 = y_mean_scale*60
  )

df_test <- df_test %>% left_join(neighbour_data %>% select(fid, ends_with("_pca"), ends_with("_scale"), lat1:rot60_y, x_mean2,y_mean2), by="fid")

#### COMPUTE field_tile_size
df_test <- df_test %>% mutate(field_tile_size = field_tile_count * tile_width * tile_height)

#### WRITE FINAL CSV
fwrite(df_test, file=output_file, row.names=FALSE)
cat("[INFO] Final_Test.csv written successfully to:", output_file, "\n")
cat("=== TEST DATA GENERATION COMPLETED ===\n")
