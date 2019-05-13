library(data.table)

if(Sys.info()["nodename"] == "broce") {                                          
  sourcefolder <- "~/git-repos/mdnngp/R/"
  datafolder <- "/storage/chris/mdnngp/"
} else {                                                                         
  stop("Unknown host detected. Please declare folder locations.")                
}

misc_folder <- paste0(datafolder, "Misc/")
modis_folder <- paste0(datafolder, "MODISTEXT/")
modis_times_folder <- paste0(datafolder, "MODISTIME/")
nldas_folder <-paste0(datafolder, "NLDASTEXT/")

modis_files_h11 <- list.files(modis_folder, pattern = "h11")
modis_times_files_h11 <- list.files(modis_times_folder, pattern = "h11")
modis_files_h12 <- list.files(modis_folder, pattern = "h12")
modis_times_files_h12 <- list.files(modis_times_folder, pattern = "h12")
nldas_files <- list.files(nldas_folder)


modis_lls <- fread(paste0(misc_folder, "MODIS_revised.txt"), header = TRUE)
nldas_lls <- fread(paste0(misc_folder, "NLDAS_revised.txt"), header = TRUE)

bbox_lat <- c(37.1, 37.6)
bbox_long <- c(-80.3, -79.8)

nldas_lls$BBOX <- (nldas_lls$Latitude > bbox_lat[1]) &
  (nldas_lls$Latitude < bbox_lat[2]) &
  (nldas_lls$Longitude > bbox_long[1]) &
  (nldas_lls$Longitude < bbox_long[2])

modis_lls$BBOX <- (modis_lls$Latitude > bbox_lat[1]) &
  (modis_lls$Latitude < bbox_lat[2]) &
  (modis_lls$Longitude > bbox_long[1]) &
  (modis_lls$Longitude < bbox_long[2])

nldas_data <- nldas_lls[nldas_lls$BBOX, -4]
modis_data <- modis_lls[modis_lls$BBOX, -c(4, 5)]
modis_time <- modis_lls[modis_lls$BBOX, -c(4, 5)]

for(i in nldas_files) {
  nldas <- fread(paste0(nldas_folder, i), header = TRUE, na.strings = "9.999e+20")$`464 224`
  nldas_data <- cbind(nldas_data, nldas[nldas_lls$BBOX])
}

for(j in 1:length(modis_files_h11)) {
  modis_h11 <- fread(paste0(modis_folder, modis_files_h11[j]))$V3
  modis_h12 <- fread(paste0(modis_folder, modis_files_h12[j]))$V3
  modis <- c(modis_h11, modis_h12)
  modis_data <- cbind(modis_data, modis[modis_lls$BBOX])
  
  modis_times_h11 <- fread(paste0(modis_times_folder, modis_times_files_h11[j]), na.strings = "255")$V3
  modis_times_h12 <- fread(paste0(modis_times_folder, modis_times_files_h12[j]), na.strings = "255")$V3
  modis_times <- c(modis_times_h11, modis_times_h12) * 0.1
  modis_times <- modis_times[modis_lls$BBOX]
  modis_hours <- floor(modis_times)
  modis_h0s <- rep("", length(modis_hours))
  modis_h0s[nchar(modis_hours) == 1] <- "0"
  modis_minutes <- round((modis_times %% 1) * 60)
  modis_m0s <- rep("", length(modis_minutes))
  modis_m0s[nchar(modis_minutes) == 1] <- "0"
  
  modis_dtimes <- paste0(modis_h0s, modis_hours, ":", modis_m0s, modis_minutes)
  modis_time <- cbind(modis_time, modis_dtimes)
}

nldas_data <- as.data.frame(nldas_data)
modis_data <- as.data.frame(modis_data)
modis_time <- as.data.frame(modis_time)

nldas_data[, -c(1, 2, 3)] <- nldas_data[, -c(1, 2, 3)] - 273.15
modis_data[, -c(1, 2, 3)] <- modis_data[, -c(1, 2, 3)] * 0.02 - 273.15
modis_data[modis_data == -273.15] <- NA
modis_time[modis_time == "NA:NA"] <- NA

nldas_dates <- as.POSIXct(substr(nldas_files, start = 19, stop = 31), format = "%Y%m%d.%H%M", tz = "GMT")
modis_dates <- as.Date(substr(modis_files_h11, start = 10, stop = 16), format = "%Y%j")

colnames(nldas_data)[-c(1, 2, 3)] <- as.character(nldas_dates)
colnames(modis_data)[-c(1, 2, 3)] <- as.character(modis_dates)
colnames(modis_time)[-c(1, 2, 3)] <- as.character(modis_dates)

write.table(nldas_data, paste0(datafolder, "nldas_data.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
write.table(modis_data, paste0(datafolder, "modis_data.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
write.table(modis_time, paste0(datafolder, "modis_time.csv"), sep = ",", row.names = FALSE, col.names = TRUE)
