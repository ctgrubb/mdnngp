# This file takes inputted start and end times, and creates a list of files to download from MODIS and NLDAS ftp servers.

args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

listFiles <- function(start, end) {
  p1 <- "https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/NLDAS_FORA0125_H.002/"
  p3 <- "/NLDAS_FORA0125_H.A"
  p5 <- ".002.grb"
  
  start <- strptime(as.character(start), format = "%Y%m%d%H%M", tz = "GMT")
  end <- strptime(as.character(end), format = "%Y%m%d%H%M", tz = "GMT")
  
  times <- seq.POSIXt(start, end, by = "hour")
  
  p2 <- strftime(times, format = "%Y/%j", tz = "GMT")
  p4 <- strftime(times, format = "%Y%m%d.%H%M", tz = "GMT")
  
  NLDASfiles <- paste(p1, p2, p3, p4, p5, sep = "")
  
  years <- strftime(seq.POSIXt(seq(start, length=2, by="-36 days")[2], seq(end, length=2, by="36 days")[2], by = "year"), format = "%Y", tz = "GMT")
  days <- seq(1, 365, by = 8)
  
  times <- unlist(lapply(years, function(x) {
    strftime(seq.POSIXt(from = strptime(paste(x, "-01-01", sep = "", tz = "GMT"), format = "%Y-%m-%d", tz = "GMT"), 
                        to = strptime(paste(x, "-12-31", sep = ""), format = "%Y-%m-%d", tz = "GMT"), 
                        by = "8 days"), format = "%Y%j", tz = "GMT")
    }))
  
  times <- times[which(seq(start, length=2, by="-36 days")[2] <= 
                         strptime(times, format = "%Y%j", tz = "GMT") & 
                         strptime(times, format = "%Y%j", tz = "GMT") <= 
                         seq(end, length=2, by="36 days")[2])]
  
  MODISfiles <- character(2*length(times))
  
  p1 <- "https://e4ftl01.cr.usgs.gov/MOLT/MOD11A2.006/"
  
  for(i in 1:length(times)) {
    p2 <- strftime(strptime(times[i], format = "%Y%j", tz = "GMT"), format = "%Y.%m.%d", tz = "GMT")
    p3 <- system(paste("curl -s https://e4ftl01.cr.usgs.gov/MOLT/MOD11A2.006/", 
                       p2, 
                       "/ | grep -o MOD11A2.A", 
                       times[i], 
                       ".h11v05.006.[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].hdf | head -n 1", 
                       sep = ""), intern = TRUE)
    p4 <- system(paste("curl -s https://e4ftl01.cr.usgs.gov/MOLT/MOD11A2.006/", 
                       p2, 
                       "/ | grep -o MOD11A2.A", 
                       times[i], 
                       ".h12v05.006.[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].hdf | head -n 1", 
                       sep = ""), intern = TRUE)
    MODISfiles[i] <- paste("https://e4ftl01.cr.usgs.gov/MOLT/MOD11A2.006/", p2, "/", p3, sep = "")
    MODISfiles[i+length(times)] <- paste("https://e4ftl01.cr.usgs.gov/MOLT/MOD11A2.006/", p2, "/", p4, sep = "")
  }
  
  cat(NLDASfiles, file = "filelistNLDAS.txt", sep = "\n")
  cat(MODISfiles, file = "filelistMODIS.txt", sep = "\n")
}

listFiles(start, end)