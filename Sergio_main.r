#clear workspace
rm(list=ls(all=TRUE))
debug = F

if (debug){
  minClusterSize=15
  pt.root = "/Users/sergiomarconi/Documents/PhD/Projects/scalingHiper/data/NIST_data_20170120"
  proj = "+init=epsg:32617 +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  conWin <- 5
  NumComps <- 20
  plot.side <- 80
}

main <- function(pt.root = NULL, proj = NULL, minClusterSize=15, conWin = 5, NumComps = 10, plot.side = 80, buffersize = 40){
  setwd("/Users/sergiomarconi/Documents/PhD/Projects/scaleUpOSBS/")
  #load external libraries
  library(neonAOP)
  library(itcSegment)
  library(raster)
  library(readr)
  library(rgdal)
  library(RStoolbox)
  
  #load project's functions
  source(paste(getwd(), '/src/itc_src.r', sep=""))
  source(paste(getwd(), '/src/mnf.r', sep=""))
  
  #loading vector dataset
  trees <- loadITC(pt.root, proj)
  polys.df <- get.plot.extent(plots = read.centroids("OSBS_diversity_plot_centroids"), buffersize)
  
  #defining bad bands
  bad_bands <- read.table("~/Documents/PhD/Projects/scaleUpOSBS/inputs/bad_bands.csv", header = FALSE)
  
  if(debug){j <- as.character(polys.df$id[27])}
  
  
  for (j in as.character(polys.df$id[27])) {
    rasters <- loadIMG(pt.root, proj, img.lb = j)
    #remove bad bands
    if(!file.exists(paste('./outputs/filtered/', j, "_bb350_2512.tif", sep= ""))){
      bb <- paste('OSBS_001_nm350_2512', t(bad_bands), sep=".")
      rasters$hsp <- filter.bad(rasters$hsp, bb,  img.lb = j)
    }else{
      rasters$hsp <- brick(paste('./outputs/filtered/', j, "_bb350_2512.tif", sep= ""))
    }
    if(!file.exists('./outputs/pca_plot1.tif')){
        pca <- rasterPCA(rasters$hsp)
    }
    pca <- brick('./outputs/pca_plot1.tif')
    
    # las given itc. here comes the tuning part
    xyz <- read_csv("~/Documents/PhD/Projects/scaleUpOSBS/inputs/xyz.txt", 
                    col_names = FALSE)
    colnames(xyz) <- c("X", "Y", "Z")
    lasITC <- itcLiDAR(X = xyz$X, Y = xyz$Y, Z = xyz$Z, epsg = 32617, resolution = 0.9, 
                       MinSearchFilSize = 3, MaxSearchFilSize = 7, TRESHSeed = .8, 
                       TRESHCrown = 0.7, minDIST = 5, maxDIST = 60, HeightThreshold = 2)
    poly.list <- itcHhps(pca,NumComps, lasITC)
    itc.plot(trees, rasters$rgb, lasITC, poly.list, F)
  }
}
main(pt.root = "/Users/sergiomarconi/Documents/PhD/Projects/scalingHiper/data/NIST_data_20170120", 
     proj = "+init=epsg:32617 +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0", minClusterSize=15,
     conWin = 5, NumComps = 10, plot.side = 80, buffersize = 40)


