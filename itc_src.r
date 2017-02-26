loadITC <- function(path, proj){
  path.vect <- paste(path, "/vector/", sep="")
  tree.tg <- readOGR(path.vect, "TaggedTrees", stringsAsFactors=F)
  tree.untg <- readOGR(path.vect, "UntaggedTrees", stringsAsFactors=F)
  plot.cnt <- read.centroids("OSBS_diversity_plot_centroids")
  utm.tg <- spTransform(tree.tg, CRS(proj))
  utm.untg <- spTransform(tree.untg, CRS(proj))
  return(list(tg = utm.tg, untg  = tree.untg))
}

read.centroids <- function(file){
  library(rgdal)
  plot.centr <- readOGR("./inputs/", file, stringsAsFactors = FALSE)
  return(plot.centr)
}

get.plot.extent <- function(plots, buffersize){
  centroids <- coordinates(plots)
  xPlus <- centroids[,2]+buffersize
  yPlus <- centroids[,1]+buffersize
  xMinus <- centroids[,2]-buffersize
  yMinus <- centroids[,1]-buffersize  
  ext.pts <- cbind(xMinus, yMinus, xMinus, yPlus, xPlus, yPlus, xPlus, yMinus, xMinus, yMinus)
  ID = plots$Plot_ID
  
  #credits: http://stackoverflow.com/questions/26620373/spatialpolygons-creating-a-set-of-polygons-in-r-from-coordinates
  polys <- SpatialPolygons(mapply(function(poly, id) {
    xy <- matrix(poly, ncol=2, byrow=TRUE)
    Polygons(list(Polygon(xy)), ID=id)
  }, split(ext.pts, row(ext.pts)), ID))
  
  # Create SPDF
  polys.df <- SpatialPolygonsDataFrame(polys, data.frame(id=ID, row.names=ID))
  plot(polys.df, col=rainbow(50, alpha=0.5))
  return(polys.df)
}

loadIMG <- function(path, proj, img.lb = 1){
  path.raster <- paste(path, "/raster/", sep="")
  chm.sample <- raster(paste(path.raster, "chm/", img.lb, "_chm.tif", sep = ""))
  hs.sample <- brick(paste(path.raster, "hs_fullBands/", img.lb, "_nm350_2512.tif", sep = ""))
  #you want to include RGB too, later
  rgb.sample <- raster(paste(path.raster, "camera/", img.lb, "_camera.tif", sep = ""))
  return(list(cmh = chm.sample, hsp = hs.sample, rgb = rgb.sample))
}

filter.bad <- function(img, bb, img.lb){
  img <- dropLayer(img, bb)
  writeRaster(img, paste('./outputs/filtered/', img.lb, "_bb350_2512.tif", sep= ""), format = 'GTiff')
  return(brick(paste('./outputs/filtered/', img.lb, "_bb350_2512.tif", sep= "")))
}

itc.plot <- function(trees, rasters, just.mnf, hybrid = NULL, f=FALSE){
  
  # setwd("/Users/sergiomarconi/Documents/PhD/Projects/scalingHiper/Outputs/")
  #pdf(paste('./outputs/crowns/', f,'_itc.pdf',sep=""),height=12,width=8)
  trees$untg <- spTransform(trees$untg, CRS("+init=epsg:32617"))
  trees$tg <- spTransform(trees$tg, CRS("+init=epsg:32617"))
  just.mnf <- spTransform(just.mnf, CRS("+init=epsg:32617"))
  
  eval.extent <-  union(trees$tg, trees$untg)
  eval.extent <- eval.extent[just.mnf,]
  tst.hsp <- just.mnf[eval.extent,]
  rasters <- crop(rasters, tst.hsp)
  
  if(!f){
    plot(rasters, col = gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
  }else{
    plotRGB(brick(rasters), r=1,g=2,b=3, stretch="Lin")
  }
  plot(tst.hsp,axes=T, border="blue", add=TRUE, lwd = 2)
  plot(eval.extent,axes=T, border="red3", add=TRUE, lwd = 2)
  if(!is.null(hybrid)){
    hybrid <- spTransform(hybrid, CRS("+init=epsg:32617"))
    hybrid <- hybrid[eval.extent,]
    plot(hybrid,axes=T, border="green", add=TRUE, lwd = 1.5)
  }
  #dev.off()
}

itcHhps <- function(hyp.clust, nC, lasITC, par=T, fast=T, verbose = F){
  library(dynamicTreeCut)
  library(cluster)
  for( i in 1:length(lasITC@polygons)){
    p.id <- hyp.clust
    if(verbose){
      plot(lasITC)
      plot(SpatialPolygons(list(lasITC@polygons[[i]])), border = 'red',add=T)
    }
    p.id <- r.mask <- crop(hyp.clust, (SpatialPolygons(list(lasITC@polygons[[i]]))))
    if(verbose){
      plot(p.id[[1]])
    }
    ext <- extent(p.id$pca_plot1.1)
    proj <- p.id@crs
    p.id <- as.array(p.id)
    p.id <- p.id[,,1:nC]
    p.dim <- dim(p.id)
    dim(p.id) <- c(p.dim[1] * p.dim[2], p.dim[3])
    # Ward Hierarchical Clustering
    d <- dist(p.id, method = "manhattan") # distance matrix
    fit <- hclust(d, method="ward.D") 

    d <- daisy(p.id, metric="manhattan")
    groups <- cutreeDynamic(fit, minClusterSize=10, distM = as.matrix(d))
    
    dim(groups) <- c(p.dim[1], p.dim[2])
    dim(p.id) <- c(p.dim[1], p.dim[2], p.dim[3])
    groups[p.id[,,1] ==-1L] <- NA
    if(max(groups, na.rm = TRUE)==0){groups <-groups + 1}
    r <- raster(groups)
    extent(r) <- ext
    crs(r) <- proj
    mask <- rasterize(SpatialPolygons(list(lasITC@polygons[[i]])), r.mask)
    extent(mask) <- ext
    r <- mask(r, mask)
    r2 <- gapfill(r, mask, ext)
    # plot(r2)
    # r2 <- gapfill(r2, mask, ext)
    
    if(length(r2)>0){
      r <- rasterToPolygons(r2, dissolve = T, n = 16)
      if(i ==1){
        poly.list <- r
      }else{
        poly.list <- rbind(poly.list, r)
      }
    }
  }
  return(poly.list)
}

get.k <- function(p.id, par, fast){
  library(pvclust)
  if(fast){
    fit <- pvclust(t(p.id), method.hclust="ward.D",
                   method.dist="manhattan", parallel = par, r = 1)
    pp <-pvpick(fit, pv="bp", alpha =.2)
  }else{
    fit <- pvclust(t(p.id), method.hclust="ward.D",
                   method.dist="manhattan", parallel = par, r = c(0.4,0.8))
    pp <-pvpick(fit, alpha =.95)
  }
  return(length(pp$clusters))
}

gapfill <- function(r,mask, ext){
  # http://stackoverflow.com/questions/24465627/clump-raster-values-depending-on-class-attribute
  # extend raster, otherwise left and right edges are 'touching'
  r <- extend(r, c(1,1))
  mask <- extend(mask, c(1,1))
  # get al unique class values in the raster
  clVal <- unique(r)
  # remove '0' (background)
  clVal <- clVal[!clVal==0]
  # create a 1-value raster, to be filled in with NA's
  r.NA <- setValues(raster(r), 1)
  # set background values to NA
  r.NA <- mask(r, mask, updatevalue = 0)
  r.NA[r==0]<- NA
  # loop over all unique class values
  for (ii in clVal) {
    # create & fill in class raster
    r.class <- setValues(raster(r), NA)
    r.class[r == ii]<- 1
    # clump class raster
    clp <- clump(r.class, directions = 4)
    # calculate frequency of each clump/patch
    cl.freq <- as.data.frame(freq(clp))
    # store clump ID's with frequency 1
    rmID <- cl.freq$value[which(cl.freq$count <= 4)]
    # assign NA to all clumps whose ID's have frequency 1
    r.NA[clp %in% rmID] <- 0
  } 
  
  # multiply original raster by the NA raster
  r2 <- r * r.NA
  
  # Create focal weights matrix.
  nbrhd <- matrix(c(0,1,0,1,1,1,0,1,0), ncol = 3, nrow = 3)
  r3 <- focal(x=r2, w=nbrhd, fun = fill.na)
  r3 <- r3 + r2
  r3 <- crop(r3, ext)
  return(r3)
}

# idea from http://stackoverflow.com/questions/27906111/r-focal-raster-package-how-to-apply-filter-to-subsets-of-background-data
fill.na <- function(x) {
  center <- x[ceiling(length(x)/2)]
  if(!is.na(center)){
    if( center==0 ) {
      fr <-x
      fr <- fr[fr!= 0]
      fr <- fr[!is.na(fr)]
      fr <- as.integer(names(which(table(fr) == max(table(fr)))[1]))
      if(length(fr)>0)
        return((fr))
    }
  }
}

# http://stackoverflow.com/questions/22785475/does-a-function-like-dist-rdist-exist-which-handles-nas
rdist.w.na <- function(X,Y)
{
  library(pdist)
  
  if (!is.matrix(X)) 
    X = as.matrix(X)
  if (!is.matrix(Y)) 
    Y = as.matrix(Y)
  distances <- matrix(pdist(X,Y)@dist, ncol=nrow(X), byrow = TRUE)
  #count NAs
  na.count <- sapply(1:nrow(X),function(i){rowSums(is.na(Y) | is.na(X[i,]))})
  #scaling to number of cols
  distances * sqrt(ncol(X)/(ncol(X) - na.count))
}
