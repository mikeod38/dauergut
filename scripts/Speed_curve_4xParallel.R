# roam/dwell test using wormlab data
# roam/dwell test using wormlab data
packages = c("ggplot2","dplyr","tidyr","reshape2","viridis","plotly", "magrittr", "pryr", "devtools", "data.table")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

devtools::install_github("hadley/lineprof")
library(lineprof)

##### setup  make sure to change working dir to location of files - only 1 of each per folder#####

#### to run parallel #####
library(parallel)

# get list of folders - set wd only folders with wormlab data should be present
folder.list <- list.files()
folder.list
directory <- getwd()
# Calculate the number of cores
no_cores <- min(detectCores() - 2,4) # uses a ton of RAM right now

# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist = c("folder.list", "directory"))
library(lineprof)

# change num.tracks
#prof <- lineprof(
system.time(
  #lapply(folder.list, function(folder,num.tracks) {
parLapply(cl,folder.list, function(folder,num.tracks) {
  position.path <- file.path(folder,list.files(path = file.path(folder),pattern = c("osition",".csv")))
  direction.path <- file.path(folder,list.files(path = file.path(folder),pattern = c("irection",".csv")))
  speed.path <- file.path(folder,list.files(path = file.path(folder),pattern = c("peed",".csv")))
  bin.length <- 10
  frame.rate <- 3
  time <- 5400
  xend <- 100
  yend <- 200
  slope <- 200/100
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(magrittr)
  ##### include functions: #####
  ######fxn to melt WL position data to long format:
  WL.pos.long <- function(data, num.tracks) {
    subset.long <- data[,1:(num.tracks*2 + 2)] %>% melt(id.vars = c(1,2)) %>%
      separate(variable, sep = "\\.", c("worm", "pos")) %>% dcast(Frame + Time + worm ~ pos)
    return(subset.long)
  }
  ######fxn to calculate 3 point curvature based on law of cosines:
  curve.angle <- function(del.x1, del.y1, del.x2, del.y2) {
    values <- list(del.x1, del.y1, del.x2, del.y2)
    if(anyNA(values)) {
      "NA"
    } else {
      x <- c(del.x1, del.y1)
      y <- c(del.x2, del.y2)
      dot.prod <- x%*%y 
      norm.x <- as.numeric(svd(x)[1]) # faster to use svd and index (which for 1x1 vec is all norm does)
      #norm(x,type="2") # length of vector
      norm.y <- as.numeric(svd(y)[1])
      #norm(y,type="2")
      theta <- acos(dot.prod / (norm.x * norm.y))
      as.numeric(theta)
    }
  }
  
  ######main fxn to analyze and plot the data
  WL.roam.data.vectorized <- function(position.path, speed.path, direction.path, bin.length, frame.rate, num.tracks) {

    ####Setting up data#########
    speed.data <-read.csv(speed.path, skip = 4)
    if(missing(num.tracks)) {
      num.tracks <- length(speed.data) - 2
    } else {
      num.tracks = num.tracks
    }
    vid.length <- max(speed.data$Frame)
    bin.length <- bin.length # bin length in s
    frame.rate <- frame.rate # usually 2 or 3 fps
    bin.size <- bin.length*frame.rate
    n.bins <- vid.length/bin.size
    ############################
   
    #### get speed data ####
    WL.speed <- speed.data[,1:(num.tracks + 2)] %>% 
      melt(id.vars = c(1,2)) %>% dplyr::filter(!is.na(value)) %>%
      separate(variable, sep = "\\.", c("worm", "stuffer")) %>% dplyr::select(-stuffer) %>%
      rename(speed = value)
    rm(speed.data)
    
    #### get position data ##########
    position <- read.csv(position.path, skip = 4)
    WL.centroid <- WL.pos.long(position, num.tracks) %>% 
      dplyr::filter(!is.na(x)) %>% mutate(type = "centroid")
    rm(position)
    
    #### get direction data #### 
    direction <- read.csv(direction.path, skip = 4)
    WL.head.dir <- direction[,1:(num.tracks + 2)] %>%
      melt(id.vars = c(1,2)) %>% dplyr::filter(!is.na(value)) %>%
      separate(variable, sep = "\\.", c("worm", "stuffer")) %>% dplyr::select(-stuffer) %>%
      rename(head.dir = value)
    rm(direction)
    
    
    #### merge data ####
    WL.alldata <- dplyr::bind_cols(list(WL.centroid, WL.speed, WL.head.dir)) %>% 
      subset(., select=which(!duplicated(names(.))))
    rm(WL.pos.long, WL.speed, WL.head.dir)
    # used bind_cols instead of merge or join to speed up
    # thus, cannot alter row ordering before this step. 
    
    # WL.alldata <- list(WL.centroid,WL.speed,WL.head.dir) %>% combine() %>% arrange(worm, Time) %>% mutate(stuffer, = NULL) ### this takes longest
    #    Reduce(function(...) dplyr::full_join(...), .) #%>% arrange(worm, Time) %>% mutate(stuffer = NULL) ### this takes longest
    ##################################
    WL.alldata %<>% group_by(worm) %>% 
      mutate(track.frame = row_number()) %>% # make 'track.frame' variable, = count every 30 frames by track
      mutate(time.bin = ceiling(track.frame/bin.size)) %>% # round up to integer ie 0.1 = 1, 1.1 = 2
      mutate(del.y2 =  y - lag(y), # change from previous point (t-1) to (t0)
             del.x2 = x - lag(x),
             del.x1 = lag(x) - lag(x, n=2), #vector from t(-2) to t(-1) for curve angle
             del.y1 = lag(y) - lag(y, n=2), 
             curve.ang = as.numeric(mapply(curve.angle, del.x1, del.y1, del.x2, del.y2))*180/pi) %>% group_by(worm, time.bin) %>%
      mutate(bin.speed = abs(mean(speed, na.rm=TRUE)), bin.ang.vel = mean(curve.ang, na.rm=TRUE)) %>% filter(!is.na(curve.ang)) 
    ###################################
    
    return(WL.alldata)
  }
  plot_density <- function(data, xend, yend) {
    truncated <- data %>% group_by(worm) %>% summarise(n = n()) %>% dplyr::filter(n<10)
    data[!data$worm %in% truncated$worm,] %>% group_by(worm,time.bin) %>%
      dplyr::filter(bin.speed < 500) %>%
      summarize(mean.speed = mean(bin.speed),mean.angle = mean(bin.ang.vel)) %>% 
      ggplot(aes(x = mean.angle, y = mean.speed)) + 
      stat_density2d(geom="raster", aes(fill = ..density..), contour = FALSE)  + 
      viridis::scale_fill_viridis(option = "inferno", begin = 0.05, end = 0.9) + 
      labs(title = "Density plot") +
      #geom_point(alpha = 0.1) + 
      coord_cartesian(xlim = c(0,150),ylim = c(0,250)) +
      geom_segment(aes(x=0, y=0, xend = xend, yend = yend), colour = "red") + theme_classic()
  }
  plot_tracks <- function(data, time) {
    data %>% dplyr::filter(Time < time, bin.ang.vel < 150) %>%
      ggplot(aes(x = x, y = y)) + 
      geom_point(aes(colour = Time), alpha = 0.01) +
      labs(title = "All tracks") +
      viridis::scale_color_viridis(option = "inferno") +theme_classic()#+ facet_wrap(~worm) #to plot each track
  }
  plot_scatter <- function(data, xend, yend) {
    truncated <- data %>% group_by(worm) %>% summarise(n = n()) %>% dplyr::filter(n<10)
    data[!data$worm %in% truncated$worm,] %>% group_by(worm,time.bin) %>%
      dplyr::filter(bin.speed < 500) %>%
      summarize(mean.speed = mean(bin.speed),mean.angle = mean(bin.ang.vel)) %>% 
      ggplot(aes(x = mean.angle, y = mean.speed)) + 
      labs(title = "Scatter plot by time bin") +
      geom_point(alpha = 0.1) + 
      coord_cartesian(xlim = c(0,150),ylim = c(0,250)) +
      geom_segment(aes(x=0, y=0, xend = xend, yend = yend), colour = "red") + theme_classic()
  }
  plot_position_changes <- function(data) {
    truncated <- data %>% group_by(worm) %>% summarise(n = n()) %>% dplyr::filter(n<10)
    p1<-data[!data$worm %in% truncated$worm,] %>% group_by(worm) %>%
      summarise(min = min(abs(del.y1)), mean = mean(abs(del.y1))) %>% 
      ggplot() + 
      geom_histogram(aes(x=mean), bins=1000) +
      labs(x = "change in y (microns)") + 
      theme_classic()# + coord_cartesian(xlim = c(0,1))
    
    p2<-data[!data$worm %in% truncated$worm,] %>% group_by(worm) %>%
      summarise(min = min(abs(del.x1)), mean = mean(abs(del.x1))) %>% 
      ggplot() + 
      geom_histogram(aes(x=mean), bins=1000) +
      labs(x = "change in x (microns)") + 
      theme_classic()# + coord_cartesian(xlim = c(0,1))
    
    gridExtra::grid.arrange(p2,p1, ncol = 2, 
                            top = "Histograms of position changes by frame")
  }
  
  #### wrapper for all plot functions above ####
  plot_all_vect <- function(folder,data, time, xend, yend) {
    p1 <- plot_density(data, xend, yend)
    p2 <- plot_scatter(data, xend, yend)
    p3 <- plot_position_changes(data)
    p4 <- plot_tracks(data, time)
    p<-gridExtra::grid.arrange(p1,p2,p3,p4, ncol = 2, nrow =2)
    ggsave(p, filename =file.path(folder,"plots.pdf"), width = 11, height = 8.5, units = "in",device = "pdf")
  }
  
  #### fxn to get pct roam ####
  roam.pct <- function(data, slope) {
    truncated <- data %>% group_by(worm) %>% summarise(n = n()) %>% dplyr::filter(n<10)
    roam <- data[!data$worm %in% truncated$worm,] %>% group_by(worm,time.bin) %>% 
      summarize(mean.speed = mean(bin.speed),mean.angle = mean(bin.ang.vel)) %>% 
      mutate(ratio = mean.speed/mean.angle) %>% dplyr::filter(ratio >= slope) %>% nrow()
    dwell <- data[!data$worm %in% truncated$worm,] %>% group_by(worm,time.bin) %>% 
      summarize(mean.speed = mean(bin.speed),mean.angle = mean(bin.ang.vel)) %>% 
      mutate(ratio = mean.speed/mean.angle) %>% dplyr::filter(ratio < slope) %>% nrow()
    pct.roam <- roam/(dwell+roam)
    print(pct.roam)
    return(data.frame(roam = roam, dwell = dwell, pct.roam = pct.roam))
  }

    
  
  #### run scripts ####   
    if(missing(num.tracks)) {
      data <- WL.roam.data.vectorized(position.path = position.path,
                                      direction.path = direction.path,
                                      speed.path = speed.path,
                                      bin.length, frame.rate)
    } else {
      num.tracks = num.tracks
      data <- WL.roam.data.vectorized(position.path = position.path,
                                      direction.path = direction.path,
                                      speed.path = speed.path,
                                      bin.length, frame.rate, num.tracks)
    }
  
  #### save data ####
    data.table::fwrite(data, file.path(folder,"all_track_data.csv"))
    plot_all_vect(folder,data,time,xend,yend)
    summary <- roam.pct(data, slope)
    data.table::fwrite(summary, file.path(folder,"roam_dwell.csv"))
  }
)
)

stopCluster(cl)
