# roam/dwell test using wormlab data
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(viridis)
library(plotly)



WL.roam.data <- function(bin.length, frame.rate, num.tracks) {
  #### import data ###
  message("Choose position data")
  position <- read.csv(file.choose(), skip=4)
  message("Choose Direction data")
  direction <- read.csv(file.choose(), skip=4)
  message("Choose speed data")
  speed <- read.csv(file.choose(), skip=4)
  
  ####Setting up parameters####
  if(missing(num.tracks)) {
    num.tracks <- length(speed) - 2
  } else {
    num.tracks = num.tracks
  }
  vid.length <- max(position$Frame)
  bin.length <- bin.length # bin length in s
  frame.rate <- frame.rate # usually 2 or 3 fps
  bin.size <- bin.length*frame.rate
  n.bins <- vid.length/bin.size
  ##############################
  
  
  ######fxn to melt WL data to long format#########
  WL.pos.long <- function(data, num.tracks) {
    subset.long <- data[,1:(num.tracks*2 + 2)] %>% melt(id.vars = c(1,2)) %>%
      separate(variable, sep = "\\.", c("worm", "pos")) %>% dcast(Frame + Time + worm ~ pos)
    return(subset.long)
  } 
  
  ######fxn to calculate 3 point curvature based on law of cosines ####
  curve.angle <- function(del.x1, del.y1, del.x2, del.y2) {
    values <- list(del.x1, del.y1, del.x2, del.y2)
    if(anyNA(values)) {
      "NA"
    } else {
      x <- c(del.x1, del.y1)
      y <- c(del.x2, del.y2)
      dot.prod <- x%*%y 
      norm.x <- norm(x,type="2")
      norm.y <- norm(y,type="2")
      theta <- acos(dot.prod / (norm.x * norm.y))
      as.numeric(theta)
    }
  }
  #################################################
  
  
  ##### merge and get data ##########
  WL.centroid <- WL.pos.long(position, num.tracks) %>% mutate(type = "centroid")
  
  WL.speed <- speed[,1:(num.tracks + 2)] %>% melt(id.vars = c(1,2)) %>%  separate(variable, sep = "\\.", c("worm", "stuffer")) %>% rename(speed = value)
  
  WL.head.dir <- direction[,1:(num.tracks + 2)] %>% melt(id.vars = c(1,2)) %>%  separate(variable, sep = "\\.", c("worm", "stuffer")) %>% rename(head.dir = value)
  
  WL.alldata <- list(WL.centroid,WL.speed,WL.head.dir) %>%
    Reduce(function(...) merge(..., all = T), .) %>% arrange(worm, Time) %>% mutate(stuffer = NULL)
  
  WL.alldata <- WL.alldata %>% group_by(worm) %>%
    mutate(del.y =  y - lag(y), # change from previous point
           del.x = x - lag(x),
           del.x1 = lag(x) - lag(x, n=2), #vector from t(-2) to t(-1) for curve angle
           del.y1 = lag(y) - lag(y, n=2), 
           del.x2 = x - lag(x,n=2), #vector from t(-2) to t(0) for curve angle 
           del.y2 = y - lag(y,n=2),
           time.bin = ntile(Time, n.bins), 
           curve.ang = as.numeric(mapply(curve.angle, del.x1, del.y1, del.x2, del.y2))*180/pi) %>% group_by(worm, time.bin) %>%
    mutate(bin.speed = mean(abs(speed), na.rm=TRUE), bin.ang.vel = mean(curve.ang, na.rm=TRUE)) %>% filter(!is.na(curve.ang))
  ###################################
  return(WL.alldata)
}

# get roaming data
system.time(WL.alldata<-WL.roam.data(bin.length = 10, frame.rate = 3))

# to plot all points - takes a while
# WL.alldata %>% ggplot(aes(x = bin.ang.vel, y = abs(bin.speed))) +
#   geom_point(aes(colour = Time), alpha = 0.01) +
#   scale_color_viridis(option = "inferno")#to plot each track

#plot denisty map of points: 
WL.alldata %>% ggplot(aes(x = bin.ang.vel, y = abs(bin.speed))) + stat_density2d(geom="raster", aes(fill = ..density..), bins=100,contour = FALSE)  + 
  scale_fill_viridis(option = "inferno") + coord_cartesian(xlim = c(0,90), ylim = c(0,250)) + geom_segment(aes(x=0, y=0, xend = 75, yend = 126), colour = "red")

# 
#+ facet_wrap(~ntile(Time,4)) #to plot each track



#plot all selected tracks 
ggplotly(WL.alldata %>% dplyr::filter(Time < 1200) %>% ggplot(aes(x = x, y = y)) + geom_point(aes(colour = ..density..), alpha = 0.2) +
  scale_color_viridis(option = "inferno") + facet_wrap(~ntile(Time,4))) #to plot each track

p<-WL.alldata %>% ggplot(aes(x=ntile(y, 800)))
grid.arrange(p+stat_summary(fun.y = "mean", geom = "point", aes(y=abs(speed))),
             p+stat_summary(fun.y = "mean", geom = "point", aes(y=abs(curve.ang))),
             p+stat_summary(fun.y = "mean", geom = "point", aes(y=bin.speed)),
             p+stat_summary(fun.y = "mean", geom = "point", aes(y=bin.ang.vel)))

ggplotly(WL.alldata %>% dplyr::filter(Time < 30, ang.vel < 100) %>% ggplot(aes(x = x, y = y)) + geom_point(aes(colour = Time), alpha = 0.5) +
               scale_color_viridis(option = "inferno")) #+ facet_wrap(~worm) #to plot each track


apply(points[3,c('del.x1','del.y1', 'del.x2', 'del.y2')], 1, function(y) curve.angle(points['del.x1'],points['del.y1'],points['del.x2'], points['del.y2']))

points1 <- data.frame(x = c(2,3,3), y = c(3,4,5), worm = "X1") %>% mutate(del.x1 = lag(x) - lag(x, n=2), del.y1 = lag(y) - lag(y, n=2), del.x2 = x - lag(x, n=2), del.y2 = y - lag(y,2))
points2 <- data.frame(x = c(3,7,2), y = c(3,3,3), worm = "X2") %>% mutate(del.x1 = lag(x) - lag(x, n=2), del.y1 = lag(y) - lag(y, n=2), del.x2 = x - lag(x, n=2), del.y2 =NA)


points <- merge(points1, points2, all = TRUE)

#mapply(curve.angle, points$del.x1, points$del.y1, points$del.x2, points$del.y2)

points %>% group_by(worm) %>% mutate(theta = mapply(curve.angle, del.x1, del.y1, del.x2, del.y2))

points %>% group_by(worm) %>% mutate(theta = do.call( function(del.x1, del.y1, del.x2, del.y2) curve.angle(...)))

theta <- do.call( function(del.x1, del.y1, del.x2, del.y2) curve.angle(del.x1, del.y1, del.x2, del.y2), points)

