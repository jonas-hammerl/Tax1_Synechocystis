# This script will integrate radial intensity data obtained with the custom
# ImageJ macro using the Radial Profile Angle Plugin and the single cell
# tracking data obtained using the TrackMate Plugin. In the end it will identify
# the directions of displacement as well as the direction of highest
# fluorescence intensity per cell. 

# Load the required packages. If not installed you can do so with 
# install.packages("name of missing package")
library(tidyverse)
library(circular)
library(openair)
library(viridis)
library(readxl)
library(jpeg)
library(gridExtra)
library(RcppRoll)
library(RANN)
library(vroom)
library(activity)


### RADIAL INTENSITY INFORMATION

# Set working directory of the .txt files from the Radial Profile Angle macro
# in ImageJ. Example directory.
# typically the  hierarchy I use is:
# > ...
#   > analysis
#     > specific-video_name
#     'file: float-correction'
#       > radial-intensity
#       > jpeg
#       'file: tracking'
#       'file: info'
setwd("//myr/home/jhammerl/Documents/microscopy/pixG-localization/20231026_pixG/analysis/20231026_pixG-eyfp_k3_1_R255_1_eyfp/radial-intensity")

# Read all result .txt files and bind them into a single table "intens". The 
# first column is only the numeration and thus omitted.
results <- list.files(pattern = ".txt")
intens <- vroom(results)
intens <- intens[,2:length(intens)]


# Set the working directory to the path where the "info.csv" from ImageJ, the 
# motility lab sheet with the tracking information (from TrackMate) and the 
# folder with the .jpg files (jpeg of the original timelapse stack to visualize 
# the calculated angles later) are saved. 
setwd("//myr/home/jhammerl/Documents/microscopy/pixG-localization/20231026_pixG/analysis/20231026_pixG-eyfp_k3_1_R255_1_eyfp/")

# Read the info file and rename the numeration column to "cell". Note that 
# cell is an enumeration variable here and does not link/track the same cell 
# over time. Cell 1 is simply the first particle detected by the Analyze
# Particles function of ImageJ. In the current version of the script the data
# from this info file is not used anymore, but it is appended here and could be 
# used.
info <- read.csv("info.csv")
info <- info %>% rename(cell=X.1)


# Convert the angle column to a circular number in degrees and group the data.
# The data is then filtered to obtain only the measurements that are between 
# 30% and 75% of the radius for each cell as we are interested in the intensity
# at the membrane here which should typically fall into this range (in the
# measurement the major radius from a dilated binary mask was used which would 
# correspond to 100% radius)
intens2 <- intens %>% 
  mutate(angle = circular(angle, units = "degrees")) %>%
  group_by(cell, angle, timepoint) %>%
  filter(radius > quantile(radius, .3), radius < quantile(radius,.75)) 

# The number of remaining NA observations. Typically, the number should be 0 or 
# very low as the NAs occur at the very small radii filtered out before.
sum(is.na(intens2$int.intensity.norm))

# Remove cells that still contain NAs in the "int.intensity.norm" column and 
# join the new data with the info table derived from the analyze particles 
# function in ImageJ.
intens2 <- na.omit(intens2)
intens2 <- left_join(intens2, info, by = "cell")


# Define a function to calculate the intersection angles of two circles (in this
# case two cells). This is later used to exclude the sectors where the 
# circles around two neighboring cells, in which the intensity is measured 
# originally, intersect. The reason is that signal intensity is usually very 
# high at the septum of dividing cells or at the touching interface of cells 
# that are close together. It may be that this is an artifact or really due to 
# protein localization, but it is not informative for localization with regard 
# to the light direction and thus adds a lot of noise to the data especially if 
# dividing cells are detected as two different cells.
circle_intersection <- function(x0, y0, r0, x1, y1, r1){
  d <- sqrt((x1-x0)^2 + (y1-y0)^2)
  if(d > r0 + r1){
    return(0)
  } else {
    d <- sqrt((x1-x0)^2 + (y1-y0)^2)
    a <- (r0^2-r1^2+d^2)/(2*d)
    h <- sqrt(r0^2-a^2)
    x2 <- x0+a*(x1-x0)/d
    y2 <- y0+a*(y1-y0)/d
    x3.1 <- x2+h*(y1-y0)/d
    x3.2 <- x2-h*(y1-y0)/d
    y3.1 <- y2-h*(x1-x0)/d
    y3.2 <- y2+h*(x1-x0)/d
    c(x3.1,y3.1,x3.2,y3.2)
    a1 <- abs(circular(coord2rad(x3.1-x0, y3.1-y0)*360/(2*pi), 
                       units="degrees")-360)
    a2 <- abs(circular(coord2rad(x3.2-x0, y3.2-y0)*360/(2*pi), 
                       units="degrees")-360)
    a3 <- abs(circular(coord2rad(x3.1-x1, y3.1-y1)*360/(2*pi), 
                       units="degrees")-360)
    a4 <- abs(circular(coord2rad(x3.2-x1, y3.2-y1)*360/(2*pi), 
                       units="degrees")-360)
    return(c(a1,a2,a3,a4))

  }
}

# Summarize the intensity data by adding the intensity to a total intensity that
# is now remaining in the given sectors and radius spans. X, Y and the maximum
# radius are retained. The radius is converted from pixels (given by the Radial 
# Profile Angle plugin in ImageJ) to µm. The value of 0.0965332 µm/pixel was 
# retrieved manually from ImageJ. It may need to be adjusted if the images are 
# acquired with other magnification and resolution settings.
intens2 <- intens2 %>% 
  summarize(total.int = sum(int.intensity.norm), 
            X = unique(X.x), 
            Y = unique(Y.x),
            r = max(radius) * 0.0965332)

# Initialize two new columns that will contain the intersection angles with the
# closest neighboring cell.
intens2$intersect.ang1 <- NA  
intens2$intersect.ang2 <- NA

# Summarize intens2 by timepoint. This is the data which is queried to find the 
# closest non-self cell and contains the coordinates of each cell at each 
# timepoint.
pos1 <- intens2 %>% 
  ungroup() %>% 
  group_by(cell) %>% 
  summarize(X = unique(X), 
            Y = unique(Y), 
            timepoint = unique(timepoint))

# Create an empty object to afterwards contain the results of the nearest 
# neighbor search.
nn_res <- NULL

# Search for the nearest neighbor of each cell using nn2() the from RANN 
# package. The search is looped over all timepoints. 
for (i in unique(pos1$timepoint)){
  
  # nn2() requires the data to be ungrouped. Create a subset for each timepoint. 
  pos2 <- ungroup(pos1[pos1$timepoint == i,])
  
  # nn2() requires two data sets of X and Y coordinates, here both pos2. k = 2 
  # will allow to find the two closest neighbors where the best hit is always 
  # the cell itself and the second one the closest non-self cell. Only the 
  # non-self cell index (.$nn.idx[,2]) is selected and stored in "z"
  z <- nn2(dplyr::select(pos2,c("X","Y")), 
           dplyr::select(pos2,c("X","Y")), 
           k=2, 
           searchtype = "radius",
           radius=5)$nn.idx[,2]
  
  # Replace no hit indices (0) by NA
  z[z == 0] <- NA
  
  # Create an object that contains the cell, its nearest neighbor cell and the 
  # timepoint queried.
  nn_res_t <- bind_cols(nn=pos2$cell[z],pos2[,c("cell","timepoint")])
  
  # Bind the nearest neighbor information of all timepoints together into one 
  # object. 
  nn_res <- bind_rows(nn_res,nn_res_t)
}

# Add the nearest neighbor information to the original intens2 table.
intens2 <- left_join(intens2,nn_res, by=c("cell","timepoint"))

# The following loop calculates the intersection angles of all cells and their 
# nearest neighbors if they exist. The loop is over all cells that have a 
# neighbor in the nn2().
for(i in unique(intens2 %>% 
                filter(!is.na(nn)) %>% 
                pull(cell))){
  
  # Retrieve the coordinates and radii of the cell and its nearest neighbor.
  x0 <- unique(intens2$X[intens2$cell==i])
  y0 <- unique(intens2$Y[intens2$cell==i])
  r0 <- unique(intens2$r[intens2$cell==i])
  x1 <- unique(intens2$X[intens2$cell==unique(intens2[intens2$cell==i,]$nn)])
  y1 <- unique(intens2$Y[intens2$cell==unique(intens2[intens2$cell==i,]$nn)])
  r1 <- unique(intens2$r[intens2$cell==unique(intens2[intens2$cell==i,]$nn)])
  
  # Use the circle_intersection function defined before to calculate angles, at
  # which the circles around a cell and its nearest neighbor intersect.
  int <- circle_intersection(x0,y0,r0,x1,y1,r1)
  
  # Only if the two circles intersect an angle is pasted in the before 
  # initialized columns otherwise the column value remains at NA.
  if(sum(int) > 0){
  intens2$intersect.ang1[intens2$cell==i] <- int[1]
  intens2$intersect.ang2[intens2$cell==i] <- int[2]
  }
}

# The intersection angles are rounded to the next full 10.
intens2$intersect.ang1 <- circular(ceiling(intens2$intersect.ang1/10)*10, 
                                   units = "degrees")
intens2$intersect.ang2 <- circular(floor(intens2$intersect.ang2/10)*10, 
                                   units = "degrees")

# The intensity at angles where two neighboring cells intersect is excluded from
# the further analysis. This is done here by setting the intensity for all 
# angles between the calculated intersection angles to 0.
intens2 <- intens2 %>% 
  mutate(total.int=case_when(intersect.ang1-intersect.ang2 > 0 & 
                               angle >= intersect.ang2 & 
                               angle <= intersect.ang1 ~ 0,
                             intersect.ang1-intersect.ang2 < 0 & 
                               angle >= intersect.ang2 ~ 0,
                             intersect.ang1-intersect.ang2 < 0 & 
                               angle <= intersect.ang1 ~ 0,
                             T ~ total.int))

# Calculate the intensity sum of 9 adjacent angles per cell as a new column 
# "rollsum" and the mean angle of the 9 adjacent angle sectors as 
# mean.roll.angle. To go over 0°, "total.int" is repeated twice to a length of 
# 44 (covering every starting angle only once) to allow for the combination of 
# angles starting with 350° as the maximum starting angle.
intens2 <- intens2 %>% 
  group_by(cell) %>% 
  arrange(angle, .by_group = T) %>%
  mutate(rollsum = roll_sum(rep(total.int,2)[1:44], 
                            n = 9, 
                            align = "left",
                            fill = numeric(0))) %>%
  # %% is the modulo operation (gives the remainder of a division by 360).
  mutate(mean.roll.angle=(angle+40)%%360)

# Summarize intens2 to intens3 containing the maximum intensity at each 
# timepoint and the mean angle of the 90° rollsum with the highest total 
# intensity. Also the polarity is calculated as the portion of intensity of the 
# highest 90° rollsum and the total sum of intensities of the respective cell.
intens3 <- intens2 %>% 
  group_by(cell, timepoint) %>% 
  summarize(max.int = max(rollsum), 
            max.int.angle = mean.roll.angle[which.max(rollsum)], 
            X = unique(X), 
            Y = unique(Y), 
            timepoint = unique(timepoint), 
            polarity = max.int/sum(total.int))






### TRACKING INFORMATION

# This part can also be used on its own to analyze single cell tracking data.

# Import motility lab sheet from cell tracking. The files are named in the 
# pattern date_strain_clone_plate_light_replicate. The name will be broken at 
# the "_" to extract the respective information automatically from the file 
# name. Light is here typically given as a one letter color and the respective 
# photon flux, e.g. R75 for red and 75 µmol. 
# Example name: 20231026_pixG-eyfp_k3_1_R255_1_eyfp
# The tracking file should be in the same folder as the "info.csv" file
# otherwise a new working directory has to be set.
files <- list.files(pattern = "*.xlsx", full.names = T)

# If required read the *.xlsx file to correct for floating displacement that 
# contains the X and Y coordinates at timepoint 0 and the last timepoint, as 
# well as the name of the file to be corrected. Columns are:
# name	dx	dy	x1	x2	y1	y2
# where dx = x2 - x1 and dy = y2 - y1. 
# Usually the positions of at least 3 non-motile cells are manually measured per 
# video/file at the first and last timepoints to calculate the floating 
# distance dx and dy. 
# In order not to be read by the previous command as *.xlsx, the file is usually
# not saved in the same folder as the tracking file.
# Floating displacement happens if the agarose block to which the cells adhere 
# glides slightly on the glass bottom dish, thus uniformly biasing all tracks.
# The tracks are corrected for this uniform movement at a constant rate at each
# timepoint. The bias is particularly strong if the distances traveled by the 
# cells by twitching motility are relatively small in comparison to the float.
# It is not swimming of cells in a liquid stream which usually does not occur in 
# these experiments. If it does, the sample is not analyzed further.
float <- read_excel("//myr/home/jhammerl/Documents/microscopy/pixG-localization/20231026_pixG/analysis/float-correction.xlsx")

# Calculate the mean X and Y float per video and initialize a new column
# that will be used to find the index of the corresponding file name in the 
# "files" variable. This works also with many tracking files simultaneously, but
# here (i.e. together with the radial intensity) it is only done with one file
# at a time.
float <- float %>% 
  group_by(name) %>%
  summarize(dx=mean(dx),dy=mean(dy)) %>%
  mutate(f.index=NA)

# Find the index of the name from the "float" in the "files" variable.
for (k in float$name){
  nu <- match(T,str_detect(files,k))
  float$f.index[float$name == k] <- nu
}

# If there is no float correction available for a given name, the entry is
# filtered out.
float <- float %>% filter(!is.na(f.index))


# Initialize variables that will store all tracking information.
tracks <- NULL
tracks.total <- NULL

# Loop to read all Motility Lab tracking files (from TrackMate) in the working 
# directory.
for(i in files) {
  # Read files. The first row is omitted (it contains the column names) and 
  # columns are renamed.
  ML <- read_excel(i, col_names = c("track", "timepoint", 
                                    "time (sec)", "x", "y", "z"), 
                   skip=1)
  
  # NA tracks are omitted.
  ML <- ML %>% filter(!is.na(track))
  
  # Calculate cell density
  cd <- ML %>% 
    group_by(timepoint) %>%
    tally() %>%
    pull(n) %>%
    mean()
  
  # Group by track
  ML <- ML %>%
    group_by(track)
  
  # Correct X and Y positions for floating if present. If the respective name
  # is not indicated in "float", X and Y are not changed. Note: Even if a
  # particular file is not corrected "float" needs to be initialized. Otherwise
  # the next command has to be skipped not to throw the error 
  # "Error: object 'float' not found".
  # Subtract the per timepoint float distance from X and Y coordinates.
  if(match(i,files) %in% float$f.index){
    ML$x <- ML$x - (ML$timepoint*float$dx[float$f.index == match(i,files)]
                    /max(ML$timepoint))
    ML$y <- ML$y - (ML$timepoint*float$dy[float$f.index == match(i,files)]
                    /max(ML$timepoint))
  }
  
  # Calculate distance moved and angle (a) of the movement vector between 
  # consecutive frames for single tracks. The angles are converted from radian 
  # to degrees manually.
  ML <- ML %>%
    mutate(distance = sqrt((x-lag(x))^2 + (y-lag(y))^2),
           a = coord2rad(x-lag(x), y-lag(y))) %>%
    mutate(a = case_when(is.na(distance) ~ NA_real_,
                         distance == 0 ~ NA_real_,
                         distance > 0 ~ as.numeric(a))) %>%
    mutate(a = circular(a*180/pi, units="degrees"))
  
  
  
  # Summarize track properties. This is for the full-length tracks, i.e. over
  # all timepoints unlike the previous command which calculated between each two 
  # consecutive timepoints.
  ML.total <- ML %>%
    summarize("tracklength [frames]" = max(timepoint) - min (timepoint),
              "tracklength [sec]" = max(`time (sec)`) - min (`time (sec)`),
              "cell density" = cd,
              "distance" = sum(distance, na.rm = T),
              "displacement" =  sqrt((last(x)-first(x))^2 + 
                                       (last(y)-first(y))^2),
              "displacement angle" = coord2rad(last(x)-first(x), 
                                               last(y)-first(y))) %>%
    mutate("displacement angle" = circular(`displacement angle`*360/(2*pi), 
                                           units="degrees"),
           "speed" = distance/`tracklength [sec]`) %>%
    mutate("displacement angle" = abs(`displacement angle` - 360))
  

  
  
  # Filter tracks. The filter properties of f are then used to classify tracks
  # as valid or not, whereas the filter properties of f2 are used to classify
  # the tracks as either motile or non-motile. The filters may be adjusted 
  # according to the observed data and desired/necessary stringency.
  f <- ML.total %>% 
    filter(displacement >= 0) %>%
    filter(distance >= 0) %>%
    filter(`tracklength [frames]` >= 5)
  
  f2 <- ML.total %>% 
    filter(displacement >= 0) %>%
    filter(distance >= 0) %>%
    filter(speed >= 0.02)
  
  # Add the logical valid track and motile classifications to the single 
  # timepoint (ML) and summary (ML.total) data frames.
  ML <- ML %>% 
    mutate(valid.track = track %in% f$track, motile=track %in% f2$track)
  ML.total <- ML.total %>% 
    mutate(valid.track = track %in% f$track, motile=track %in% f2$track)
  
  
  # Add strain, clone, plate, replicate, color and intensity from the file name 
  # to ML and ML.total. The file name is broken at "_" characters. The light 
  # is split after the first character into color and intensity.
  strain = unlist(str_split(i, "_"))[2]
  clone = unlist(str_split(i, "_"))[3]
  plate = unlist(str_split(i, "_"))[4]
  replicate = unlist(str_split(i, "_"))[6]
  color = substring(unlist(str_split(i,"_"))[5],1,1)
  intensity = substring(unlist(str_split(i,"_"))[5],2)

  ML <- ML %>% 
    mutate(strain = strain, clone = clone, replicate = replicate, color = color, 
           intensity = intensity, plate = plate)
  ML.total <- ML.total %>% 
    mutate(strain = strain, clone = clone, replicate = replicate, color = color, 
           intensity = intensity, plate = plate)
  
  # Combine (multiple) tracking file(s) into one summary and one single 
  # timepoint resolved tracks object defined in the beginning before the loop.
  tracks <- bind_rows(tracks, ML)
  tracks.total <- bind_rows(tracks.total, ML.total)
  
}


# Group tracks
tracks <- tracks %>% group_by(strain, clone, replicate, track)
tracks.total <- tracks.total %>% group_by(strain, clone, replicate, track)

# The tracks can be saved as an *.RData object here. This is recommended in 
# particular if many tracking files have been processed simultaneously, as this
# can take some time in the for loop above. The object can then be loaded again
# when R is restarted without having to go through the loop again.
save(tracks, tracks.total, file="./tracks.RData")

# In the following the tracks are visualized as a lineplot in X and Y to give an 
# overview of movement directions. For ease of interpretation the tracks are 
# shifted into one common starting position.
p <- tracks %>% filter(valid.track==T) %>% mutate(x=x-first(x), y=y-first(y))

# Motile tracks are counted.
m_count <- tracks.total %>% group_by(strain, clone) %>% count(motile)

# Create a line plot in X and Y as a ggplot object.
p1 <- ggplot(p, aes(x = x, y = y)) +
  geom_path(aes(color = as.factor(track)), linewidth = 1.5, alpha = 0.3) +
  theme_light() +
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  coord_fixed(ratio=1) +
  labs(x = "x [µm]", y = "y [µm]")+
  scale_color_viridis(discrete = T, option="turbo") +
  theme_bw() +
  theme(legend.position="none")

# The above created plot can be split by different variables, either two 
# (facet_grid()) or only one (facet_wrap()). The number of motile and non-motile
# tracks is pasted in the plot.
p1 + facet_grid(clone ~ strain) +
  geom_text(data=filter(m_count, motile==T), aes(label=paste("motile:", n)), 
            x=-10, y=-7.5, size=2.5, hjust = 0) +
  geom_text(data=filter(m_count, motile==F), aes(label=paste("non-motile:", n)),
            x=-10, y=-9, size=2.5, hjust = 0)

p1 + facet_wrap(~clone)


# The following creates a summary object encompassing different statistic 
# measures for each strain.
# Group and filter valid and motile tracks 
p <- tracks.total %>% group_by(strain) %>% filter(valid.track==T)
tc <- p %>% summarize('ntracks' = length(track)) 
p <- p %>% filter(motile==T)

# Create a summary table for each strain. To assess the directionality, a 
# rayleigh test is performed. A specific direction alternative hypothesis, e.g.
# 0 degrees may be defined as:
## 'rayleigh test' = unlist(rayleigh.test(`displacement angle`, 
##                                        mu=circular(0, units="degrees"))[1]
# where [1] is the r statistic from the test. In the case of a defined mu, r 
# ranges from -1 to +1 where +1 is unidirectionality of all tracks at 0° and -1
# unidirectional movement at 180°. If no mu is defined, r ranges from 0 to +1
# where 0 is no directional bias and +1 is unidirectional movement in any 
# direction.
anno <- p %>% 
  summarize('rayleigh test' = unlist(rayleigh.test(`displacement angle`)[1]),
            'motile tracks' = sum(motile),
            'mean speed' = mean(speed),
            'mean distance' = mean(distance),
            'mean displacement' = mean(displacement),
            'SD speed' = sd(speed),
            'cell density' = mean(`cell density`)) %>% 
  left_join(tc) %>%
  mutate('non-motile tracks' = ntracks - `motile tracks`)

# Print the summary table and save it as a *.csv file to the working directory.
print(anno)
write.csv(anno, "results.csv")


# Plot windrose diagram (circular histogram) to visualize the distribution of 
# movement directions using the openair package. It is intended originally for 
# the visualization of wind directions (wd) and wind speeds (ws). Therefore, 
# speed and displacement angle are renamed so they are automatically recognized 
# by the windRose() command. The ws here is the single-cell speed but it may be 
# replaced by a different variable.
p <- p %>% transmute(ws=speed, wd=`displacement angle`)

p2 <- windRose(p, paddle=F,
               type = "strain",
               main="Frequencies of counts by movement direction (%)",
               ws.int = plyr::round_any(max(p$ws)/10, 0.01, f = ceiling),
               breaks = 11,
               key.header = "speed [µm/s]",
               key.footer = "[µm]",
               dig.lab = 2,
               key.position = "right",
               grid.line = list(value = 25, lty = 2, col = "grey"),
               cols = c("dodgerblue4", "white", "firebrick"),
               #border="black",
               annotate = F,
               #par.settings=list(fontsize=list(text=36))
               
)



# The following creates a windrose plot as above but with different colors
# indicating whether movement is towards or away from the directional light 
# source. 
#
# First p is reset, since it was changed for the previous windrose plot.
p <- tracks.total %>% 
  group_by(strain, color, intensity) %>% 
  filter(valid.track==T)
tc <- p %>% summarize('ntracks' = length(track)) 
p <- p %>% filter(motile == T)


p <- p %>% transmute(intensity = as.factor(intensity), 
                     ws = abs(`displacement angle`-180), 
                     wd = (`displacement angle`),
                     strain = strain)

# The following are filters that can be used to change the order of strains and
# intensities in the plot afterwards, but are not necessary. Strain names and 
# intensities are only exemplary.
#
## levels(p$strain) <- c("wt","pSK-pixG","pSK-pixG-DA","delta-pixD",
##                       "delta-pixD-x-pixG","delta-pixD-x-pixG-DA",
##                       "pSK-pixH","pSK-pixH-DA")
## levels(p$intensity)<- c("0","40","75")
## p <- p %>% filter(intensity == "40")
## p <- p %>% filter(color=="R",
##                   strain %in% c("WT","deltapixGH/pixGD326AH"))

# The following is to create a legend for the color scheme that indicates
# movement deviation from the incident light direction in 18° increments.
# The windRose() command requires it to be a list() to be accepted as a "key" 
# argument.
en <- seq(from=18, to =180, by=18)
st <- seq(from=0, to=162, by =18)
ll <- NULL 
for (i in 1:10){
  ll[i] <- print(paste0(st[i], " to ", en[i]))
}  
l <- list(labels=rev(ll))

# Plot windrose as above. Note that "type" is split into two variables to 
# create a 2D grid of plots by color and intensity, but it may as well be split
# by other variables such as strain. A maximum displayed frequency can be set
# which will change the proportions of the plots in the panels.
p2 <- windRose(p, paddle=F,
               type = c("color","intensity"),
               main="Frequencies of counts by movement direction (%)",
               ws.int = plyr::round_any(max(p$ws)/10, 1.0, f = ceiling),
               breaks = 11,
               key.header = "movement deviation from incident light",
               key.footer = "degree",
               #dig.lab = 2,
               #offset = 50,
               key.position = "bottom",
               grid.line = list(value = 30, lty = 2, col = "grey"),
               cols = c("grey30", "white", "darkgoldenrod2"),
               border="grey15",
               annotate = F,
               par.settings = list(fontsize=list(text=16)),
               key = l,
               #max.freq = 35
               
)



### COMBINATION OF RADIAL INTENSITY AND TRACKING DATA

# Identify corresponding cells in tracking and intensity data sets. 
# This is necessary because the two datasets come from two different plugins in 
# ImageJ, each of which does its own cell detection. Thus, the radial intensity 
# and tracking information have to be brought together manually here.
# The data sets are duplicated simply to retain the original data frames. 
tracks2 <- tracks
intens3.2 <- intens3

# An empty object into which the results will be printed afterwards
luest <- NULL

# Conceptually this is the same command as for the identification of closest 
# neighbors that has been performed above but here the query is the intensity
# data set and the target is the single cell tracking data set. 
#
# Note, that in the current state, it is possible that two cells from the 
# intens data frame are mapped onto the same track from the tracks data frame.
# Effectively, this overwrites the first matched intensity by the second matched
# one. This problems occurs for instance if the tracking plugin recognizes a 
# dividing cell as a single track but it is detected as two different cells 
# during the radial intensity measurements. This introduces some additional 
# noise currently.

for (i in 0:max(intens3.2$timepoint)) {
  query <- ungroup(intens3.2) %>% filter(timepoint == i)
  target <- ungroup(tracks2) %>% filter(timepoint == i)
  neigh <- nn2(target[,c("x","y")], query[,c("X","Y")],
               k = 1, 
               searchtype = "radius", 
               radius=2.5)$nn.idx
  neigh[neigh == 0] <- NA
  to.append <- data.frame(track = target$track[neigh], 
                          cell = query$cell, 
                          max.int = query$max.int, 
                          polarity = query$polarity, 
                          max.int.angle = query$max.int.angle,
                          timepoint = query$timepoint)
  luest <- bind_rows(luest, to.append)
}

tracks2 <- left_join(tracks2,luest,by=c("timepoint","track"))

# Omit tracks that have no matching in the intensity data, e.g. tracks that are
# on the edge of the image are excluded in the intensity measurement but not the
# single cell tracking.
tracks2 <- tracks2[!is.na(tracks2$max.int.angle),]


# This next section is not strictly necessary, it only serves as a means to 
# visualize the data that has been produced so far and look at it graphically as 
# to assess whether it seems sensible or if some parameters may need adjustment  
# in the analysis. It will create an image that is one slice of the original 
# time lapse series overlaid with the recognized cells and their respective 
# highest fluorescence intensity angle.
# 
# Load the .jpeg files of each single slice of the stack that have been saved
# before using ImageJ. It is useful to save the jpeg frames from ImageJ using a 
# LUT that allows easy differentiation of intensity, e.g. physics and adjusting 
# the contrast.
l <- list.files("./jpeg",pattern = "*.jpg",full.names = T)

# Read one .jpeg. timepoint 0 should correspond to l[1], timepoint 1 to l[2] and
# so on if all slices are saved in the folder.
t1 <- readJPEG(l[1])

# Create the raster image graphical object to use for the plot from the .jpeg. 
# Basically this is the background onto which everything else will be plotted.
t2 <- grid::rasterGrob(t1, 
                       width = unit(1,"npc"), 
                       height = unit(1,"npc"), 
                       interpolate = FALSE)

# Create a subset of the original intens table at the selected timepoint (match
# with the selected .jpeg) that will be used to draw circles indicating the ring
# between 30% and 75% radius. As before the radius is converted into 
# micrometers with the µm/pixel ratio from ImageJ. The subset intg is then 
# summarized to intg1 containing the respective radii.
intg <- intens %>% 
  filter(timepoint == 0) %>% 
  group_by(cell, angle) %>% 
  mutate(radius = radius * 0.0965332)

intg1 <- intg %>% 
  summarize(X=unique(X), 
            Y=unique(Y), 
            large = quantile(radius, .75), 
            small = quantile(radius, .3))

# Add the maximum intensity angle from the intens3 data set.
intg2 <- left_join(intg1,
                   dplyr::select(intens3,c("cell","max.int.angle")),by="cell")

# Create the plot on the jpeg with indicated directions
tracks2 %>% ggplot(aes(x=x, y=-y)) +  
  # This is the jpeg background
  annotation_custom(t2, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  # This is the single cell tracking. All positions of the tracks are displayed 
  # as color coded dots
  geom_point(data=tracks, aes(x=x, y=-y), 
             color = as.factor(tracks$track), alpha=0.5) +
  # This is the 30% radius line, that goes with an arbitrary trial and error 
  # scaling factor to match the respective line as observed in the radial angle
  # profile plugin because the conversion to micrometers does not scale 
  # correctly 
  geom_point(data = intg1, aes(x = X, y = -Y), 
             pch = 21, size = 7.8*intg1$small,color="white") + 
  # 75% radius
  geom_point(data = intg1, aes(x = X, y = -Y), 
             pch = 21, size = 7.8*intg1$large,color="white") +
  # Add an arrow indicating the maximum fluorescence angle
  geom_spoke(data = intg2,  aes(x = X, y=-Y, 
                                angle = max.int.angle*pi/180, radius = 2.7 ), 
             arrow = arrow(length = unit(0.2,"cm")),linewidth=1,color="white") +
  # Add two lines to indicate the 90° sector which has been assigned the highest
  # total intensity
  geom_spoke(data = intg2, aes(x = X, y = -Y, 
                               angle = (max.int.angle+45)*pi/180, radius = 2.7), 
             color = "grey", linewidth = 0.7) +
  geom_spoke(data = intg2, aes(x = X, y = -Y, 
                               angle = (max.int.angle-45)*pi/180, radius = 2.7), 
             color = "grey", linewidth = 0.7) +
  # Scale the plot so that the objects align correctly onto the jpeg. This is 
  # also, inelegantly, derived from a trial and error process. Perhaps there is 
  # an easy way to do it, I did not succeed in figuring it out here.
  coord_fixed(xlim=c(2.2,47.22),ylim=c(-47.22,-2.2)) +
  # The size of the image in the background.
  scale_x_continuous(breaks = c(0,49.42)) + 
  scale_y_continuous(breaks = c(0,-49.42))

  

# This last section contains a lot of rotating operations for the calculated 
# angles. The purpose is to match the later created circular histograms to the 
# actual microscopy images (e.g. the cells move down in the microscopy image but  
# right in the plot). This has different reasons, for example that
# different packages in R and different plugins in ImageJ use different starting
# angles (e.g. 0 degrees are up or toward the right) or rotate angles in 
# different directions (counter clockwise or clockwise), which is contributed by
# the fact that ImageJ has the coordinate system origin at the top left of the 
# image.
# The transformations of the angles have been carefully checked using 
# minimalist mock images (black images with single brighter circles with 
# single white spots in defined directions) that were subjected to the entire 
# analysis. 
# 
# The mean fluorescence angle is averaged over different time spans to see 
# whether bleaching has a strong effect on the distribution plotted afterwards.
#
# The maximum timepoints to which the maximum intensity angles should be 
# averaged
times <- c(5, 10, 20, 30, 60)

# Initialize the objects to contain the results later depending on the cutoff 
# timepoint. They are only initialized the first time and need to be reset 
# manually when switching to a different strain. The way the analysis was done
# here is that for each strain and condition (e.g. pixG-eyfp in the dark) all
# videos were analyzed to create the tracks.total.full objects for given times.
# Then they were saved and deleted from the environment before proceeding to 
# the next strain/condition. 
for (v in 1:length(times)) {
  if(exists("tracks.total.full60") == F){
  eval(parse(text = paste('tracks.total.full', times,' = NULL', sep ='')[v]))
  }
}

# The looping variable for the following "for()" loop
if(exists("da") == F){
  da <- c( "tracks.total.full60",
           "tracks.total.full10",
           "tracks.total.full20",
           "tracks.total.full30",
           "tracks.total.full5")
}

for(k in c(da)){
  # Rotate and mirror the mean intensity angles to match them to the microscopy 
  # images. The angles over the time determined by "da" are averaged as they 
  # fluctuate over timepoints although the general direction of the maximum  
  # fluorescence does not show obvious rapid changes. 
  # Also, the polarity in the first frame is calculated to give an idea of it 
  # before any bleaching.

  tracks2.test <- tracks2 %>% 
    filter(timepoint < unlist(str_split(k, "full"))[2]) %>% 
    group_by(track) %>% 
    summarize(mean.intens.angle=(mean(max.int.angle-90)+360)%%360,
              mean.polarity=mean(polarity),
              polarity_start=first(polarity)) %>%
    mutate(mirror_fl = (180 - ((mean.intens.angle)%%180))*2) %>% 
    group_by(track)%>%
    mutate(mean.intens.angle=sum(c(mean.intens.angle,mirror_fl))%%360) %>%
    mutate(mean.intens.angle=(mean.intens.angle))
  # Bind the newly calculated values to tracks.total2
  tracks.total2 <- left_join(tracks.total,tracks2.test, by = "track")
  # Same angle transformations as before for the intensity angles but for the
  # displacement angles. No mean is calculated here as the displacement angle is
  # already the sum of all single movements made between all timepoints, which
  # is essentially an average of the direction weighted by the step size.
  tracks.total2 <- tracks.total2 %>% 
    filter(!is.na(mean.intens.angle)) %>% 
    mutate(`displacement angle`= ((`displacement angle`-90)+360)%%360) %>%
    mutate(mirror_movement = (180 - ((`displacement angle`)%%180))*2) %>% 
    group_by(track) %>%
    mutate(`displacement angle`= sum(c(`displacement angle`,
                                       mirror_movement))%%360) %>%
    mutate(mean.intens.angle=(mean.intens.angle)) %>% 
    # Calculate the angular deviation as the difference between displacement 
    # and highest fluorescence angle such that 0° means the both angles are 
    # equal. 
    mutate(ang.dev=((`displacement angle`- mean.intens.angle)+360)%%360)
  
  # If there is already data saved in the tracks.total.full object, the track
  # number is changed so that there are not two different entries with the same
  # track index.
  if(length(get(k))>0){
    tracks.total2$track <- tracks.total2$track+(max(get(k)$track)+1)
  }
  # The above calculated data is stored in the tracks.total.full object (k)
  assign(k,bind_rows(get(k),tracks.total2)) 
}

# This is the data only for the current replicate.
tracks.total.full_m3 <- tracks.total2

### At this point you can return to the very beginning of the script and read 
### the data from the next video. In the end the calculated angles will be 
### appended to the tracks.total.full data frames. After all replicates of 
### one strain and condition have been compiled, I usually save the image of the
### current environment so I can retrieve all the variables created for this 
### particular strain and condition later again, before continuing to the next
### condition. You can however also continue with single replicates and look at
### the distributions of single replicates to identify if single replicates
### show any weird anomalies.



# Save either the entire environment (save.image()) or only single objects
# (Save()). Just examplary paths.
## save(tracks.total.full, file="//myr/home/jhammerl/Documents/microscopy/pixG-localization/20231026_pixG/analysis/tracks-total-full.RData")
## save.image(file="//myr/home/jhammerl/Documents/microscopy/pixG-localization/20231026_pixG/analysis/20231026_pixG-eyfp_k3_1_R255_3_eyfp/G255_environment.RData")



# The following steps, including the very long loop, serve the purpose to   
# manually bin the data for visualization as a circular histogram as there were  
# problems binning over the 0° angle when done with ggplot automatically. 
# Unfortunately, I did not succeed in solving the problems and thus resorted to 
# manual binning for now. The code to plot the same data in ggplot is appended 
# at the end. The plot looks the same as the one from manual binning but the bin 
# over 0° cannot be displayed.
#
# The data to be binned. Here it is tracks.total.full5 but it can be changed to 
# any of the other tracks.total.full data objects.
dog <- tracks.total.full60 %>%
  filter(valid.track==T) 
# Bins for the polarity.
pony <- seq(0.25,0.75,0.025)

# The data frame to contain the bin counts afterwards. "angle.bins" is the mean
# angle of the respective bin that is repeated for the number of possible 
# polarity bins. The classifications to be binned into are the displacement 
# angle (move.cnts), the angle of highest fluorescence intensity (fluor.cnts),
# and the angular deviation of displacement and fluorescence angle. 
hist.df <- data.frame(angle.bins = rep(seq(from = 0, to = 345,by = 15),
                                       each = length(pony)-1),
                      polarity.start = rep(pony[-length(pony)],24),
                      polarity.end = rep(pony[2:length(pony)],24),
                      move.cnts = NA,
                      fluor.cnts = NA,
                      ang.dev.cnts = NA)
  
seq <- seq(from = 7.5, to = 352.5, by = 15)

# Loop over the angles between which to bin (i) and the polarity cutoffs (j), 
# the latter of which are currently not displayed in the plot afterwards.
# The dummy is length()-1 because the lower border is used for the binning, thus
# the dummy does not take the last value in the sequence itself.
for(i in 1:(length(seq)-1)){
  for (j in 1:(length(pony)-1)) {
    # Binning displacement angles 
    mc <- sum(dog$`displacement angle` > seq[i] & 
              dog$`displacement angle` <= seq[i+1] & 
              dog$mean.polarity >= pony[j] & 
              dog$mean.polarity < pony[j+1])
    hist.df[hist.df$angle.bins==mean(c(seq[i],seq[i+1])) & 
            hist.df$polarity.start == pony[j],]$move.cnts <- mc
    # This if() is to bin over 0° if the polarity dummy is not at its max.
    if(i==length(seq)-1 & j!=length(pony)-1){
      hist.df[hist.df$angle.bins==0 & 
              hist.df$polarity.start==pony[j],]$move.cnts <- 
        sum((dog$`displacement angle`>seq[i+1] & 
             dog$mean.polarity >= pony[j] & 
             dog$mean.polarity < pony[j+1]) | 
            (dog$`displacement angle` <= seq[1] & 
             dog$mean.polarity >= pony[j] & 
             dog$mean.polarity < pony[j+1]))
    }
    # This if() is to count in the case where i and j reach the maximum value, 
    # i.e. binning over the 0° border at the highest possible polarity.
    if(i==length(seq)-1 & j==length(pony)-1){
      hist.df[hist.df$angle.bins==0 & 
              hist.df$polarity.start==pony[j],]$move.cnts <- 
        sum((dog$`displacement angle`>seq[i+1] & 
             dog$mean.polarity >= pony[j] & 
             dog$mean.polarity <= pony[j+1]) | 
            (dog$`displacement angle` <= seq[1] & 
             dog$mean.polarity >= pony[j] & 
             dog$mean.polarity <= pony[j+1]))
    }
    # Binning fluorescence angles; same principles as for displacement.
    fc <- sum(dog$mean.intens.angle>seq[i] & 
              dog$mean.intens.angle <= seq[i+1] & 
              dog$mean.polarity >= pony[j] & 
              dog$mean.polarity < pony[j+1])
    hist.df[hist.df$angle.bins==mean(c(seq[i],seq[i+1]))& 
            hist.df$polarity.start == pony[j],]$fluor.cnts <- fc
    
    if(i==length(seq)-1 & j!=length(pony)-1){
      hist.df[hist.df$angle.bins==0 & 
              hist.df$polarity.start==pony[j],]$fluor.cnts <- 
        sum((dog$mean.intens.angle>seq[i+1] & 
              dog$mean.polarity >= pony[j] & 
              dog$mean.polarity < pony[j+1]) | 
            (dog$mean.intens.angle <= seq[1] & 
              dog$mean.polarity >= pony[j] & 
              dog$mean.polarity < pony[j+1]))
    }
    if(i==length(seq)-1 & j==length(pony)-1){
      hist.df[hist.df$angle.bins==0 & 
              hist.df$polarity.start==pony[j],]$fluor.cnts <- 
        sum((dog$mean.intens.angle>seq[i+1] & 
             dog$mean.polarity >= pony[j] &
             dog$mean.polarity < pony[j+1]) |
            (dog$mean.intens.angle <= seq[1] & 
             dog$mean.polarity >= pony[j] & 
             dog$mean.polarity <= pony[j+1]))
    }
    # Binning angular deviation angles; same principles as for displacement.
    ac <- sum(dog$ang.dev > seq[i] &
              dog$ang.dev <= seq[i+1] & 
              dog$mean.polarity >= pony[j] & 
              dog$mean.polarity < pony[j+1])
    
    hist.df[hist.df$angle.bins==mean(c(seq[i],seq[i+1])) & 
            hist.df$polarity.start == pony[j],]$ang.dev.cnts <- ac
    if(i==length(seq)-1 & j!=length(pony)-1){
      hist.df[hist.df$angle.bins==0 & 
              hist.df$polarity.start==pony[j],]$ang.dev.cnts <- 
        sum((dog$ang.dev>seq[i+1] & 
             dog$mean.polarity >= pony[j] & 
             dog$mean.polarity < pony[j+1]) |
            (dog$ang.dev <= seq[1] & 
             dog$mean.polarity >= pony[j] & 
             dog$mean.polarity < pony[j+1]))
    }
    if(i==length(seq)-1 & j==length(pony)-1){
      hist.df[hist.df$angle.bins==0 & 
              hist.df$polarity.start==pony[j],]$ang.dev.cnts <- 
        sum((dog$ang.dev>seq[i+1] & 
             dog$mean.polarity >= pony[j] & 
             dog$mean.polarity < pony[j+1]) |
            (dog$ang.dev <= seq[1] & 
             dog$mean.polarity >= pony[j] & 
             dog$mean.polarity <= pony[j+1] ))
    }
  }
}

# Calculate percentages from the counts binned above.
hist.df <- hist.df %>% 
  mutate(move.per = move.cnts/sum(move.cnts), 
         fluor.per = fluor.cnts/sum(fluor.cnts),
         ang.dev.per = ang.dev.cnts/sum(ang.dev.cnts)) %>% 
  group_by(angle.bins)

# Plot a circular histogram overlaying the displacement angles and the mean 
# maximum intensity angles. It is possible depending on the data, that bins are
# omitted because they exceed the set y scale limits. In such a case, the limits
# need to be adjusted. 
ggplot(data = hist.df %>% group_by(angle.bins))+
  # Polar coordinates to make the histogram circular.
  coord_polar(theta = "x", start = -0.1305) +
  # X and Y axis definition
  scale_x_continuous(breaks = seq(from = 0, to = 360, by = 15), 
                     expand = c(.002, 0)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2),
                     limits = c(0, 0.2))+
  # The displacement angle histogram (barplot with manual bins).
  geom_bar(data=hist.df %>% summarize(move.per = sum(move.per)),
           stat = "identity", 
           fill = "salmon1",
           aes(x = angle.bins, y = move.per),
           alpha = 0.9, 
           colour = "black", 
           size = 0.75) +
  # The fluorescence angle histogram.
  geom_bar(data=hist.df %>% summarize(fluor.per = sum(fluor.per)),
           stat = "identity",
           fill = "aquamarine2",
           aes(x = angle.bins,y = fluor.per),
           alpha = 0.65,
           colour = "black", 
           size= 0.75) +
  theme_bw() +
  theme(panel.grid = element_line(color = "grey75",size = 0.75,linetype = 1),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20), 
        axis.text.x = element_text(hjust = 10),
        panel.grid.major = element_line(linewidth = 0.5))  


# For each strain and replicate the final tracks.total.full5 data frame is
# stored under the respective new name so that afterwards the data is available 
# for additional plots or statistical analysis. The naming is for example: 
# pol_DA_R255: pixG-DA-eyfp with directional red light.
pol_DA_R255 <- tracks.total.full5
pol_DA_R0 <- tracks.total.full5
pol_G_R255 <- tracks.total.full5
pol_G_R0 <- tracks.total.full5

# Here certain replicates where omitted due to floating of the agarose block, 
# and thus the cells, that could not be well corrected with the implemented  
# floating correction during the single cell tracking analysis. It can be done 
# like this:
pol_DA_R0.f <- pol_DA_R0 %>% 
  filter(valid.track==T, !replicate %in% as.character(c(2,3,6))) 
pol_DA_R255.f <- pol_DA_R255 %>% 
  filter(valid.track==T) 
pol_G_R0.f <- pol_G_R0 %>% 
  filter(valid.track==T, replicate != "10") 
pol_G_R255.f <- pol_G_R255 %>% 
  filter(valid.track==T, !replicate %in% as.character(c(1,6))) 

# Bind the data into a single final data frame 
all.data.final <- bind_rows(pol_DA_R0.f,pol_DA_R255.f)
all.data.final <- bind_rows(all.data.final,pol_G_R0.f)
all.data.final <- bind_rows(all.data.final,pol_G_R255.f)

# Watson two sample test to test circular distributions. Below an example for 
# the displacement angle distribution of pixG-DA-eyfp in light versus in dark.
watson.two.test(circular(pol_DA_R0.f$`displacement angle`,units = "degrees"), 
                circular(pol_DA_R255.f$`displacement angle`, units = "degrees"))

# Compare overlap of circular distributions using the activity package. Below an
# example for pixG-eyfp under directional light with 999 bootstrap iterations.
compareCkern(fitact(pol_G_R255.f$`displacement angle`*(pi/180)),
             fitact(pol_G_R255.f$mean.intens.angle*(pi/180)), 
             reps = 999)



# The ggplot version of the circular histogram, without manual binning
# Note that the bin that spans 0° is missing in comparison to the manually 
# binned version
ggplot(tracks.total.full60%>%filter(valid.track==T), 
       aes(`displacement angle`)) +
  geom_histogram(aes(y = (..count..)/sum(..count..),fill="1"), 
                 alpha=0.9,
                 boundary= -7.5,
                 bins= 24, 
                 colour = "black", 
                 size = 0.75) + 
  geom_histogram(aes(x = mean.intens.angle,y = (..count..)/sum(..count..),
                     fill= "2"),
                 alpha=0.65,
                 colour = "black",
                 boundary= -7.5, 
                 bins = 24, 
                 size= 0.75)+
  coord_polar() +
  scale_x_continuous(limits=c(0,360),breaks = seq(0,360,by=45))+
  theme_light()+
  scale_fill_manual(values=c("salmon1","aquamarine2"), 
                    labels = c("movement","fluorescence")) +
  theme(panel.grid.minor = element_line(linewidth = 0.5), 
        panel.grid.major = element_line(linewidth = 1),
        text = element_text(size=20), 
        axis.text.x = element_text(hjust = 10)) #+
# Can be added to compare different replicates
## facet_wrap(~as.integer(replicate))