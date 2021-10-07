
# Setup -------------------------------------------------------------------
require(dplyr)
require(ggplot2)
require(lme4)
require(lmerTest)

theme_update(strip.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             plot.background = element_blank(),
             axis.text=element_text(size=30, colour = "black"),
             axis.title=element_text(size=18, face = "bold"),
             axis.title.x = element_text(vjust = 0),
             legend.title = element_text(size = 24, face = "bold"),
             legend.text = element_text(size = 20),
             plot.title = element_text(size = 26, face = "bold", vjust = 1),
             panel.margin = unit(1.0, "lines"), 
             plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
             axis.line = element_line(colour = "black", size = 2),
             axis.ticks = element_line(color = 'black', size = 3),
             axis.ticks.length = unit(.25, 'cm')
)

theme_black = function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "white"),  
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
    )
}
se = function(x) {return(sd(x, na.rm = T) / sqrt(sum(!is.na(x))))}
as.string.vector = function(x) {
  return(strsplit(x,',')[[1]])
}
as.string = function(x) {
  return(paste(x, collapse = ','))
}
dodge <- position_dodge(width=0.9)

# Only works in RStudio -- otherwise you have to set the path manually
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data ---------------------------------------------------------------
df.demo.raw = read.csv('demo.csv', stringsAsFactors = F) %>% arrange(subject) %>% mutate(total_time_real = total_time / 60000)
df.s1.raw = read.csv('s1.csv', stringsAsFactors = F) %>% arrange(subject)
df.s2.raw = read.csv('s2.csv', stringsAsFactors = F) %>% arrange(subject)


# Do filtering ---------------------------------------------------------

df.demo = df.demo.raw
df.s1 = df.s1.raw %>% filter(practice == 0) %>% mutate(subject = factor(subject))
df.s2 = df.s2.raw

# Check out data ----------------------------------------------------------

hist(df.demo$total_time_real)

df.s1.subj = df.s1 %>% group_by(subject) %>%
  summarize(total.time = sum(rt) / 60000)
hist(df.s1.subj$total.time)

## Stage 1 choices
atts = c('Number of Bedrooms','Size of Garage','Amount of Crime in Neighborhood','Proximity to Parks','Proximity to Waterfront/Beaches',
'Proximity to Cafes/Restaurants','Noise Pollution','Reputation of Closest School','Amount of Natural Light','Age of Building','Washer/Dryer','Size of Yard','Fireplace',
'Central AC','Climate of Area','Hardwood Floors','Freshly Painted Exterior','Size of Home')

atts.opt1 = paste0(atts,'.opt1')
atts.opt2 = paste0(atts,'.opt2')
atts.opt1.enclosed = paste0('`',atts.opt1,'`')
atts.opt2.enclosed = paste0('`',atts.opt2,'`')

df.s1[,atts.opt1] = NA
df.s1[,atts.opt2] = NA

for (i in 1:nrow(df.s1)) {
  cur.atts = as.string.vector(df.s1$attributes[i])
  cur.opt1.vals = as.string.vector(df.s1$opt1_values[i])
  cur.opt2.vals = as.string.vector(df.s1$opt2_values[i])
  
  #att.nums = numeric(length(cur.atts))
  for (j in 1:length(cur.atts)) {
    #att.nums[j] = which(cur.atts[j] == atts)
    df.s1[i,paste0(cur.atts[j],'.opt1')] = cur.opt1.vals[j]
    df.s1[i,paste0(cur.atts[j],'.opt2')] = cur.opt2.vals[j]
  }
  
  #df.s1$att.nums[i] = as.string(att.nums)
}

# clean up data
scale1 = c('None', 'Small', 'Medium', 'Large')
scale2 = c('Very Low', 'Low', 'Moderate', 'High', 'Very High')
scale3 = c('Very Bad', 'Bad', 'Moderate', 'Good', 'Very Good')
scale4 = c('Old', 'Medium', 'New')
scale5 = c('No', 'Yes')
scale6 = c('No natural light', 'A little natural light', 'Moderate natural light', 'A lot of natural light')
scale7 = c('Low', 'Medium', 'High')

for (i in 1:length(atts)) {
  cur.att.opt1 = paste0(atts[i],'.opt1')
  cur.att.opt2 = paste0(atts[i],'.opt2')
  
  if (atts[i] %in% c('Number of Bedrooms')) {
    df.s1[,cur.att.opt1] = as.numeric(df.s1[,cur.att.opt1])
    df.s1[,cur.att.opt2] = as.numeric(df.s1[,cur.att.opt2])
  } else if (atts[i] %in% c('Size of Home', 'Proximity to Parks', 'Proximity to Waterfront/Beaches', 'Proximity to Cafes/Restaurants')) {
    df.s1[,cur.att.opt1] = as.numeric(sub('\\ .*', '', df.s1[,cur.att.opt1]))
    df.s1[,cur.att.opt2] = as.numeric(sub('\\ .*', '', df.s1[,cur.att.opt2]))
  } else if (atts[i] %in% c('Amount of Crime in Neighborhood')) {
    df.s1[,cur.att.opt1] = as.numeric(factor(df.s1[,cur.att.opt1], scale2))
    df.s1[,cur.att.opt2] = as.numeric(factor(df.s1[,cur.att.opt2], scale2))
  } else if (atts[i] %in% c('Central AC', 'Fireplace', 'Hardwood Floors', 'Freshly Painted Exterior', 'Washer/Dryer')) {
    df.s1[,cur.att.opt1] = factor(df.s1[,cur.att.opt1], scale5)
    df.s1[,cur.att.opt2] = factor(df.s1[,cur.att.opt2], scale5)
  } else if (atts[i] %in% c('Size of Garage', 'Size of Yard')) {
    df.s1[,cur.att.opt1] = as.numeric(factor(df.s1[,cur.att.opt1], scale1))
    df.s1[,cur.att.opt2] = as.numeric(factor(df.s1[,cur.att.opt2], scale1))
  } else if (atts[i] %in% c('Age of Building')) {
    df.s1[,cur.att.opt1] = as.numeric(factor(df.s1[,cur.att.opt1], scale4))
    df.s1[,cur.att.opt2] = as.numeric(factor(df.s1[,cur.att.opt2], scale4))
  } else if (atts[i] %in% c('Amount of Natural Light')) {
    df.s1[,cur.att.opt1] = as.numeric(factor(df.s1[,cur.att.opt1], scale6))
    df.s1[,cur.att.opt2] = as.numeric(factor(df.s1[,cur.att.opt2], scale6))
  } else if (atts[i] %in% c('Reputation of Closest School', 'Climate of Area')) {
    df.s1[,cur.att.opt1] = as.numeric(factor(df.s1[,cur.att.opt1], scale3))
    df.s1[,cur.att.opt2] = as.numeric(factor(df.s1[,cur.att.opt2], scale3))
  } else if (atts[i] %in% c('Noise Pollution')) {
    df.s1[,cur.att.opt1] = as.numeric(factor(df.s1[,cur.att.opt1], scale7))
    df.s1[,cur.att.opt2] = as.numeric(factor(df.s1[,cur.att.opt2], scale7))
  }
}

# run analysis!
formula = paste0('choice ~ ',
                 paste(atts.opt1.enclosed, collapse = " + "),
                 #' + ',
                 #paste(atts.opt2.enclosed, collapse = " + "),
                 ' + (1 | subject)')
m1 = glmer(as.formula(formula), data = df.s1, family = 'binomial')