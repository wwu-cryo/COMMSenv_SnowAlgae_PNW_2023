#### pigment analysis ####
rm(list=ls())

library(ggplot2)
library(ggpubr)
library(readr)
library(lubridate)
library(tidyr)

#import data
setwd("~/OneDrive - Western Washington University/Thesis/Manuscripts/Bagley UAV Paper")

bb <- read.csv("All Bagley 2021 data.csv")
bb$Date <- as.POSIXct(as.character(bb$Date), format = "%m/%d/%Y")
bb$Date <- as.Date(bb$Date)
bb$Date <- bb$Date %m+% years(2000)

bb<-dplyr::select(bb, -dplyr::contains("RSD"))
bb <- bb[1:54,1:40]


names(bb)[names(bb) == "Astaxanthin.Conc...ug.L."] <- "astaxanthin"
names(bb)[names(bb) == "Lutein.Conc...ug.L."] <- "lutein"
names(bb)[names(bb) == "Chl.B.Conc...ug.L."] <- "chl.b"
names(bb)[names(bb) == "Chl.A.Conc...ug.L."] <- "chl.a"
names(bb)[names(bb) == "Beta.Carotene..ug.L."] <- "bc"
names(bb)[names(bb) == "algae.concentration..cells.mL."] <- "cell"

options(digits = 9)

bb$astaxanthin <- as.character(bb$astaxanthin)
bb$lutein <- as.character(bb$lutein)
bb$chl.b <- as.character(bb$chl.b)
bb$chl.a <- as.character(bb$chl.a)
bb$bc <- as.character(bb$bc)

bb$astaxanthin <- as.numeric(bb$astaxanthin)
bb$lutein <- as.numeric(bb$lutein)
bb$chl.b <- as.numeric(bb$chl.b)
bb$chl.a <- as.numeric(bb$chl.a)
bb$bc <- as.numeric(bb$bc)
bb <- bb[1:54,1:13]


bb$ast.chla <- bb$astaxanthin/bb$chl.a
bb$lut.chla <- bb$lutein/bb$chl.a
bb$chlb.chla <- bb$chl.b/bb$chl.a
bb$bc.chla <- bb$bc/bb$chl.a

bb.702 <- subset(bb, bb$Date=="2021-07-02")
bb.730 <- subset(bb, bb$Date=="2021-07-30")

#bb.mica <- rbind(bb.702,bb.730)
bb.mica.long.702 <- gather(bb.702, pigment, concentration, ast.chla:bc.chla, factor_key=TRUE)
bb.mica.long.730 <- gather(bb.730, pigment, concentration, ast.chla:bc.chla, factor_key=TRUE)
bb.mica.long <- rbind(bb.mica.long.702,bb.mica.long.730)

tapply(bb.mica.long.702$concentration, bb.mica.long.702$pigment, max)
tapply(bb.mica.long.730$concentration, bb.mica.long.730$pigment, max)
tapply(bb.mica.long$concentration, bb.mica.long$pigment, mean)

bb.mica.long$Date <- as.factor(bb.mica.long$Date)

ggplot(data=bb.mica.long, aes(x=pigment, y=concentration))+
  geom_boxplot(aes(fill = Date))+
  xlab("") + ylab(expression(paste("Ratio to Chlorophyll ",italic("a"))))+
  scale_x_discrete(labels=c("Astaxanthin", "Lutein",expression(paste("Chlorophyll ", italic("b"))), "Beta Carotene"))+
  scale_fill_manual(labels=c("2 July", "30 July"),
                    values = c("2021-07-02"="darksalmon", "2021-07-30"="darkturquoise")) +
  theme_minimal(base_size = 20)

ggsave(filename = "pigment_ratios.jpg", device='jpg', dpi=2100)

# plot pigments over time
par(mfrow=c(2,2))
cell <- ggplot(data=bb, aes(x=factor(Date),y=cell)) +
  geom_boxplot() +
  labs(title= "Cell Concentration",
       x="",
       y="Concentration (cells/mL)",
       scale_x_discrete(guide = guide_axis(angle = 90))) +
  #ylim(0,60)+
  theme_classic()

ast <- ggplot(data=bb, aes(x=factor(Date),y=astaxanthin)) +
  geom_boxplot() +
  labs(title= "Astaxanthin",
       x="",
       y=expression(paste("Concentration (", mu,"g/L)")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ylim(0,60)+
  theme_classic()

lut <- ggplot(data=bb, aes(x=factor(Date),y=lutein)) +
  geom_boxplot() +
  labs(title= "Lutein",
       x="",
       y=expression(paste("Concentration (", mu,"g/L)")))+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ylim(0,60)+
  theme_classic()

chl.b<- ggplot(data=bb, aes(x=factor(Date),y=chl.b)) +
  geom_boxplot() +
  labs(title= "Chlorophyll B",
       x="Date",
       y=expression(paste("Concentration (", mu,"g/L)")))+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylim(0,60)+
  theme_classic()

chl.a <- ggplot(data=bb, aes(x=factor(Date),y=chl.a)) +
  geom_boxplot() +
  labs(title= "Chlorophyll A",
       x="Date",
       y=expression(paste("Concentration (", mu,"g/L)")))+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylim(0,60)+
  theme_classic()

bc <- ggplot(data=bb, aes(x=factor(Date),y=bc)) +
  geom_boxplot() +
  labs(title= "Beta Carotene",
       x="Date",
       y=expression(paste("Concentration (", mu,"g/L)")))+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylim(0,60)+
  theme_classic()

#create figure 1
ggarrange(cell+ rremove("x.text"), ast+ rremove("x.text"), lut+ rremove("x.text"), 
          chl.a, chl.b, bc, 
          labels = c("a", "b", "c", "d", "e","f"),
          ncol = 3, nrow = 2)

ggsave(filename = "fig_1.png", device="png", dpi=500, width = 8, height = 5.4, units = "in")

#### PCA Analysis ####
setwd("~/OneDrive - Western Washington University/Thesis/Manuscripts/Bagley UAV Paper/PCA")

PCA.702 <- read.csv("PCA_20210702_accuracy_points.csv")
PCA.730 <- read.csv("PCA_20210730_accuracy_points.csv")

##### 20210702 #####
# plot PC1 vs PC2
my.plot.df <- data.frame(
  class = as.factor(PCA.702$class),
  cell = PCA.702$cell.concentration..cells.mL.,
  cell.ln =log(PCA.702$cell.concentration..cells.mL.),
  PC1 = PCA.702$b1_PCA_202,
  PC2 = PCA.702$b2_PCA_202,
  PC3 = PCA.702$b3_PCA_202
)

levels(my.plot.df$class)
levels(my.plot.df$class) <- c(levels(my.plot.df$class), "other")

my.plot.df["class"][my.plot.df["class"] == "rock"] <- "other"
my.plot.df["class"][my.plot.df["class"] == "water"] <- "other"
my.plot.df["class"][my.plot.df["class"] == "vegetation"] <- "other"

summary(my.plot.df)

# make ggplot  
p1 <- ggplot(my.plot.df, aes(x = PC1, y = PC2, color=class, shape=class)) +
  geom_point(size = 2) +
  geom_hline(yintercept = -2500, linetype = "dashed") +
  geom_vline(xintercept = -90000, linetype="dashed") +
  #geom_abline(aes(intercept = 2500, slope =.1), linetype="dashed") +
  scale_color_manual(name="Class", labels = c("Snow", "Snow Algae", "Other"), 
                     values = c("snow" = "darkturquoise",
                                "snow algae"="darksalmon","other"="black")) +
  scale_shape_manual(name="Class", labels = c("Snow", "Snow Algae", "Other"), 
                     values = c("snow" = 16,
                                "snow algae"=17,"other"=3))+
  theme_minimal(base_size = 20)

p1 


mean(my.plot.df["PC1"][my.plot.df["class"] == "snow"])
mean(my.plot.df["PC2"][my.plot.df["class"] == "snow"])

##### 20210730 #####
# plot PC1 vs PC2
my.plot.df.730 <- data.frame(
  class = as.factor(PCA.730$class),
  #conc = PCA.702$cell.concentration..cells.mL.,
  PC1 = PCA.730$b1_PCA_202,
  PC2 = PCA.730$b2_PCA_202,
  PC3 = PCA.730$b3_PCA_202
)

levels(my.plot.df.730$class)
levels(my.plot.df.730$class) <- c(levels(my.plot.df.730$class), "other")

my.plot.df.730["class"][my.plot.df.730["class"] == "rock"] <- "other"
my.plot.df.730["class"][my.plot.df.730["class"] == "water"] <- "other"
my.plot.df.730["class"][my.plot.df.730["class"] == "vegetation"] <- "other"

summary(my.plot.df.730)

# make ggplot  
p1.1 <- ggplot(my.plot.df.730, aes(x = PC1, y = PC3, color=class, shape=class)) +
  geom_point(size = 2) +
  #geom_abline(aes(intercept = -1000, slope =.02), linetype="dashed") +
  geom_hline(yintercept = -500, linetype = "dashed") +
  geom_vline(xintercept = 50000, linetype="dashed") +
  scale_color_manual(name="Class", labels = c("Snow", "Snow Algae", "Other"), 
                     values = c("snow" = "darkturquoise",
                                "snow algae"="darksalmon","other"="black")) +
  scale_shape_manual(name="Class", labels = c("Snow", "Snow Algae", "Other"), 
                     values = c("snow" = 16,
                                "snow algae"=17,"other"=3))+
  theme_minimal(base_size = 20)
p1.1 

ggarrange(p1, p1.1,
          labels = c("a", "b"),
          font.label=list(color="black",size=18),
          common.legend = TRUE,
          legend = "bottom",
          ncol = 2, nrow = 1)

ggsave(filename = "fig_3.png", device='png', dpi=500,width = 11, height=6, units = "in")

#### MicaSense Analysis ####
setwd("~/OneDrive - Western Washington University/Thesis/Manuscripts/Bagley UAV Paper/micasense")
mica.702 <- read.csv("micasense_20210702_accuracy.csv")
mica.730 <- read.csv("micasense_20210730_accuracy.csv")

#### 20210702 ####

mica.730 <- mica.730 %>%
  mutate(class = recode(class, rock = 'other', vegetation = 'other',
                        water = 'other'))
mica.702.df <- data.frame(
  class = as.factor(mica.702$class),
  date = rep(c("2021-07-02"), times=107),
  b1 = mica.702$b1,
  b2 = mica.702$b2,
  b3 = mica.702$b3,
  b4 = mica.702$b4,
  b5 = mica.702$b5,
  b6 = mica.702$b6,  
  b7 = mica.702$b7,
  b8 = mica.702$b8,
  b9 = mica.702$b9,
  b10 = mica.702$b10
)
mica.730.df <- data.frame(
  class = as.factor(mica.730$class),
  date = rep(c("2021-07-30"), times=106),
  b1 = mica.730$b1_micasen,
  b2 = mica.730$b2_micasen,
  b3 = mica.730$b3_micasen,
  b4 = mica.730$b4_micasen,
  b5 = mica.730$b5_micasen,
  b6 = mica.730$b6_micasen,  
  b7 = mica.730$b7_micasen,
  b8 = mica.730$b8_micasen,
  b9 = mica.730$b9_micasen,
  b10 = mica.730$b10_micase
)

mica.df <- rbind(mica.702.df, mica.730.df)
str(mica.df)

summary(my.plot.df)
mica.702.df$b5 <- mica.702.df$b5/32768
mica.702.df$b3 <- mica.702.df$b3/32768

p702<-ggplot(mica.702.df, aes(x = b3, y = b5, color=class, shape=class)) +
  geom_point(size = 2) +
  geom_abline(aes(intercept = .015, slope =1.02), linetype="dashed") +
  xlab("Band 3 Reflectance") +
  ylab("Band 5 Reflectance") +
  ylim(0,1) +
  xlim(0,1)+
  scale_color_manual(name="Class", labels = c("Snow", "Snow Algae", "Other"), 
                     values = c("snow" = "darkturquoise",
                                "snow algae"="darksalmon","other"="black")) +
  scale_shape_manual(name="Class", labels = c("Snow", "Snow Algae", "Other"), 
                     values = c("snow" = 16,
                                "snow algae"=17,"other"=3))+
  
  theme_minimal(base_size = 20)
p702

ggsave(filename = "b5b3_20210702.png", device='png', dpi=700)



#### 20210730 ####

my.plot.df.730 <- data.frame(
  class = as.factor(mica.730$class),
  #conc = mica.702$cell.concentration..cells.mL.,
  b1 = mica.730$b1_micasen/32768,
  b2 = mica.730$b2_micasen/32768,
  b3 = mica.730$b3_micasen/32768,
  b4 = mica.730$b4_micasen/32768,
  b5 = mica.730$b5_micasen/32768,
  b6 = mica.730$b6_micasen/32768,  
  b7 = mica.730$b7_micasen/32768,
  b8 = mica.730$b8_micasen/32768,
  b9 = mica.730$b9_micasen/32768,
  b10 = mica.730$b10_micase/32768
)
levels(my.plot.df.730$class) <- c(levels(my.plot.df.730$class), "other")

my.plot.df.730["class"][my.plot.df.730["class"] == "rock"] <- "other"
my.plot.df.730["class"][my.plot.df.730["class"] == "water"] <- "other"
my.plot.df.730["class"][my.plot.df.730["class"] == "vegetation"] <- "other"

mica.730.df$b5 <- mica.730.df$b5/32768
mica.730.df$b3 <- mica.730.df$b3/32768
summary(mica.730.df)

p730<-ggplot(mica.730.df, aes(x = b3, y = b5, color=class, shape=class)) +
  geom_point(size = 2) +
  geom_abline(aes(intercept = .015, slope =1.02), linetype="dashed") +
  xlab("Band 3 Reflectance") +
  ylab("Band 5 Reflectance") +
  ylim(0,1)+
  xlim(0,1)+
  scale_color_manual(name="Class", labels = c("Snow", "Snow Algae", "Other"), 
                     values = c("snow" = "darkturquoise",
                                "snow algae"="darksalmon","other"="black")) +
  scale_shape_manual(name="Class", labels = c("Snow", "Snow Algae", "Other"), 
                     values = c("snow" = 16,
                                "snow algae"=17,"other"=3))+
  theme_minimal(base_size = 20)
p730

ggsave(filename = "b5b3_20210730.png", device='png', dpi=700)

ggarrange(p702,p730,
          labels = c("a", "b"),
          font.label=list(color="black",size=18),
          common.legend = TRUE,
          legend="bottom",
          ncol = 2, nrow = 1)
ggsave(filename = "fig_4.png", device='png', dpi=500,width = 11, height=6, units = "in")
#### SNICAR reflectance plot ####
setwd("~/OneDrive - Western Washington University/Thesis/Manuscripts/Bagley UAV Paper/micasense")

##### 20210702 #####
mica_702 <- read.csv("micasense_20210702_accuracy.csv")

library(tidyr)
data_long_702 <- gather(mica_702, band, reflectance, b1:b10, factor_key=TRUE)
head(data_long_702)

#data_long_702 <- data_long_702[, 3:8]

library(dplyr)
data_long_702 <- data_long_702 %>%
  mutate(class = recode(class, rock = 'other', vegetation = 'other',
                        water = 'other'))

data_long_702 <- data_long_702 %>%
  mutate(band = recode(band, b1_micasense_bagley_20210702_WGS84_ortho = '1',
                       b2_micasense_bagley_20210702_WGS84_ortho= '2',
                       b3_micasense_bagley_20210702_WGS84_ortho = '3',
                       b4_micasense_bagley_20210702_WGS84_ortho = '4',
                       b5_micasense_bagley_20210702_WGS84_ortho = '5',
                       b6_micasense_bagley_20210702_WGS84_ortho = '6',
                       b7_micasense_bagley_20210702_WGS84_ortho = '7',
                       b8_micasense_bagley_20210702_WGS84_ortho = '8',
                       b9_micasense_bagley_20210702_WGS84_ortho = '9',
                       b10_micasense_bagley_20210702_WGS84_ortho = '10',))

data_long_702$wavelength <- data_long_702$band
data_long_702 <- data_long_702 %>%
  mutate(wavelength = recode(wavelength, 'b1' = '444', 'b2' = '475',
                             'b3' = '531','b4' = '560',
                             'b5' = '650','b6' = '668',
                             'b7' = '705','b8' = '717',
                             'b9' = '740','b10' = '842',))
data_long_702$wavelength <- as.numeric(as.character(data_long_702$wavelength))
data_long_702$wavelength.factor <- as.factor(data_long_702$wavelength)

summary(data_long_702)

data_long_702$reflectance <- data_long_702$reflectance/32768

p.702 <- ggplot() +
  geom_point(position=position_dodge(width=15)) +
  geom_smooth(data=data_long_702, method="loess",aes(x=wavelength, y=reflectance, color=class, group=class)) +
  #geom_smooth(data=snicar, method="loess", aes(x=wavelength.nm, y=albedo, fill="SNICAR output"), color="darkred")+
  scale_color_manual(name="Class", labels = c("Snow", "Snow Algae", "Other"), 
                     values = c("snow" = "darkturquoise",
                                "snow algae"="darksalmon","other"="azure4")) +
  ylim(0,1)+
  xlab("Wavelength (nm)")+ ylab("Reflectance") +
  theme_minimal(base_size = 18)
p.702

#get summary stats by group
tapply(data_long_702$reflectance, data_long_702$wavelength, summary)


data_long_702_low <- subset(data_long_702, data_long_702$wavelength < 600) 
data_long_702_high <- subset(data_long_702, data_long_702$wavelength > 600) 
tapply(data_long_702_low$reflectance, data_long_702_low$class, sd)
tapply(data_long_702_high$reflectance, data_long_702_high$class, sd)


#### 20210730 ####
mica_730 <- read.csv("micasense_20210730_accuracy.csv")

library(tidyr)
data_long_730 <- gather(mica_730, band, reflectance, b1_micasen:b10_micase, factor_key=TRUE)
head(data_long_730)

library(dplyr)
data_long_730 <- data_long_730 %>%
  mutate(class = recode(class, rock = 'other', vegetation = 'other',
                        water = 'other'))

data_long_730 <- data_long_730 %>%
  mutate(band = recode(band, b1_micasen = '1', b2_micasen = '2',
                       b3_micasen = '3',b4_micasen = '4',
                       b5_micasen = '5',b6_micasen = '6',
                       b7_micasen = '7',b8_micasen = '8',
                       b9_micasen = '9',b10_micase = '10',))

data_long_730$wavelength <- data_long_730$band
data_long_730 <- data_long_730 %>%
  mutate(wavelength = recode(wavelength, '1' = '444', '2' = '475',
                             '3' = '531','4' = '560',
                             '5' = '650','6' = '668',
                             '7' = '705','8' = '717',
                             '9' = '740','10' = '842',))
data_long_730$wavelength <- as.numeric(as.character(data_long_730$wavelength))

summary(data_long_730)

data_long_730$reflectance <- data_long_730$reflectance/32768

p.730 <- ggplot(data_long_730, aes(x=wavelength, y=reflectance, group=class,color=class)) +
  #geom_point(position=position_dodge(width=15)) +
  geom_smooth() +
  #geom_boxplot()+
  ylim(0,1)+
  scale_color_manual(name="Class", labels = c("Snow", "Snow Algae", "Other"), 
                     values = c("snow" = "darkturquoise",
                                "snow algae"="darksalmon","other"="azure4")) +
  xlab("Wavelength (nm)")+ ylab("Reflectance")+
  theme_minimal(base_size = 18)
p.730

ggarrange(p.702, p.730,
          labels = c("a", "b"),
          common.legend = TRUE,
          legend = "bottom",
          ncol = 2, nrow = 1)

tapply(data_long_730$reflectance, data_long_730$class, summary)


data_long_730_low <- subset(data_long_730, data_long_730$wavelength < 600) 
data_long_730_high <- subset(data_long_730, data_long_730$wavelength > 600) 
tapply(data_long_730_low$reflectance, data_long_730_low$class, sd)
tapply(data_long_730_high$reflectance, data_long_730_high$class, sd)

#### cell count correlation ####
library(tidyr)
library(dplyr)
setwd("C:/Users/healys2/OneDrive - Western Washington University/Thesis/Manuscripts/Bagley UAV Paper/micasense")

index.702 <- read.csv("micasense_20210702_accuracy.csv")
index.702 <- index.702[,3:16]
#index.702 <- index.702[ which(index.702$class=='snow algae'), ]
index.702$cell <- index.702$cell.concentration..cells.mL.

index.702$b5b3 <- (index.702$b5-0.015)/(1.02*index.702$b3)
index.702$cell.log <- log(index.702$cell)

index.702 <- index.702 %>%
  mutate(class = recode(class, rock = 'other', vegetation = 'other',
                        water = 'other'))

ggplot(index.702, aes(x=b3, y = b5, color = class)) +
  geom_point() +
  theme_minimal()

ggplot(index.702, aes(x=b5b3, y = cell.log)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_minimal()

index.730 <- read.csv("micasense_20210730_accuracy.csv")
index.730 <- index.730 %>%
  mutate(class = recode(class, rock = 'other', vegetation = 'other',
                        water = 'other'))
index.730$cell <- index.730$Cell.concentration..cells.mL.

index.730$b5b3 <- (index.730$b5-0.015)/(1.02*index.730$b3)
index.730$cell.log <- log(index.730$cell)

index.702 <- index.702[ which(index.702$class=='snow algae'), ]
index.730 <- index.730[ which(index.730$class=='snow algae'), ]
index.702.df <- data.frame(class = index.702$class,
                           sample = index.702$CID,
                           cell = index.702$cell,
                           index = index.702$b5b3,
                           cell.ln = index.702$cell.log,
                           date = rep(c("2021-07-02"), times=9))
index.730.df <- data.frame(class = index.730$class,
                           sample = rep(c("bb_20210730_"), times=10),
                           cell = index.730$cell,
                           index = index.730$b5b3,
                           cell.ln = index.730$cell.log,
                           date = rep(c("2021-07-30"), times=10))
index <- rbind(index.702.df, index.730.df)

write.csv(index,"mica_index.csv", row.names = FALSE)
index <- read.csv("mica_index.csv")
index <- index[ which(index$class=='snow algae'), ]

ggplot(index, aes(x=index, y = cell.ln)) +
  geom_smooth(method="lm", se=FALSE,color="black", linetype="dashed") +
  geom_smooth(method="lm", aes(color=date)) +
  geom_point(aes(color=date, shape=date), size=2) +
  xlab("Optimized Red/Green Index") + ylab("ln(Cell Concentration (cells/mL))")+
  scale_color_manual(name="Date", labels=c("7/2/21"="2 July","7/30/21"="30 July"), values=c("7/2/21"="darksalmon","7/30/21"="darkturquoise"))+
  scale_shape_manual(name="Date", labels=c("7/2/21"="2 July","7/30/21"="30 July"), values=c("7/2/21"=16,"7/30/21"=17))+
  theme_minimal(base_size = 20)


ggsave(filename = "Figure_5.png", device='png', dpi=2100)

log.cell.reg<- lm(cell.log~index, index)
summary(log.cell.reg)

ggplot(index, aes(x=index, y = cell.log, color=date)) +
  geom_point() +
  geom_smooth(method="lm") +
  xlab("MicaSense Index") + ylab("ln(Cell Concentration (cells/mL))")+
  scale_color_discrete(name="Date", labels=c("July 2, 2021","July 30, 2021"))+
  theme_minimal()
ggsave(filename = "log_cell_conc_vs_index_date.png", device='png', dpi=700)

index.702.df <- index.702.df[ which(index.702.df$class=='snow algae'), ]
index.730.df <- index.730.df[ which(index.730.df$class=='snow algae'), ]
log.cell.reg.702<- lm(cell.log~index, index.702.df)
summary(log.cell.reg.702)
log.cell.reg.730<- lm(cell.log~index, index.730.df)
summary(log.cell.reg.730)

#pigments
index <- read.csv("mica_index.csv")
index$carotenoids <- index$ast +index$lut+index$bc
index$carotenoids.chla <- index$carotenoids/index$chla
summary(index$carotenoids.chla)

ggplot(index, aes(x=index, y=ast.chla)) +
  geom_point(aes(color=date))+
  geom_smooth(method = 'lm', color="darkred")+
  #guides(color=FALSE)+
  xlab("MicaSense Index")+
  ylab("Astaxanthin/Chlorophyll A") +
  scale_color_discrete(name="Date", labels=c("July 2, 2021", "July 30, 2021"))+
  theme_minimal()
ggsave(filename = "ast_chla_vs_index.png", device='png', dpi=700)

ast.reg<- lm(ast.chla~index, index)
summary(ast.reg)

setwd("C:/Users/healys2/OneDrive - Western Washington University/Thesis/Manuscripts/Bagley UAV Paper/micasense")
mica.index <- read.csv("algae_micasense_index.csv")

library(tidyr)
library(ggplot2)

mica.index.long <- gather(mica.index, index, value, algae_index_702_b3b5_b5:mica_0730_b3b5, factor_key=TRUE)
head(mica.index.long)

ggplot(mica.index.long, aes(x=value, y=cell.conc, color=index)) +
  geom_point()+
  geom_smooth()+
  theme_minimal()

ggplot(mica.index, aes(x=index_b3b5, y=cell.conc)) +
  geom_point()+
  geom_smooth(method = 'lm',aes(color="darksalmon"))+
  guides(color=FALSE)+
  xlab("Index Value")+
  ylab("Cell Concentration (cells/mL)") +
  theme_minimal()

ggsave(filename = "cell_conc_vs_index.png", device='png', dpi=700)

pigments.index <- read.csv("pigments_micasense_index.csv")

pigments.index.long <- gather(pigments.index, pigment, value, ast:bc.chla, factor_key=TRUE)
head(pigments.index.long)

ggplot(pigments.index.long, aes(x=index_b3b5, y=value, color=pigment)) +
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,aes(color=pigment))+
  #guides(color=FALSE)+
  xlab("Index Value")+
  ylab("Pigment") +
  theme_minimal()

ggplot(pigments.index, aes(x=index_b3b5, y=ast, color=Date)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  #guides(color=FALSE)+
  xlab("Index Value")+
  ylab("Pigment") +
  theme_minimal()

pigments.index$carotenoids <- pigments.index$ast +pigments.index$lut + pigments.index$bc

pigments.index$caro.chla <- pigments.index$carotenoids/pigments.index$chla

p.pig<-ggplot(pigments.index, aes(x=index_b3b5, y=carotenoids, color=Date)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  #guides(color=FALSE)+
  xlab("Index Value")+
  ylab(expression( "Carotenoid Concentration (" ~mu~ "g/mL)")) +
  scale_color_manual(name="Date", labels = c("July 2, 2021", "July 30, 2021"), 
                     values = c("7/2/2021" = "darkturquoise",
                                "7/30/2021"="darksalmon")) +
  theme_minimal()
head(mica.index$Date)

ggsave(filename = "carotenoid_conc_vs_index_date.png", device='png', dpi=700)

p.cell<-ggplot(mica.index, aes(x=index_b3b5, y=cell.conc, color= Date)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  #guides(color=FALSE)+
  xlab("Index Value")+
  ylab(expression( "Cell Concentration (cells/mL)")) +
  scale_color_manual(name="Date", labels = c("July 2, 2021", "July 30, 2021"), 
                     values = c("7/2/2021 0:00" = "darkturquoise",
                                "7/30/2021 0:00"="darksalmon")) +
  theme_minimal()

ggsave(filename = "cell_conc_vs_index_date.png", device='png', dpi=700)

library(ggpubr)
ggarrange(p.cell, p.pig,
          labels = c("A", "B"),
          common.legend = TRUE,
          legend = "bottom",
          ncol = 2, nrow = 1)

ggsave(filename = "cell_pig_conc_vs_index_date.png", device='png', dpi=700)

cell.702 <- mica.index[1:9,]
cell.730 <-mica.index[10:19,]
pig.730 <- pigments.index[10:19,]

p.cell730<-ggplot(cell.730, aes(x=index_b3b5, y=cell.conc)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab("Index Value")+
  ylab(expression( "Cell Concentration (cells/mL)")) +
  theme_minimal()

p.pig730<-ggplot(pig.730, aes(x=index_b3b5, y=carotenoids)) +
  geom_point()+
  geom_smooth(method = 'lm')+
  xlab("Index Value")+
  ylab(expression( "Carotenoid Concentration (" ~mu~ "g/mL)")) +
  theme_minimal()

ggarrange(p.cell730, p.pig730,
          labels = c("A", "B"),
          common.legend = TRUE,
          legend = "bottom",
          ncol = 2, nrow = 1)

ggsave(filename = "cell_pig_conc_vs_index_730.png", device='png', dpi=700)


#### 730 cell regression ####
cell730.reg<- lm(cell.conc~index_b3b5, cell.730)
print(summary(cell730.reg))

cell702.reg<- lm(cell.conc~index_b3b5, cell.702)
print(summary(cell702.reg))

cell.reg<- lm(cell.conc~index_b3b5, mica.index)
print(summary(cell.reg))


#### SNICAR output ####
snicar <- read.csv("SNICAR output.csv")

snicar$wavelength.nm <- snicar$wavelength*1000
snicar$date  <- ifelse(grepl("702", snicar$sample), "2021-07-2", "2021-07-30")
summary(snicar)

snicar$date <- as.factor(snicar$date)
snicar_low <- subset(snicar, snicar$wavelength.nm < 600) 
tapply(snicar_low$albedo, snicar_low$date, sd)

#install.packages("remotes")
#remotes::install_github("coolbutuseless/ggpattern")
library(ggplot2)
library(ggpattern)
snicar_sub <- subset(snicar, snicar$wavelength.nm > 350) 
snicar_sub <- subset(snicar_sub, snicar_sub$wavelength.nm < 1500)

ggplot(snicar_sub, aes(x=wavelength.nm, y=albedo, color=cell.concentration))+
  geom_rect(aes(xmin=430, xmax=458, ymin=0, ymax=Inf), fill="lightgrey", color="darkgrey",alpha=0.25) +
  geom_rect(aes(xmin=459, xmax=491, ymin=0, ymax=Inf), fill="lightgrey", color="darkgrey",alpha=0.25) +
  geom_rect(aes(xmin=524, xmax=538, ymin=0, ymax=Inf), fill="darkgrey", color="black") +
  geom_rect(aes(xmin=546, xmax=574, ymin=0, ymax=Inf), fill="lightgrey", color="darkgrey",alpha=0.25) +
  geom_rect(aes(xmin=642, xmax=658, ymin=0, ymax=Inf), fill="darkgrey", color="black") +
  geom_rect(aes(xmin=661, xmax=675, ymin=0, ymax=Inf), fill="lightgrey", color="darkgrey",alpha=0.25) +
  geom_rect(aes(xmin=700, xmax=710, ymin=0, ymax=Inf), fill="lightgrey", color="darkgrey",alpha=0.25) +
  geom_rect(aes(xmin=711, xmax=723, ymin=0, ymax=Inf), fill="lightgrey", color="darkgrey",alpha=0.25) +
  geom_rect(aes(xmin=731, xmax=749, ymin=0, ymax=Inf), fill="lightgrey", color="darkgrey",alpha=0.25) +
  geom_rect(aes(xmin=813, xmax=871, ymin=0, ymax=Inf), fill="lightgrey", color="darkgrey",alpha=0.25) +
  geom_segment(aes(x=531, y=-0.07, xend=531, yend=-0.008),color="black", size=.9, lineend = "round", linejoin="round", arrow=arrow(length = unit(0.15, "inches"))) +
  geom_segment(aes(x=650, y=-0.07, xend=650, yend=-0.008),color="black", size=.9, lineend = "round", linejoin="round", arrow=arrow(length = unit(0.15, "inches"))) +
  geom_line(aes(group=sample, linetype=date), size=1)+
  coord_cartesian(ylim = c(0, 1), clip="off") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(350, 1500), expand = c(0,0), oob = scales::squish) +
  scale_color_continuous(high = "darkred", low = "lightpink", name = "Cell concentration\n(cells/mL)",
                         breaks = c(100000,500000,1000000), labels = c("100,000", "500,000","1,000,000"))+
  scale_linetype(name="Survey Date", labels = c("2 July", "30 July"))+
  xlab("Wavelength (nm)")+
  ylab("Hemispheric Albedo")+
  theme_minimal(base_size = 20)


ggsave(filename = "snicar_out_micabands_2.png", device='png', dpi=2100)
