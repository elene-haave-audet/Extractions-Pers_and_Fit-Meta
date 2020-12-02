# Cook 2011 Rec 62
# Phenotypic extractions from figures

# Load libraries====
library(here)

# Set wd====
dir<-here()

# Extractions====
c62.temp<-read.csv("Data/Cook2011Temp.csv")
cor.test(c62.temp$TempScore, c62.temp$Preg) #r=-0.1064453

c62.chute<-read.csv("Data/Cook2011Chute.csv")
cor.test(c62.chute$ChuteScore, c62.chute$Preg) #r=-0.06475307

c62.exit<-read.csv("Data/Cook2011Exit.csv")
cor.test(c62.exit$ExitScore, c62.exit$Preg) #r=-0.1063492