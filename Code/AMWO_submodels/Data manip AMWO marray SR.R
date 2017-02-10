#------------------------------------------------------------------------------
# Data manipulation of banding data to create m-arrays for band recovery models
#------------------------------------------------------------------------------
rm(list=ls())
# Read in dataset
#raw.data<-read.csv(file.choose())
raw<-read.csv("~/Google Drive/AMWO IPM/Datasets/M array table AMWO CSV.csv")

########################################################################
#cleaning data
########################################################################

#bring in B.month, convert to season
clean<-matrix(NA,nrow=length(raw$B.Month),ncol=1)
clean<-data.frame(clean)
clean[raw$B.Month<=6,]<-1
clean[raw$B.Month>6,]<-2

#Bring in B.year
clean[,2]<-raw$B.Year

#bring in recovery year and account for recoveries occuring in Jan-March
clean[,3]<-NA
clean[raw$R.Month>=4,3]<-raw[raw$R.Month>=4,"R.Year"]
#adjust if you don't want to include birds recovered in March
clean[raw$R.Month<4,3]<-raw[raw$R.Month<4,"R.Year"]-1

#bring in B region
clean[,4]<-0
clean[raw$B.Flyway==1,4]<-1
clean[raw$B.Flyway%in%2:3,4]<-2
clean[raw$B.Flyway==6&raw$BRegion..STA%in%c("QC","NS","NB","PE","NF","PQ"),4]<-1
clean[raw$B.Flyway==6&raw$BRegion..STA%in%c("ONT"),4]<-2

#bring in R region, this is only to exlude region crossers
clean[,5]<-0 #specify different number from previous step to flag it in the next step
clean[raw$B.Flyway==1,5]<-1
clean[raw$B.Flyway%in%2:3,5]<-2
clean[raw$B.Flyway==6&raw$BRegion..STA%in%c("QC","NS","NB","PE","NF","PQ"),5]<-1
clean[raw$B.Flyway==6&raw$BRegion..STA%in%c("ONT"),5]<-2

# pull out places you don't care about
raw<-raw[clean$V4!=0|clean$V5!=0,]
clean<-clean[clean$V4!=0|clean$V5!=0,]
clean<-clean[,1:4] #remove R.state becuase it is redundant 

#bring in age
# local = 1, hatch year = 2, adult = 3
clean[,5]<-NA
clean[raw$Age..VAGE=="After Hatch Year",5]<-3
clean[raw$Age..VAGE=="After Second Year",5]<-3
clean[raw$Age..VAGE=="After Third Year",5]<-3
clean[raw$Age..VAGE=="Second Year",5]<-3
clean[raw$Age..VAGE=="Unknown",5]<-NA
clean[raw$Age..VAGE=="Hatch Year",5]<-2
clean[raw$Age..VAGE=="Local",5]<-1
#remove unknowns
raw<-raw[!is.na(clean[,5]),]
clean<-clean[!is.na(clean[,5]),]

# get rid of hact years in months 5 and 6
clean <- clean[!(raw$Age..VAGE=="Hatch Year"&raw$B.Month%in%c(5:6)),]
raw <- raw[!(raw$Age..VAGE=="Hatch Year"&raw$B.Month%in%c(5:6)),]

#bring in sex and convert to age class
# 1=local, 2=juv, 3=male, 4=female
clean[,6]<-NA
clean[clean[,5]%in%1:2,6]<-clean[clean[,5]%in%1:2,5]
clean[raw$Sex..VSEX%in%c("Male","Male; from subsequent encounter")&clean[,5]==3,6]<-3
clean[raw$Sex..VSEX%in%c("Female","Female; from subsequent encounter")&clean[,5]==3,6]<-4
#remove unknown adults
raw<-raw[!(is.na(clean[,6])&clean[,5]==3),]
clean<-clean[!(is.na(clean[,6])&clean[,5]==3),]

# remove unwanted banding ????? I'm guessing this is a thing
clean<-clean[!raw$B.Month%in%c(10:12,1:3),]
raw<-raw[!raw$B.Month%in%c(10:12,1:3),]

clean[,7]<-1
colnames(clean)<-c("bSeason","bYear","rYear","region","age","class","dummy")
########################################################################
#create the marray
########################################################################
Year<-unique(clean$bYear)
NYear<-length(Year)
Season<-unique(clean$bSeason)
NSeason<-length(Season)
Class<-unique(clean$class)
NClass<-length(Class)
Region<-unique(clean$region)
NRegion<-length(Region)

awc<-array(NA,dim=c(NYear,NYear,NSeason,NClass,NRegion),
           dimnames =list(Year, Year, c("spring","not_spring"),
                       c("local","Hatch_Year","Adult_Male","Adult_Female"),
                       c("Eastern","Central")))
for (s in 1:NSeason){
  for (cc in 1:NClass){
    for (i in 1:NRegion){
      for (b in 1:NYear){
        for (r in 1:NYear){
                awc[b,r,s,cc,i]<-sum(clean[clean$bYear==Year[b]&clean$rYear==Year[r]&clean$class==Class[cc]&clean$region==Region[i],7])
        }}}}}


