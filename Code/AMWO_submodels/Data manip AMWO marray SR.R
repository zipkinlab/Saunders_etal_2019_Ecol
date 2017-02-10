#------------------------------------------------------------------------------
# Data manipulation of banding data to create m-arrays for band recovery models
#------------------------------------------------------------------------------
rm(list=ls())
# Read in dataset
#raw.data<-read.csv(file.choose())
raw<-read.csv("AMWO recoveries.csv")  #reading in CSV from GitHub-linked timberdoodle folder

########################################################################
#cleaning data
########################################################################
#SS added following conditions to subset raw data:

#only use status 3 birds
raw<-subset(raw,Status==3)

#only use how obtained category 1 (shot)
raw<-subset(raw,How.Obt==1)

#only use B.Year from 1963 onwards
raw<-subset(raw,B.Year>=1963)

#subset B.Month between 4 and 9 to cover our 2 seasons
raw<-subset(raw,B.Month>=4&B.Month<=9)

#take out recoveries in March
raw<-subset(raw,R.Month!=3)

#bring in B.month, convert to season
clean<-matrix(NA,nrow=length(raw$B.Month),ncol=1)
clean<-data.frame(clean)
clean[raw$B.Month>=4&raw$B.Month<=6,]<-1  
clean[raw$B.Month>=7&raw$B.Month<=9,]<-2

#Bring in B.year
clean[,2]<-raw$B.Year  

#bring in recovery year and account for recoveries occurring in Jan-Feb
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
clean[raw$R.Flyway==1,5]<-1
clean[raw$R.Flyway%in%2:3,5]<-2
clean[raw$R.Flyway==6&raw$RRegion..STA%in%c("QC","NS","NB","PE","NF","PQ"),5]<-1
clean[raw$R.Flyway==6&raw$RRegion..STA%in%c("ONT"),5]<-2

# pull out places you don't care about
#raw<-raw[clean$V4!=0|clean$V5!=0,]
#test<-clean[clean$V4==0|clean$V5==0,]
# removes lines of region crossers
raw<-raw[clean$V4==clean$V5,]
clean<-clean[clean$V4==clean$V5,]
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

# get rid of hatch years in months 4, 5 and 6
clean <- clean[!(raw$Age..VAGE=="Hatch Year"&raw$B.Month%in%c(4:6)),]
raw <- raw[!(raw$Age..VAGE=="Hatch Year"&raw$B.Month%in%c(4:6)),]

#bring in sex and convert to age class so this is more like a sex-ageclass column
# 1=local, 2=juv, 3=male, 4=female
clean[,6]<-NA
clean[clean[,5]%in%1:2,6]<-clean[clean[,5]%in%1:2,5]     
clean[raw$Sex..VSEX%in%c("Male","Male; from subsequent encounter")&clean[,5]==3,6]<-3
clean[raw$Sex..VSEX%in%c("Female","Female; from subsequent encounter")&clean[,5]==3,6]<-4

#remove unknown adults for now--only losing 20 individuals if we don't include them
raw<-raw[!(is.na(clean[,6])&clean[,5]==3),]
clean<-clean[!(is.na(clean[,6])&clean[,5]==3),]

#adding dummy column to use sum function below for marray
clean[,7]<-1
colnames(clean)<-c("bSeason","bYear","rYear","region","age","class","dummy")

########################################################################
#create the marray
########################################################################
Year<-unique(clean$bYear)
Year<-sort(Year)          #SS sorted
NYear<-length(Year)
Season<-unique(clean$bSeason)
NSeason<-length(Season)
Class<-unique(clean$class)
Class<-sort(Class)        #SS sorted
NClass<-length(Class)
Region<-unique(clean$region)
Region<-sort(Region)       #SS sorted
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

#take a look at subset of giant marray--basically this is 16 marrays
awc[1:20,1:20,1,3,1]

save(awc, file="AMWO_Marray.rda")

#Creating separate M-Array that is more compatible with Todd's example
#Will need to be finished at a later date
awc.todd <- array(NA, dim=c(2*NYear, NYear, NClass, NRegion),
                  dimnames =list(Year, Year,
                                 c("local","Hatch_Year","Adult_Male","Adult_Female"),
                                 c("Eastern","Central")))

  for (cc in 1:NClass){
    for (i in 1:NRegion){
      for (b in 1:NYear){
        for (r in 1:NYear){
          awc.todd[b,r,cc,i]<-sum(clean[clean$bYear==Year[b]&clean$rYear==Year[r]&clean$class==Class[cc]&clean$region==Region[i],7])
        }}}}

save(awc.todd, file="AMWO_Todd_Marray.rda")

#---------------------------------------------------------------------------
#need to add last column of unrecovered individuals to marray
#THIS IS NOT FINISHED YET EITHER
#---------------------------------------------------------------------------
#bring in bandings file
bands<-read.csv("AMWO bandings.csv")

#need to summarize bandings according to: banding year, region, class, season

#only use status 3 birds
bands<-subset(bands,Status==3)

#only use B.Year from 1963 onwards
bands<-subset(bands,B.Year>=1963)

#bring in B.month, convert to season
clean<-matrix(NA,nrow=length(bands$B.Month),ncol=1)
clean<-data.frame(clean)
clean[bands$B.Month<=6,]<-1  #shouldn't we change this to between 4 and 6?
clean[bands$B.Month>6,]<-2   #shouldn't we change this to between 7 and 9?

#Bring in B.year
clean[,2]<-bands$B.Year  

#bring in B region
clean[,3]<-0
clean[bands$B.Flyway==1,3]<-1
clean[bands$B.Flyway%in%2:3,3]<-2
clean[bands$B.Flyway==6&bands$BRegion..STA%in%c("QC","NS","NB","PE","NF","PQ"),3]<-1
clean[bands$B.Flyway==6&bands$BRegion..STA%in%c("ONT"),3]<-2

# pull out places you don't care about
bands<-bands[clean$V3!=0,]
clean<-clean[clean$V3!=0,]
clean<-clean[,1:3] #remove R.state becuase it is redundant 

#bring in age
# local = 1, hatch year = 2, adult = 3
clean[,4]<-NA
clean[bands$Age..VAGE=="After Hatch Year",4]<-3
clean[bands$Age..VAGE=="After Second Year",4]<-3
clean[bands$Age..VAGE=="After Third Year",4]<-3
clean[bands$Age..VAGE=="Second Year",4]<-3
clean[bands$Age..VAGE=="Unknown",4]<-NA
clean[bands$Age..VAGE=="Hatch Year",4]<-2
clean[bands$Age..VAGE=="Local",4]<-1
#remove unknowns
bands<-bands[!is.na(clean[,4]),]
clean<-clean[!is.na(clean[,4]),]

# get rid of hatch years in months 5 and 6
clean <- clean[!(bands$Age..VAGE=="Hatch Year"&bands$B.Month%in%c(5:6)),]
bands <- bands[!(bands$Age..VAGE=="Hatch Year"&bands$B.Month%in%c(5:6)),]

#bring in sex and convert to age class
# 1=local, 2=juv, 3=male, 4=female




#remove unknown adults for now? Can treat unknowns via mixtures acc. to Todd


# remove unwanted banding periods (Oct-Dec) for now


