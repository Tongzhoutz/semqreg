### lll
 #load(file="~/Dropbox/Research_Proposal/data_center/rda/FAM1999ER.rda",verbose=T)
# FAM98 <- as.data.table(x)
#load(file="~/Dropbox/Research_Proposal/data_center/rda/FAM2007ER.rda")


#load(file="~/Dropbox/Research_Proposal/data_center/rda/FAM2001ER.rda")
#FAM100 <- as.data.table(x)

#load(file="~/Dropbox/Research_Proposal/data_center/rda/FAM2003ER.rda")
#FAM102 <- as.data.table(x)

#load(file="~/Dropbox/Research_Proposal/data_center/rda/FAM2005ER.rda")
#FAM104 <- as.data.table(x)
#FAM106 <- as.data.table(x)
#load(file="~/Dropbox/Research_Proposal/data_center/rda/FAM2009ER.rda")
#FAM108 <- as.data.table(x)

#rm(x)
#fam2005 <- as.data.table(import("~/Dropbox/Research_Proposal/data_center/rda/FAM2005ER.rda"))
 data_path <- "~/Dropbox/Research_Proposal/data_center/rda/"
 files <- c("FAM1999ER.rda","FAM2001ER.rda","FAM2003ER.rda","FAM2005ER.rda","FAM2007ER.rda","FAM2009ER.rda","FAM2011ER.rda","FAM2013ER.rda","FAM2015ER.rda","FAM2017ER.rda","FAM2019ER.rda")

fam <- files %>% 
  map( ~as.data.table(import_list( file.path(data_path,.)  )) )


 
dt = data.table::data.table(fam)
rm(fam)
names(dt[[1]]) <- c("fam98","fam100","fam102","fam104","fam106","fam108","fam110",
                    "fam112","fam114","fam116","fam118")

DT1 <- data.table( x = rnorm(100), y = rnorm(100))
DT1[, x %*% t(x)]
DT2 <- data.table(z = rt(100,3))

DT3 <- data.table::data.table( a = list(DT1,DT2))
DT3[, DT1$x %*% t(DT2$z)]
