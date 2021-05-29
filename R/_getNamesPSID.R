fsize <- getNamesPSID("ER13009",cwf, years= all_year)
id <- getNamesPSID("ER17002",cwf,years = all_year) 
age <- getNamesPSID("ER13010",cwf,years = all_year)
sex <- getNamesPSID("ER13011",cwf,years = all_year)
agew <- getNamesPSID("ER13012",cwf,years = all_year)
kids <- getNamesPSID("ER13013",cwf,years = all_year)
marit <- getNamesPSID("ER13021",cwf,years = all_year)
marit_generated <- getNamesPSID("ER16423",cwf,years = all_year)
race1 <- getNamesPSID("ER15928",cwf,years = all_year)
race2 <- getNamesPSID("ER15929",cwf,years = all_year)
race3 <- getNamesPSID("ER15930",cwf,years = all_year)
race4 <- getNamesPSID("ER15931",cwf,years = all_year)
wrace1 <- getNamesPSID("ER15836",cwf,years = all_year)
wrace2 <- getNamesPSID("ER15837",cwf,years = all_year)
wrace3 <- getNamesPSID("ER15838",cwf,years = all_year)
wrace4 <- getNamesPSID("ER15839",cwf,years = all_year)
outkid <- getNamesPSID("ER14976",cwf,years = all_year)
smsa <- getNamesPSID("ER16432",cwf,years = all_year)
house <- getNamesPSID("ER13041",cwf,years = all_year)
weight <- getNamesPSID("ER16518",cwf,years = all_year)
weight_CS <- getNamesPSID("ER16519",cwf,years = c(1999,2001,2003,2017,2019))
## home insurance 
homeinsure <- getNamesPSID("ER13043",cwf,years = all_year)
## Asset Stock
cash <- getNamesPSID("ER43586",cwf,years = all_year)
cash[11,2] = "ER73848"
bond <- getNamesPSID("ER43607",cwf,years = all_year)
stocks <- getNamesPSID("ER43558",cwf, years = all_year)
real_estate <- getNamesPSID("ER43544",cwf,years = all_year)
real_estate[8:11,2] = getNamesPSID("ER61723",cwf,years = c(2013,2015,2017,2019))[,2]

carval <- getNamesPSID("ER43548",cwf,years = all_year)
busval <- getNamesPSID("ER43553",cwf,years = all_year)
penval <- getNamesPSID("ER43580",cwf,years = all_year)
mortgage1 <- getNamesPSID("ER42043",cwf,years = all_year)
mortgage2 <- getNamesPSID("ER42062",cwf,years = all_year)
other_debt <- getNamesPSID("ER43612",cwf,years = all_year)
#### health expenditures
hinsurance <- getNamesPSID("ER15780",cwf,years = all_year)
hinsurance[8:11,2] = getNamesPSID("ER70683",cwf,years = c(2013,2015,2017,2019))[,2]

nurse <- getNamesPSID("ER15781",cwf,years = all_year)
###getNamesPSID("ER70689",cwf,years = all_year)   times 2? 

doctor <- getNamesPSID("ER15787",cwf,years = all_year)
prescription <- getNamesPSID("ER15793",cwf,years = all_year)
totalhealth <- getNamesPSID("ER15799",cwf,years = all_year)

#### Utilities
electric <- getNamesPSID("ER42114",cwf,years = all_year)
heating <- getNamesPSID("ER42112",cwf,years = all_year)
water <- getNamesPSID("ER42118",cwf,years = all_year)
miscutils_t <- getNamesPSID("ER42124",cwf,years = all_year)
telephone <- getNamesPSID("ER42120",cwf,years = c(2005,2007,2009,2011,2013,2015,2017,2019))
### Car expenditure 
carinsList <- getNamesPSID("ER13191",cwf, years=all_year)
carrepairList <- getNamesPSID("ER13195",cwf, years=all_year)

gasolineList <- getNamesPSID("ER13196",cwf, years=all_year)
parkingList <- getNamesPSID("ER13197",cwf, years=all_year)
busfareList <- getNamesPSID("ER13198",cwf, years=all_year)
taxifareList <- getNamesPSID("ER13199",cwf, years=all_year)
othertransList <- getNamesPSID("ER13200",cwf, years=all_year)
### Education expenditures
tuitionList <- getNamesPSID("ER13202",cwf, years=all_year)
otherschoolList <- getNamesPSID("ER13204",cwf, years=all_year)
### Child Care
childcareList <- getNamesPSID("ER14232",cwf, years=all_year)

### state  (only appear after 2001)
stateList <- getNamesPSID("ER17004",cwf,years = all_year)
### fchg   (only appear after 2001)
fchgList <- getNamesPSID("ER17007",cwf,years = all_year)

###  empst for head
empst1List <- getNamesPSID("ER42140",cwf,years = all_year)
empst2List <- getNamesPSID("ER42141",cwf,years = all_year)
empst3List <- getNamesPSID("ER42142",cwf,years = all_year)
wempst1List <- getNamesPSID("ER42392",cwf,years = all_year)
wempst2List <- getNamesPSID("ER42393",cwf,years = all_year)
wempst3List <- getNamesPSID("ER42394",cwf,years = all_year)

### education
eduList <- getNamesPSID("ER46981",cwf,years = all_year)
weduList <- getNamesPSID("ER46982",cwf,years = all_year)

## food variables
fs_jan <- getNamesPSID("ER14258",cwf,years = all_year)
fs_feb <- getNamesPSID("ER14259",cwf,years = all_year)
fs_mar <- getNamesPSID("ER14260",cwf,years = all_year)
fs_apr <- getNamesPSID("ER14261",cwf,years = all_year)
fs_may <- getNamesPSID("ER14262",cwf,years = all_year)
fs_jun <- getNamesPSID("ER14263",cwf,years = all_year)
fs_jul <- getNamesPSID("ER14264",cwf,years = all_year)
fs_aug <- getNamesPSID("ER14265",cwf,years = all_year)
fs_sep <- getNamesPSID("ER14266",cwf,years = all_year)
fs_oct <- getNamesPSID("ER14267",cwf,years = all_year)
fs_nov <- getNamesPSID("ER14268",cwf,years = all_year)
fs_dec <- getNamesPSID("ER14269",cwf, years = all_year)

## Asset income variables
getNamesPSID("ER14447",cwf,years = all_year)
