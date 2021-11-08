
all.idat.files <- list.files(idat.pooled.files.dir,pattern = ".idat")

all.idat.files.clean <- (gsub("_Grn.idat|_Red.idat","",all.idat.files))
all.idat.files.clean <- all.idat.files.clean[!duplicated(all.idat.files.clean)]

all.idat.files.clean %>% length() ##30

tumor.idat.files  <- all.idat.files.clean[all.idat.files.clean %in% tumor$Sample_Name]

normal.idat.files <- all.idat.files.clean[all.idat.files.clean %in% normal$Sample_Name]

weird.idats <- all.idat.files.clean[!all.idat.files.clean %in% tumor.idat.files & !all.idat.files.clean %in% normal.idat.files]

weird.idats %>% length() ##10

"6164655053_R06C01" 
"6164655053_R06C02" 
"6264496109_R01C01"
"6285650057_R02C02"
"7796806055_R01C02"
"7796806075_R06C01"
"7796806095_R02C01"
"7796806101_R02C02"
"7796806101_R05C01"
"8795194084_R06C02"
