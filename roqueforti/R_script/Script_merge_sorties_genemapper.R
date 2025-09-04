setwd("/Users/alodie/Desktop/Manuela_recombinaison/2eme essai_YingChu")

MP01 <- read.table("Proq12.txt", header=TRUE)
MP01$Sample_File_new <- gsub("MP01_", "", MP01$Sample_File)
MP01$Sample_Name_new <- gsub("MP01_", "", MP01$Sample_Name)

Ch01f03b <- data.frame(subset(MP01, Marker=="Ch01f03b")$Sample_File_new, subset(MP01, Marker=="Ch01f03b")$Sample_Name_new, subset(MP01, Marker=="Ch01f03b")$Marker, subset(MP01, Marker=="Ch01f03b")$Allele_1, subset(MP01, Marker=="Ch01f03b")$Allele_2)
names(Ch01f03b)[names(Ch01f03b)=="subset.MP01..Marker.....Ch01f03b...Sample_File_new"] <- "Sample_File"
names(Ch01f03b)[names(Ch01f03b)=="subset.MP01..Marker.....Ch01f03b...Sample_Name_new"] <- "Sample_Name"
names(Ch01f03b)[names(Ch01f03b)=="subset.MP01..Marker.....Ch01f03b...Allele_1"] <- "Ch01f03b_A"
names(Ch01f03b)[names(Ch01f03b)=="subset.MP01..Marker.....Ch01f03b...Allele_2"] <- "Ch01f03b_B"
Ch01f03b$subset.MP01..Marker.....Ch01f03b...Marker <- NULL

Ch01h10 <- data.frame(subset(MP01, Marker=="Ch01h10")$Sample_File_new, subset(MP01, Marker=="Ch01h10")$Sample_Name_new, subset(MP01, Marker=="Ch01h10")$Marker, subset(MP01, Marker=="Ch01h10")$Allele_1, subset(MP01, Marker=="Ch01h10")$Allele_2)
names(Ch01h10)[names(Ch01h10)=="subset.MP01..Marker.....Ch01h10...Sample_File_new"] <- "Sample_File"
names(Ch01h10)[names(Ch01h10)=="subset.MP01..Marker.....Ch01h10...Sample_Name_new"] <- "Sample_Name"
names(Ch01h10)[names(Ch01h10)=="subset.MP01..Marker.....Ch01h10...Allele_1"] <- "Ch01h10_A"
names(Ch01h10)[names(Ch01h10)=="subset.MP01..Marker.....Ch01h10...Allele_2"] <- "Ch01h10_B"
Ch01h10$subset.MP01..Marker.....Ch01h10...Marker <- NULL

Ch01h01 <- data.frame(subset(MP01, Marker=="Ch01h01")$Sample_File_new, subset(MP01, Marker=="Ch01h01")$Sample_Name_new, subset(MP01, Marker=="Ch01h01")$Marker, subset(MP01, Marker=="Ch01h01")$Allele_1, subset(MP01, Marker=="Ch01h01")$Allele_2)
names(Ch01h01)[names(Ch01h01)=="subset.MP01..Marker.....Ch01h01...Sample_File_new"] <- "Sample_File"
names(Ch01h01)[names(Ch01h01)=="subset.MP01..Marker.....Ch01h01...Sample_Name_new"] <- "Sample_Name"
names(Ch01h01)[names(Ch01h01)=="subset.MP01..Marker.....Ch01h01...Allele_1"] <- "Ch01h01_A"
names(Ch01h01)[names(Ch01h01)=="subset.MP01..Marker.....Ch01h01...Allele_2"] <- "Ch01h01_B"
Ch01h01$subset.MP01..Marker.....Ch01h01...Marker <- NULL

merged_MP01_01 <- merge (Ch01f03b,Ch01h10)
merged_MP01 <- merge (merged_MP01_01,Ch01h01)
write.table (merged_MP01, file="merged_MP01.txt")

#MP02
MP02 <- read.table("MP02_r?colte2014_genotypes.txt", header=TRUE)
MP02$Sample_File_new <- gsub("MP02_", "", MP02$Sample_File)
MP02$Sample_Name_new <- gsub("MP02_", "", MP02$Sample_Name)

NZ05g08 <- data.frame(subset(MP02, Marker=="NZ05g08")$Sample_File_new, subset(MP02, Marker=="NZ05g08")$Sample_Name_new, subset(MP02, Marker=="NZ05g08")$Marker, subset(MP02, Marker=="NZ05g08")$Allele_1, subset(MP02, Marker=="NZ05g08")$Allele_2)
names(NZ05g08)[names(NZ05g08)=="subset.MP02..Marker.....NZ05g08...Sample_File_new"] <- "Sample_File"
names(NZ05g08)[names(NZ05g08)=="subset.MP02..Marker.....NZ05g08...Sample_Name_new"] <- "Sample_Name"
names(NZ05g08)[names(NZ05g08)=="subset.MP02..Marker.....NZ05g08...Allele_1"] <- "NZ05g08_A"
names(NZ05g08)[names(NZ05g08)=="subset.MP02..Marker.....NZ05g08...Allele_2"] <- "NZ05g08_B"
NZ05g08$subset.MP02..Marker.....NZ05g08...Marker <- NULL

CH05f06 <- data.frame(subset(MP02, Marker=="CH05f06")$Sample_File_new, subset(MP02, Marker=="CH05f06")$Sample_Name_new, subset(MP02, Marker=="CH05f06")$Marker, subset(MP02, Marker=="CH05f06")$Allele_1, subset(MP02, Marker=="CH05f06")$Allele_2)
names(CH05f06)[names(CH05f06)=="subset.MP02..Marker.....CH05f06...Sample_File_new"] <- "Sample_File"
names(CH05f06)[names(CH05f06)=="subset.MP02..Marker.....CH05f06...Sample_Name_new"] <- "Sample_Name"
names(CH05f06)[names(CH05f06)=="subset.MP02..Marker.....CH05f06...Allele_1"] <- "CH05f06_A"
names(CH05f06)[names(CH05f06)=="subset.MP02..Marker.....CH05f06...Allele_2"] <- "CH05f06_B"
CH05f06$subset.MP02..Marker.....CH05f06...Marker <- NULL

Ch02d08 <- data.frame(subset(MP02, Marker=="Ch02d08")$Sample_File_new, subset(MP02, Marker=="Ch02d08")$Sample_Name_new, subset(MP02, Marker=="Ch02d08")$Marker, subset(MP02, Marker=="Ch02d08")$Allele_1, subset(MP02, Marker=="Ch02d08")$Allele_2)
names(Ch02d08)[names(Ch02d08)=="subset.MP02..Marker.....Ch02d08...Sample_File_new"] <- "Sample_File"
names(Ch02d08)[names(Ch02d08)=="subset.MP02..Marker.....Ch02d08...Sample_Name_new"] <- "Sample_Name"
names(Ch02d08)[names(Ch02d08)=="subset.MP02..Marker.....Ch02d08...Allele_1"] <- "Ch02d08_A"
names(Ch02d08)[names(Ch02d08)=="subset.MP02..Marker.....Ch02d08...Allele_2"] <- "Ch02d08_B"
Ch02d08$subset.MP02..Marker.....Ch02d08...Marker <- NULL

Ch04e05 <- data.frame(subset(MP02, Marker=="Ch04e05")$Sample_File_new, subset(MP02, Marker=="Ch04e05")$Sample_Name_new, subset(MP02, Marker=="Ch04e05")$Marker, subset(MP02, Marker=="Ch04e05")$Allele_1, subset(MP02, Marker=="Ch04e05")$Allele_2)
names(Ch04e05)[names(Ch04e05)=="subset.MP02..Marker.....Ch04e05...Sample_File_new"] <- "Sample_File"
names(Ch04e05)[names(Ch04e05)=="subset.MP02..Marker.....Ch04e05...Sample_Name_new"] <- "Sample_Name"
names(Ch04e05)[names(Ch04e05)=="subset.MP02..Marker.....Ch04e05...Allele_1"] <- "Ch04e05_A"
names(Ch04e05)[names(Ch04e05)=="subset.MP02..Marker.....Ch04e05...Allele_2"] <- "Ch04e05_B"
Ch04e05$subset.MP02..Marker.....Ch04e05...Marker <- NULL

merged_MP02_1 <- merge (NZ05g08,CH05f06)
merged_MP02_2 <- merge (merged_MP02_1,Ch02d08)
merged_MP02 <- merge (merged_MP02_2,Ch04e05)
write.table (merged_MP02, file="merged_MP02.txt")

#MP03
MP03 <- read.table("MP03_r?colte2014_genotypes.txt", header=TRUE)
MP03$Sample_File_new <- gsub("MP03_", "", MP03$Sample_File)
MP03$Sample_Name_new <- gsub("MP03_", "", MP03$Sample_Name)

Hi02c07 <- data.frame(subset(MP03, Marker=="Hi02c07")$Sample_File_new, subset(MP03, Marker=="Hi02c07")$Sample_Name_new, subset(MP03, Marker=="Hi02c07")$Marker, subset(MP03, Marker=="Hi02c07")$Allele_1, subset(MP03, Marker=="Hi02c07")$Allele_2)
names(Hi02c07)[names(Hi02c07)=="subset.MP03..Marker.....Hi02c07...Sample_File_new"] <- "Sample_File"
names(Hi02c07)[names(Hi02c07)=="subset.MP03..Marker.....Hi02c07...Sample_Name_new"] <- "Sample_Name"
names(Hi02c07)[names(Hi02c07)=="subset.MP03..Marker.....Hi02c07...Allele_1"] <- "Hi02c07_A"
names(Hi02c07)[names(Hi02c07)=="subset.MP03..Marker.....Hi02c07...Allele_2"] <- "Hi02c07_B"
Hi02c07$subset.MP03..Marker.....Hi02c07...Marker <- NULL

Ch01f02 <- data.frame(subset(MP03, Marker=="Ch01f02")$Sample_File_new, subset(MP03, Marker=="Ch01f02")$Sample_Name_new, subset(MP03, Marker=="Ch01f02")$Marker, subset(MP03, Marker=="Ch01f02")$Allele_1, subset(MP03, Marker=="Ch01f02")$Allele_2)
names(Ch01f02)[names(Ch01f02)=="subset.MP03..Marker.....Ch01f02...Sample_File_new"] <- "Sample_File"
names(Ch01f02)[names(Ch01f02)=="subset.MP03..Marker.....Ch01f02...Sample_Name_new"] <- "Sample_Name"
names(Ch01f02)[names(Ch01f02)=="subset.MP03..Marker.....Ch01f02...Allele_1"] <- "Ch01f02_A"
names(Ch01f02)[names(Ch01f02)=="subset.MP03..Marker.....Ch01f02...Allele_2"] <- "Ch01f02_B"
Ch01f02$subset.MP03..Marker.....Ch01f02...Marker <- NULL

Ch02c11 <- data.frame(subset(MP03, Marker=="Ch02c11")$Sample_File_new, subset(MP03, Marker=="Ch02c11")$Sample_Name_new, subset(MP03, Marker=="Ch02c11")$Marker, subset(MP03, Marker=="Ch02c11")$Allele_1, subset(MP03, Marker=="Ch02c11")$Allele_2)
names(Ch02c11)[names(Ch02c11)=="subset.MP03..Marker.....Ch02c11...Sample_File_new"] <- "Sample_File"
names(Ch02c11)[names(Ch02c11)=="subset.MP03..Marker.....Ch02c11...Sample_Name_new"] <- "Sample_Name"
names(Ch02c11)[names(Ch02c11)=="subset.MP03..Marker.....Ch02c11...Allele_1"] <- "Ch02c11_A"
names(Ch02c11)[names(Ch02c11)=="subset.MP03..Marker.....Ch02c11...Allele_2"] <- "Ch02c11_B"
Ch02c11$subset.MP03..Marker.....Ch02c11...Marker <- NULL

merged_MP03_1 <- merge (Hi02c07,Ch01f02)
merged_MP03 <- merge (merged_MP03_1,Ch02c11)
write.table (merged_MP03, file="merged_MP03.txt")

#MP04
MP04 <- read.table("MP04_r?colte2014_genotypes.txt", header=TRUE)
MP04$Sample_File_new <- gsub("MP04_", "", MP04$Sample_File)
MP04$Sample_Name_new <- gsub("MP04_", "", MP04$Sample_Name)

CH04c07 <- data.frame(subset(MP04, Marker=="CH04c07")$Sample_File_new, subset(MP04, Marker=="CH04c07")$Sample_Name_new, subset(MP04, Marker=="CH04c07")$Marker, subset(MP04, Marker=="CH04c07")$Allele_1, subset(MP04, Marker=="CH04c07")$Allele_2)
names(CH04c07)[names(CH04c07)=="subset.MP04..Marker.....CH04c07...Sample_File_new"] <- "Sample_File"
names(CH04c07)[names(CH04c07)=="subset.MP04..Marker.....CH04c07...Sample_Name_new"] <- "Sample_Name"
names(CH04c07)[names(CH04c07)=="subset.MP04..Marker.....CH04c07...Allele_1"] <- "CH04c07_A"
names(CH04c07)[names(CH04c07)=="subset.MP04..Marker.....CH04c07...Allele_2"] <- "CH04c07_B"
CH04c07$subset.MP04..Marker.....CH04c07...Marker <- NULL

GD12 <- data.frame(subset(MP04, Marker=="GD12")$Sample_File_new, subset(MP04, Marker=="GD12")$Sample_Name_new, subset(MP04, Marker=="GD12")$Marker, subset(MP04, Marker=="GD12")$Allele_1, subset(MP04, Marker=="GD12")$Allele_2)
names(GD12)[names(GD12)=="subset.MP04..Marker.....GD12...Sample_File_new"] <- "Sample_File"
names(GD12)[names(GD12)=="subset.MP04..Marker.....GD12...Sample_Name_new"] <- "Sample_Name"
names(GD12)[names(GD12)=="subset.MP04..Marker.....GD12...Allele_1"] <- "GD12_A"
names(GD12)[names(GD12)=="subset.MP04..Marker.....GD12...Allele_2"] <- "GD12_B"
GD12$subset.MP04..Marker.....GD12...Marker <- NULL

CH03d07 <- data.frame(subset(MP04, Marker=="CH03d07")$Sample_File_new, subset(MP04, Marker=="CH03d07")$Sample_Name_new, subset(MP04, Marker=="CH03d07")$Marker, subset(MP04, Marker=="CH03d07")$Allele_1, subset(MP04, Marker=="CH03d07")$Allele_2)
names(CH03d07)[names(CH03d07)=="subset.MP04..Marker.....CH03d07...Sample_File_new"] <- "Sample_File"
names(CH03d07)[names(CH03d07)=="subset.MP04..Marker.....CH03d07...Sample_Name_new"] <- "Sample_Name"
names(CH03d07)[names(CH03d07)=="subset.MP04..Marker.....CH03d07...Allele_1"] <- "CH03d07_A"
names(CH03d07)[names(CH03d07)=="subset.MP04..Marker.....CH03d07...Allele_2"] <- "CH03d07_B"
CH03d07$subset.MP04..Marker.....CH03d07...Marker <- NULL

CH02c09 <- data.frame(subset(MP04, Marker=="CH02c09")$Sample_File_new, subset(MP04, Marker=="CH02c09")$Sample_Name_new, subset(MP04, Marker=="CH02c09")$Marker, subset(MP04, Marker=="CH02c09")$Allele_1, subset(MP04, Marker=="CH02c09")$Allele_2)
names(CH02c09)[names(CH02c09)=="subset.MP04..Marker.....CH02c09...Sample_File_new"] <- "Sample_File"
names(CH02c09)[names(CH02c09)=="subset.MP04..Marker.....CH02c09...Sample_Name_new"] <- "Sample_Name"
names(CH02c09)[names(CH02c09)=="subset.MP04..Marker.....CH02c09...Allele_1"] <- "CH02c09_A"
names(CH02c09)[names(CH02c09)=="subset.MP04..Marker.....CH02c09...Allele_2"] <- "CH02c09_B"
CH02c09$subset.MP04..Marker.....CH02c09...Marker <- NULL

merged_MP04_1 <- merge (CH04c07,GD12)
merged_MP04_2 <- merge (merged_MP04_1,CH03d07)
merged_MP04 <- merge (merged_MP04_2,CH02c09)
write.table (merged_MP04, file="merged_MP04.txt")

#Hi4b
Hi4b <- read.table("Hi4b_r?colte2014_genotypes.txt", header=TRUE)
Hi4b$Sample_File_new <- gsub("Hi4b_", "", Hi4b$Sample_File)
Hi4b$Sample_Name_new <- gsub("Hi4b_", "", Hi4b$Sample_Name)

CH02g01 <- data.frame(subset(Hi4b, Marker=="CH02g01")$Sample_File_new, subset(Hi4b, Marker=="CH02g01")$Sample_Name_new, subset(Hi4b, Marker=="CH02g01")$Marker, subset(Hi4b, Marker=="CH02g01")$Allele_1, subset(Hi4b, Marker=="CH02g01")$Allele_2)
names(CH02g01)[names(CH02g01)=="subset.Hi4b..Marker.....CH02g01...Sample_File_new"] <- "Sample_File"
names(CH02g01)[names(CH02g01)=="subset.Hi4b..Marker.....CH02g01...Sample_Name_new"] <- "Sample_Name"
names(CH02g01)[names(CH02g01)=="subset.Hi4b..Marker.....CH02g01...Allele_1"] <- "CH02g01_A"
names(CH02g01)[names(CH02g01)=="subset.Hi4b..Marker.....CH02g01...Allele_2"] <- "CH02g01_B"
CH02g01$subset.Hi4b..Marker.....CH02g01...Marker <- NULL

Hi23g02 <- data.frame(subset(Hi4b, Marker=="Hi23g02")$Sample_File_new, subset(Hi4b, Marker=="Hi23g02")$Sample_Name_new, subset(Hi4b, Marker=="Hi23g02")$Marker, subset(Hi4b, Marker=="Hi23g02")$Allele_1, subset(Hi4b, Marker=="Hi23g02")$Allele_2)
names(Hi23g02)[names(Hi23g02)=="subset.Hi4b..Marker.....Hi23g02...Sample_File_new"] <- "Sample_File"
names(Hi23g02)[names(Hi23g02)=="subset.Hi4b..Marker.....Hi23g02...Sample_Name_new"] <- "Sample_Name"
names(Hi23g02)[names(Hi23g02)=="subset.Hi4b..Marker.....Hi23g02...Allele_1"] <- "Hi23g02_A"
names(Hi23g02)[names(Hi23g02)=="subset.Hi4b..Marker.....Hi23g02...Allele_2"] <- "Hi23g02_B"
Hi23g02$subset.Hi4b..Marker.....Hi23g02...Marker <- NULL

merged_Hi4b <- merge (CH02g01,Hi23g02)
write.table (merged_Hi4b, file="merged_Hi4b.txt")

#Hi4a
Hi4a <- read.table("Hi4a_r?colte2014_genotypes.txt", header=TRUE)
Hi4a$Sample_File_new <- gsub("Hi4a_", "", Hi4a$Sample_File)
Hi4a$Sample_Name_new <- gsub("Hi4a_", "", Hi4a$Sample_Name)

CH02c02b <- data.frame(subset(Hi4a, Marker=="CH02c02b")$Sample_File_new, subset(Hi4a, Marker=="CH02c02b")$Sample_Name_new, subset(Hi4a, Marker=="CH02c02b")$Marker, subset(Hi4a, Marker=="CH02c02b")$Allele_1, subset(Hi4a, Marker=="CH02c02b")$Allele_2)
names(CH02c02b)[names(CH02c02b)=="subset.Hi4a..Marker.....CH02c02b...Sample_File_new"] <- "Sample_File"
names(CH02c02b)[names(CH02c02b)=="subset.Hi4a..Marker.....CH02c02b...Sample_Name_new"] <- "Sample_Name"
names(CH02c02b)[names(CH02c02b)=="subset.Hi4a..Marker.....CH02c02b...Allele_1"] <- "CH02c02b_A"
names(CH02c02b)[names(CH02c02b)=="subset.Hi4a..Marker.....CH02c02b...Allele_2"] <- "CH02c02b_B"
CH02c02b$subset.Hi4a..Marker.....CH02c02b...Marker <- NULL

CH04e03 <- data.frame(subset(Hi4a, Marker=="CH04e03")$Sample_File_new, subset(Hi4a, Marker=="CH04e03")$Sample_Name_new, subset(Hi4a, Marker=="CH04e03")$Marker, subset(Hi4a, Marker=="CH04e03")$Allele_1, subset(Hi4a, Marker=="CH04e03")$Allele_2)
names(CH04e03)[names(CH04e03)=="subset.Hi4a..Marker.....CH04e03...Sample_File_new"] <- "Sample_File"
names(CH04e03)[names(CH04e03)=="subset.Hi4a..Marker.....CH04e03...Sample_Name_new"] <- "Sample_Name"
names(CH04e03)[names(CH04e03)=="subset.Hi4a..Marker.....CH04e03...Allele_1"] <- "CH04e03_A"
names(CH04e03)[names(CH04e03)=="subset.Hi4a..Marker.....CH04e03...Allele_2"] <- "CH04e03_B"
CH04e03$subset.Hi4a..Marker.....CH04e03...Marker <- NULL

merged_Hi4a <- merge (CH02c02b,CH04e03)
write.table (merged_Hi4a, file="merged_Hi4a.txt")

#Hi5_10
Hi5_10 <- read.table("Hi5_10_r?colte2014_genotypes.txt", header=TRUE)
Hi5_10$Sample_File_new <- gsub("Hi5_10_", "", Hi5_10$Sample_File)
Hi5_10$Sample_Name_new <- gsub("Hi5_10_", "", Hi5_10$Sample_Name)

CH05e06 <- data.frame(subset(Hi5_10, Marker=="CH05e06")$Sample_File_new, subset(Hi5_10, Marker=="CH05e06")$Sample_Name_new, subset(Hi5_10, Marker=="CH05e06")$Marker, subset(Hi5_10, Marker=="CH05e06")$Allele_1, subset(Hi5_10, Marker=="CH05e06")$Allele_2)
names(CH05e06)[names(CH05e06)=="subset.Hi5_10..Marker.....CH05e06...Sample_File_new"] <- "Sample_File"
names(CH05e06)[names(CH05e06)=="subset.Hi5_10..Marker.....CH05e06...Sample_Name_new"] <- "Sample_Name"
names(CH05e06)[names(CH05e06)=="subset.Hi5_10..Marker.....CH05e06...Allele_1"] <- "CH05e06_A"
names(CH05e06)[names(CH05e06)=="subset.Hi5_10..Marker.....CH05e06...Allele_2"] <- "CH05e06_B"
CH05e06$subset.Hi5_10..Marker.....CH05e06...Marker <- NULL

CH02b12 <- data.frame(subset(Hi5_10, Marker=="CH02b12")$Sample_File_new, subset(Hi5_10, Marker=="CH02b12")$Sample_Name_new, subset(Hi5_10, Marker=="CH02b12")$Marker, subset(Hi5_10, Marker=="CH02b12")$Allele_1, subset(Hi5_10, Marker=="CH02b12")$Allele_2)
names(CH02b12)[names(CH02b12)=="subset.Hi5_10..Marker.....CH02b12...Sample_File_new"] <- "Sample_File"
names(CH02b12)[names(CH02b12)=="subset.Hi5_10..Marker.....CH02b12...Sample_Name_new"] <- "Sample_Name"
names(CH02b12)[names(CH02b12)=="subset.Hi5_10..Marker.....CH02b12...Allele_1"] <- "CH02b12_A"
names(CH02b12)[names(CH02b12)=="subset.Hi5_10..Marker.....CH02b12...Allele_2"] <- "CH02b12_B"
CH02b12$subset.Hi5_10..Marker.....CH02b12...Marker <- NULL

Hi03a10 <- data.frame(subset(Hi5_10, Marker=="Hi03a10")$Sample_File_new, subset(Hi5_10, Marker=="Hi03a10")$Sample_Name_new, subset(Hi5_10, Marker=="Hi03a10")$Marker, subset(Hi5_10, Marker=="Hi03a10")$Allele_1, subset(Hi5_10, Marker=="Hi03a10")$Allele_2)
names(Hi03a10)[names(Hi03a10)=="subset.Hi5_10..Marker.....Hi03a10...Sample_File_new"] <- "Sample_File"
names(Hi03a10)[names(Hi03a10)=="subset.Hi5_10..Marker.....Hi03a10...Sample_Name_new"] <- "Sample_Name"
names(Hi03a10)[names(Hi03a10)=="subset.Hi5_10..Marker.....Hi03a10...Allele_1"] <- "Hi03a10_A"
names(Hi03a10)[names(Hi03a10)=="subset.Hi5_10..Marker.....Hi03a10...Allele_2"] <- "Hi03a10_B"
Hi03a10$subset.Hi5_10..Marker.....Hi03a10...Marker <- NULL

MS06g03 <- data.frame(subset(Hi5_10, Marker=="MS06g03")$Sample_File_new, subset(Hi5_10, Marker=="MS06g03")$Sample_Name_new, subset(Hi5_10, Marker=="MS06g03")$Marker, subset(Hi5_10, Marker=="MS06g03")$Allele_1, subset(Hi5_10, Marker=="MS06g03")$Allele_2)
names(MS06g03)[names(MS06g03)=="subset.Hi5_10..Marker.....MS06g03...Sample_File_new"] <- "Sample_File"
names(MS06g03)[names(MS06g03)=="subset.Hi5_10..Marker.....MS06g03...Sample_Name_new"] <- "Sample_Name"
names(MS06g03)[names(MS06g03)=="subset.Hi5_10..Marker.....MS06g03...Allele_1"] <- "MS06g03_A"
names(MS06g03)[names(MS06g03)=="subset.Hi5_10..Marker.....MS06g03...Allele_2"] <- "MS06g03_B"
MS06g03$subset.Hi5_10..Marker.....MS06g03...Marker <- NULL

CH02b03b <- data.frame(subset(Hi5_10, Marker=="CH02b03b")$Sample_File_new, subset(Hi5_10, Marker=="CH02b03b")$Sample_Name_new, subset(Hi5_10, Marker=="CH02b03b")$Marker, subset(Hi5_10, Marker=="CH02b03b")$Allele_1, subset(Hi5_10, Marker=="CH02b03b")$Allele_2)
names(CH02b03b)[names(CH02b03b)=="subset.Hi5_10..Marker.....CH02b03b...Sample_File_new"] <- "Sample_File"
names(CH02b03b)[names(CH02b03b)=="subset.Hi5_10..Marker.....CH02b03b...Sample_Name_new"] <- "Sample_Name"
names(CH02b03b)[names(CH02b03b)=="subset.Hi5_10..Marker.....CH02b03b...Allele_1"] <- "CH02b03b_A"
names(CH02b03b)[names(CH02b03b)=="subset.Hi5_10..Marker.....CH02b03b...Allele_2"] <- "CH02b03b_B"
CH02b03b$subset.Hi5_10..Marker.....CH02b03b...Marker <- NULL

merged_Hi5_10_01 <- merge (CH05e06,CH02b12)
merged_Hi5_10_02 <- merge (merged_Hi5_10_01,Hi03a10)
merged_Hi5_10_03 <- merge (merged_Hi5_10_02,MS06g03)
merged_Hi5_10 <- merge (merged_Hi5_10_03,CH02b03b)

write.table (merged_Hi5_10, file="merged_Hi5_10.txt")

#Hi13b
Hi13b <- read.table("Hi13b_r?colte_2014_genotypes.txt", header=TRUE)
Hi13b$Sample_File_new <- gsub("Hi13b_", "", Hi13b$Sample_File)
Hi13b$Sample_Name_new <- gsub("Hi13b_", "", Hi13b$Sample_Name)

CH05f04 <- data.frame(subset(Hi13b, Marker=="CH05f04")$Sample_File_new, subset(Hi13b, Marker=="CH05f04")$Sample_Name_new, subset(Hi13b, Marker=="CH05f04")$Marker, subset(Hi13b, Marker=="CH05f04")$Allele_1, subset(Hi13b, Marker=="CH05f04")$Allele_2)
names(CH05f04)[names(CH05f04)=="subset.Hi13b..Marker.....CH05f04...Sample_File_new"] <- "Sample_File"
names(CH05f04)[names(CH05f04)=="subset.Hi13b..Marker.....CH05f04...Sample_Name_new"] <- "Sample_Name"
names(CH05f04)[names(CH05f04)=="subset.Hi13b..Marker.....CH05f04...Allele_1"] <- "CH05f04_A"
names(CH05f04)[names(CH05f04)=="subset.Hi13b..Marker.....CH05f04...Allele_2"] <- "CH05f04_B"
CH05f04$subset.Hi13b..Marker.....CH05f04...Marker <- NULL

GD103 <- data.frame(subset(Hi13b, Marker=="GD103")$Sample_File_new, subset(Hi13b, Marker=="GD103")$Sample_Name_new, subset(Hi13b, Marker=="GD103")$Marker, subset(Hi13b, Marker=="GD103")$Allele_1, subset(Hi13b, Marker=="GD103")$Allele_2)
names(GD103)[names(GD103)=="subset.Hi13b..Marker.....GD103...Sample_File_new"] <- "Sample_File"
names(GD103)[names(GD103)=="subset.Hi13b..Marker.....GD103...Sample_Name_new"] <- "Sample_Name"
names(GD103)[names(GD103)=="subset.Hi13b..Marker.....GD103...Allele_1"] <- "GD103_A"
names(GD103)[names(GD103)=="subset.Hi13b..Marker.....GD103...Allele_2"] <- "GD103_B"
GD103$subset.Hi13b..Marker.....GD103...Marker <- NULL

NH009b <- data.frame(subset(Hi13b, Marker=="NH009b")$Sample_File_new, subset(Hi13b, Marker=="NH009b")$Sample_Name_new, subset(Hi13b, Marker=="NH009b")$Marker, subset(Hi13b, Marker=="NH009b")$Allele_1, subset(Hi13b, Marker=="NH009b")$Allele_2)
names(NH009b)[names(NH009b)=="subset.Hi13b..Marker.....NH009b...Sample_File_new"] <- "Sample_File"
names(NH009b)[names(NH009b)=="subset.Hi13b..Marker.....NH009b...Sample_Name_new"] <- "Sample_Name"
names(NH009b)[names(NH009b)=="subset.Hi13b..Marker.....NH009b...Allele_1"] <- "NH009b_A"
names(NH009b)[names(NH009b)=="subset.Hi13b..Marker.....NH009b...Allele_2"] <- "NH009b_B"
NH009b$subset.Hi13b..Marker.....NH009b...Marker <- NULL

CH05g08 <- data.frame(subset(Hi13b, Marker=="CH05g08")$Sample_File_new, subset(Hi13b, Marker=="CH05g08")$Sample_Name_new, subset(Hi13b, Marker=="CH05g08")$Marker, subset(Hi13b, Marker=="CH05g08")$Allele_1, subset(Hi13b, Marker=="CH05g08")$Allele_2)
names(CH05g08)[names(CH05g08)=="subset.Hi13b..Marker.....CH05g08...Sample_File_new"] <- "Sample_File"
names(CH05g08)[names(CH05g08)=="subset.Hi13b..Marker.....CH05g08...Sample_Name_new"] <- "Sample_Name"
names(CH05g08)[names(CH05g08)=="subset.Hi13b..Marker.....CH05g08...Allele_1"] <- "CH05g08_A"
names(CH05g08)[names(CH05g08)=="subset.Hi13b..Marker.....CH05g08...Allele_2"] <- "CH05g08_B"
CH05g08$subset.Hi13b..Marker.....CH05g08...Marker <- NULL

merged_Hi13b_01 <- merge (CH05f04,GD103)
merged_Hi13b_02 <- merge (merged_Hi13b_01,NH009b)
merged_Hi13b <- merge (merged_Hi13b_02,CH05g08)

write.table (merged_Hi13b, file="merged_Hi13b.txt")

#Hi6
Hi6 <- read.table("Hi6_r?colte2014_genotypes.txt", header=TRUE)
Hi6$Sample_File_new <- gsub("Hi6_", "", Hi6$Sample_File)
Hi6$Sample_Name_new <- gsub("Hi6_", "", Hi6$Sample_Name)

Ch05a05 <- data.frame(subset(Hi6, Marker=="Ch05a05")$Sample_File_new, subset(Hi6, Marker=="Ch05a05")$Sample_Name_new, subset(Hi6, Marker=="Ch05a05")$Marker, subset(Hi6, Marker=="Ch05a05")$Allele_1, subset(Hi6, Marker=="Ch05a05")$Allele_2)
names(Ch05a05)[names(Ch05a05)=="subset.Hi6..Marker.....Ch05a05...Sample_File_new"] <- "Sample_File"
names(Ch05a05)[names(Ch05a05)=="subset.Hi6..Marker.....Ch05a05...Sample_Name_new"] <- "Sample_Name"
names(Ch05a05)[names(Ch05a05)=="subset.Hi6..Marker.....Ch05a05...Allele_1"] <- "Ch05a05_A"
names(Ch05a05)[names(Ch05a05)=="subset.Hi6..Marker.....Ch05a05...Allele_2"] <- "Ch05a05_B"
Ch05a05$subset.Hi6..Marker.....Ch05a05...Marker <- NULL

CH_Vf1 <- data.frame(subset(Hi6, Marker=="CH_Vf1")$Sample_File_new, subset(Hi6, Marker=="CH_Vf1")$Sample_Name_new, subset(Hi6, Marker=="CH_Vf1")$Marker, subset(Hi6, Marker=="CH_Vf1")$Allele_1, subset(Hi6, Marker=="CH_Vf1")$Allele_2)
names(CH_Vf1)[names(CH_Vf1)=="subset.Hi6..Marker.....CH_Vf1...Sample_File_new"] <- "Sample_File"
names(CH_Vf1)[names(CH_Vf1)=="subset.Hi6..Marker.....CH_Vf1...Sample_Name_new"] <- "Sample_Name"
names(CH_Vf1)[names(CH_Vf1)=="subset.Hi6..Marker.....CH_Vf1...Allele_1"] <- "CH_Vf1_A"
names(CH_Vf1)[names(CH_Vf1)=="subset.Hi6..Marker.....CH_Vf1...Allele_2"] <- "CH_Vf1_B"
CH_Vf1$subset.Hi6..Marker.....CH_Vf1...Marker <- NULL

CH03d12 <- data.frame(subset(Hi6, Marker=="CH03d12")$Sample_File_new, subset(Hi6, Marker=="CH03d12")$Sample_Name_new, subset(Hi6, Marker=="CH03d12")$Marker, subset(Hi6, Marker=="CH03d12")$Allele_1, subset(Hi6, Marker=="CH03d12")$Allele_2)
names(CH03d12)[names(CH03d12)=="subset.Hi6..Marker.....CH03d12...Sample_File_new"] <- "Sample_File"
names(CH03d12)[names(CH03d12)=="subset.Hi6..Marker.....CH03d12...Sample_Name_new"] <- "Sample_Name"
names(CH03d12)[names(CH03d12)=="subset.Hi6..Marker.....CH03d12...Allele_1"] <- "CH03d12_A"
names(CH03d12)[names(CH03d12)=="subset.Hi6..Marker.....CH03d12...Allele_2"] <- "CH03d12_B"
CH03d12$subset.Hi6..Marker.....CH03d12...Marker <- NULL

merged_Hi6_01 <- merge (Ch05a05,CH_Vf1)
merged_Hi6 <- merge (merged_Hi6_01,CH03d12)

write.table (merged_Hi6, file="merged_Hi6.txt")

merged_all_01 <- merge (merged_Hi4a, merged_Hi13b)
merged_all_02 <- merge (merged_Hi4b, merged_all_01)
merged_all_03 <- merge (merged_Hi5_10, merged_all_02)
merged_all_04 <- merge (merged_Hi6, merged_all_03)
merged_all_05 <- merge (merged_MP01, merged_all_04)
merged_all_06 <- merge (merged_MP02, merged_all_05)
merged_all_07 <- merge (merged_MP03, merged_all_06)
merged_all_08 <- merge (merged_MP04, merged_all_07)
write.table (merged_all_08, file="dataset_r?colte2014.txt")
