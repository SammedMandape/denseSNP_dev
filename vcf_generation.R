library(tidyverse)
library(data.table)

# complement funciton
complement <- function(x){
  case_when(
    x=="A" ~ "T",
    x=="T" ~ "A",
    x=="G" ~ "C",
    x=="C" ~ "G"
  )
}
fread("CEPH_chips1_4_full_data.final_results.txt") -> mydataIn
mydataIn %>% group_by(`Sample ID`) %>% group_split() -> mydataIn_1.1

group_keys(mydataIn %>% group_by(`Sample ID`))
#mydataIn_1.1[[1]] %>% head(n=1000) -> x
rm(mydataIn, mydataIn_1)
rm(mydataIn)

#mydataIn %>% group_by(Position) %>% filter(n()>4) %>% arrange(Position) -> bar

#mydataIn_1.1 %>% colnames() %>% as_tibble() -> foo

mydataIn_analy <- function (x){
SampleNameID <-  sym(x %>% select(`Sample ID`) %>% distinct() %>% pull)
mydataIn_1.2 <- x %>% 
  ungroup() %>%
  filter(!(TrueAllele1_Forward == "") | !(TrueAllele2_Forward == "")) %>%
  separate(SNP, into = c("Snp1","Snp2"), sep = "/", remove=FALSE) %>% 
  separate(Snp1, into = c(NA,"Snp1"), sep =-1) %>%
  separate(Snp2, into = c("Snp2", NA), sep = -1) %>%
  mutate(Alt=case_when(ref_base==Snp1 ~ Snp2,
                       ref_base==Snp2 ~ Snp1,
                       ref_base==complement(Snp1) ~ complement(Snp2),
                       ref_base==complement(Snp2) ~ complement(Snp1),
                       ),
         QUAL=".",FILTER=".",INFO=".",
         GT=case_when(ref_base==TrueAllele1_Forward & ref_base== TrueAllele2_Forward ~ "0/0",
                      (ref_base==TrueAllele1_Forward & Alt==TrueAllele2_Forward) | (ref_base==TrueAllele2_Forward & Alt==TrueAllele1_Forward) ~ "0/1",
                      Alt==TrueAllele1_Forward & Alt==TrueAllele2_Forward ~ "1/1",
                      TRUE ~ "Something is wrong"
                      )
                    ) %>% 
  select(Chr, Position, `SNP Name`, ref_base, Alt, QUAL, FILTER, INFO, 
         GT, SNP, TrueAllele1_Forward, TrueAllele2_Forward, `Allele1 - AB`,`Allele2 - AB`,
         `GC Score`, `GT Score`, Theta, R, X, Y, `X Raw`, `Y Raw`, `Log R Ratio`
         ) 

# extract colnames from given user input and after INFO col
# remove any white space in the col name that will go into FORMAT column
# paste or unite the extracted columns into sample id column. 
ncol(mydataIn_1.2) -> totCol
fromCol <- totCol - 8

mydataIn_1.2[,c(9:totCol)] %>% colnames() -> reqCols
reqCols %>% gsub(" ","",.)->reqCols_1

mydataIn_1.3 <- mydataIn_1.2 %>% 
  mutate(FORMAT = paste(paste(reqCols_1,collapse = ":"),"GL",sep = ":")) %>% 
  unite(!!SampleNameID, reqCols,sep=":") %>% 
  select(Chr:INFO,FORMAT,!!SampleNameID) %>%
  rename("#CHROM"="Chr","POS"="Position","ID"="SNP Name","REF"="ref_base","ALT"="Alt")
  #colnames(paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",!!SampleNameID)))
  #select() -> foo
  #mutate(SampleNameID = paste(reqCols,sep = ";")) -> foo
header<-"##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=SNP,Number=1,Type=String,Description=\"Observed snp\">
##FORMAT=<ID=TrueAllele1_Forward,Number=1,Type=Character,Description=\"Allele 1 reported on forward strand\">
##FORMAT=<ID=TrueAllele2_Forward,Number=1,Type=Character,Description=\"Allele 2 reported on forward strand\">
##FORMAT=<ID=Allele1-AB,Number=1,Type=Character,Description=\"Allele 1 A or B call\">
##FORMAT=<ID=Allele2-AB,Number=1,Type=Character,Description=\"Allele 2 A or B call\">
##FORMAT=<ID=GCScore,Number=1,Type=Float,Description=\"GenCall, a confidence measure assigned to each call\">
##FORMAT=<ID=GTScore,Number=1,Type=Float,Description=\"Measure of the cluster quality for the SNP\">
##FORMAT=<ID=Theta,Number=1,Type=Float,Description=\"Normalized Theta-value for the sample\">
##FORMAT=<ID=R,Number=1,Type=Float,Description=\"Normalized R-value for the sample\">
##FORMAT=<ID=X,Number=1,Type=Float,Description=\"Normalized intensity of the A allele\">
##FORMAT=<ID=Y,Number=1,Type=Float,Description=\"Normalized intensity of the B allele\">
##FORMAT=<ID=XRaw,Number=1,Type=Integer,Description=\"Raw intensity of the A allele\">
##FORMAT=<ID=YRaw,Number=1,Type=Integer,Description=\"Raw intensity of the B allele\">
##FORMAT=<ID=LogRRatio,Number=1,Type=Float,Description=\"Base-2 log of the normalized R value over the expected R value for the theta value\">"
#########################################cat(header, sep="\n",file=paste0(SampleNameID,".txt"))
write.table(as_tibble(header),paste0(SampleNameID,".txt"),quote=F, col.names = F, row.names = F)
write_tsv(mydataIn_1.3,paste0(SampleNameID,".txt"),append = T, col_names = T)
print("1")
}

map(mydataIn_1.1, mydataIn_analy)

# sanity checks
mydataIn_1.2 %>% filter(grepl("[^ACGT]",ref_base)) %>% nrow() # should be 0
mydataIn_1.2 %>% filter(grepl("[^ACGT]",Alt)) %>% nrow() # should be 0

# test==1 should be 0, this can be used to get the # of homo, alt homo, and hetero
mydataIn_1.2 %>% mutate(test=ifelse(GT=="0/0" & TrueAllele1_Forward==TrueAllele2_Forward & `Allele1 - AB`==`Allele2 - AB`,0,
                                    ifelse(GT=="1/1" & TrueAllele1_Forward==TrueAllele2_Forward & `Allele1 - AB`==`Allele2 - AB`,-1,
                                           ifelse(`Allele1 - AB`=="A" & `Allele2 - AB`=="B",-2,1)))) %>% filter(test==1) %>% nrow() 
#mydataIn %>% filter(`Allele1 - AB`=="B" & `Allele2 - AB`=="A")
