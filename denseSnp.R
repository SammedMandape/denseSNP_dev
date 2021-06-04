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


#mytidydata_input <- read_csv("CEPH_chips1_4_full_data.table.txt.gz_tidy_data.txt", col_names = TRUE)
#mytidydata_input <- fread("CEPH_chips1_4_full_data.table.txt.gz_tidy_data.txt")
#mytidydata_input <- fread("CEPH_chips5-8_FinalReport.txt_tidy_data.txt")
mytidydata_input <- fread("CEPH_chips_9-12_FinalReport.txt_tidy_data.txt")

#system.time({mytidydata_input <- fread("CEPH_chips1_4_full_data.table.txt.gz_tidy_data.txt")})

#myreffasta<-read_delim("regions_variants_ref_oneline.fasta", col_names = c('region_samtools','ref_base'), delim = "\t")
myreffasta<-read_delim("regions_variants_chips_9_12_ref.fastaoneline.fasta", col_names = c('region_samtools','ref_base'), delim = "\t")
setDT(myreffasta)

# a problem with fread, it skips the first row
#myreffasta<-fread("regions_variants_ref_oneline.fasta", col.names = c('region_samtools','ref_base'))
#myreffasta<-fread("regions_variants_chips_5-8_ref.fastaoneline.fasta", col.names = c('region_samtools','ref_base'))

#setkey(myreffasta)

myreffasta_1 <- myreffasta %>% 
  separate(region_samtools, into = c('foo','region_samtools'), remove = TRUE, sep = ">") %>% 
  select(-foo)

myreffasta_2 <- myreffasta_1 %>% distinct()

# following is an intermediate step as the size of the objects are too big
mytidydata_input_join <- mytidydata_input %>% 
  left_join(myreffasta_2, by=c("region_samtools"="region_samtools")) %>% 
  distinct()


# skip if you don't want intermediate outputs
#write_csv(mytidydata_input_join,path = "mytidydata_join.csv")
#mytidydata_foo <- read_delim("mytidydata_join.csv", delim=",")
#mytidydata <- fread("mytidydata_join.csv")

# following is an intermediate step as the size of the objects are too big
mytidydata <- mytidydata_input_join

# some clean up to emtpy space
rm(mytidydata_input_join,myreffasta,mytidydata_input,myreffasta_1,myreffasta_2)

################################################################################
## Testing
# mytidydata %>% head()
# 
# #setDT(mytidydata)[,]
# 
# mytidydata %>% 
#   group_by(`Sample ID`, Position, Chr) %>% 
#   summarise(dup = n()) %>% 
#   arrange(desc(dup))
# 
# mytidydata_500K %>% 
#   group_by(`Sample ID`, Position, Chr) %>% 
#   mutate(dup = n()) -> mytidydata_1
# 
# mytidydata_1 %>% 
#   filter(dup > 1) -> mytidydata_2
# 
# # to see if there are any dup snps
# mytidydata_2 %>% ungroup() %>% 
#   group_by(`Sample ID`,Position, Chr) %>% 
#   arrange(Chr, Position, desc(dup)) -> mytidydata_3
# 
# # maybe this info is useful?
# infiStrand <- read_delim("InfiniumOmni5-4v1-2_A1-b37.Ilmn.strand", delim = "\t",  col_names = F)

#################################################################################

#setDT(mytidydata)

# mytidydata %>% head(n=1000000) -> mytidydata_500K # test data set

mytidydata %>% 
  select(-region_samtools) %>%
  separate(SNP, into = c("foo1","foo2"),sep = "[/]", remove = F) %>% 
  separate(foo1, into = c(NA, "Snp1"),sep = -1) %>% 
  separate(foo2, into = c("Snp2", NA), sep=-1) -> mytidydata

# get the forward alleles comparing it to the reference allele
#system.time(
mytidydata %>% 
  #select(`SNP Name`:Position,SNP:`Plus/Minus Strand`,ref_base) %>% #filtered for test purpose
  mutate(FrwdOrRevCalculation=ifelse(ref_base == `Allele1 - Forward` | ref_base == `Allele2 - Forward`, 
                                    "Frwd",
                                    ifelse(ref_base == complement(`Allele1 - Forward`) | ref_base == complement(`Allele2 - Forward`), 
                                           "Complement",
                                            ifelse(ref_base == Snp1, 
                                                   case_when(Snp2==`Allele1 - Forward` | Snp2==`Allele2 - Forward` ~ "Frwd",
                                                             complement(Snp2)==`Allele1 - Forward` | complement(Snp2)==`Allele2 - Forward` ~ "Complement"),
                                                   ifelse(ref_base == complement(Snp1),
                                                          case_when(complement(Snp2)==`Allele1 - Forward` | complement(Snp2)==`Allele2 - Forward` ~ "Frwd",
                                                                    Snp2==`Allele1 - Forward` | Snp2==`Allele2 - Forward` ~ "Complement"),
                                                          ifelse(ref_base == Snp2,
                                                                 case_when(Snp1==`Allele1 - Forward` | Snp1==`Allele2 - Forward` ~ "Frwd",
                                                                           complement(Snp1) ==`Allele1 - Forward` | complement(Snp1)==`Allele2 - Forward` ~ "Complement"),
                                                                 ifelse(ref_base==complement(Snp2),
                                                                        case_when(complement(Snp1)==`Allele1 - Forward` | complement(Snp1)==`Allele2 - Forward` ~ "Frwd",
                                                                                  Snp1==`Allele1 - Forward` | Snp1==`Allele2 - Forward` ~ "Complement"
                                                                                  ),
                                                                        "SomethingisWrong"
                                                                        )
                                                                 )
                                                          )
                                                   )       
                                    )
                        )
                        ) %>% 
  mutate(TrueAllele1_Forward = ifelse(FrwdOrRevCalculation=="Frwd",`Allele1 - Forward`,complement(`Allele1 - Forward`)),
         TrueAllele2_Forward = ifelse(FrwdOrRevCalculation=="Frwd",`Allele2 - Forward`,complement(`Allele2 - Forward`))
         ) %>%
  select(`SNP Name`,`Sample ID`,TrueAllele1_Forward,TrueAllele2_Forward,
         ref_base,`Allele1 - AB`,`Allele2 - AB`,`GC Score`:`SNP Aux`,
         Chr:`CNV Confidence`,everything())-> mytidydata_results#)

setDT(mytidydata_results)
#fwrite(mytidydata_results, file = "mytidydata_final_results.txt", sep = "\t")
fwrite(mytidydata_results, file = "CEPH_chips_9-12_FinalReport.final_results.txt", sep = "\t", buffMB = 1024L)
