import pandas as pd
import sys

# user input
myfile = sys.argv[1]
#myfile1="CEPH_chips5-8_FinalReport.txt"
myfile1="CEPH_chips_9-12_FinalReport.txt"
# import data
#mydata = pd.read_csv('/home/data/LabData/SnpArrayData/Coriell/CEPH_chips_1_4/CEPH_chips1_4_full_data.table.txt.gz',sep='\t', skiprows=9, usecols=[1,18,19,0,10,11], dtype={'Sample ID':str,'Chr':str,'Position':int,'SNP Name':str,'Allele1 - Forward':str,'Allele2 - Forward':str} )
#mydata = pd.read_csv(myfile,sep='\t', skiprows=9, usecols=[1,18,19,0,10,11,4,20,*range(28,33),35,36], dtype={'SNP Name':str,\
mydata = pd.read_csv(myfile,sep='\t', skiprows=9, dtype={\
'SNP Name':str,\
'Sample ID':str,\
'Allele1 - Top':str,\
'Allele2 - Top':str,\
'GC Score':float,\
'Sample Name':str,\
'Sample Group':str,\
'Sample Index':str,\
'SNP Index':float,\
'SNP Aux':float,\
'Allele1 - Forward':str,\
'Allele2 - Forward':str,\
'Allele1 - Design':str,\
'Allele2 - Design':str,\
'Allele1 - AB':str,\
'Allele2 - AB':str,\
'Allele1 - Plus':str,\
'Allele2 - Plus':str,\
'Chr':str,\
'Position':int,\
'GT Score':float,\
'Cluster Sep':float,\
'SNP':str,\
'ILMN Strand':str,\
'Customer Strand':str,\
'Top Genomic Sequence':str,\
'Plus/Minus Strand':str,\
'Theta':float,\
'R':float,\
'X':float,\
'Y':float,\
'X Raw':float,\
'Y Raw':float,\
'B Allele Freq':float,\
'Log R Ratio':float,\
'CNV Value':float,\
'CNV Confidence':float} )

# checking type
type(mydata)

# for samtools faidx
mydata.loc[:,'region_samtools']= ("chr"+mydata['Chr']+":"+\
mydata['Position'].map(str)+"-"+\
mydata['Position'].map(str))

# filtering uninformative fields
mydata_filtered=mydata[(mydata.Chr !='XY') & (mydata.Chr!="0") & (mydata['Allele1 - Forward'] !="-") & (mydata.Chr !='MT')]

# Uncomment the following. Was commented out to avoid regenerating the same file again.
mydata_chrRegions= mydata_filtered['region_samtools']

# file for samtools faidx
mydata_chrRegions.to_csv(myfile1+'_regions_samtools.txt',index=False,header=False)

# output tidy file
mydata_filtered.to_csv(myfile1+'_tidy_data.txt',index=False,header=True)

# to August
#mydata_to_August = mydata_filtered[['Sample ID','Chr','Position','Allele1 - Forward','Allele2 - Forward','GC Score','GT Score','R','X','Y','X Raw','Y Raw','CNV Value','CNV Confidence']]
#mydata_to_August.to_csv(myfile1+'_tidy_data_to_August.tsv', index=False,header=True,sep="\t")


