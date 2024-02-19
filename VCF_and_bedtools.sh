###################################################################

##vcftools##

#initialise bash shell
conda init bash

#create and activate the vcf environment
conda create --name vcf bioconda::vcftools
conda activate vcf

#change into the relevent directory
cd /Users/cats/Downloads/LIFE4141_coursework_resources-20240122/

#test different window values to see when uniquely high values start to appear#

###################################################################

##5000##

#run code for window-size = 5000 and window-step = 5000
vcftools --gzvcf LAB_NEN_ODN.clean_BI.ann.3mbChr5.vcf.gz --max-missing 0.8 --maf 0.05 --weir-fst-pop LAB_population.txt --weir-fst-pop NEN_population.txt --fst-window-size 5000 --fst-window-step 5000

#sort by highest fst value
sort -nrk 5 windowsize5000 | > windowsize5000_sort

#extract top 20 rows in file
sed -n 1,20p windowsize5000_sort | > 5000_top20

#used to view the data
less 5000_top20

###################################################################

##2500##

#run code for window-size = 2500 and window-step = 2500
vcftools --gzvcf LAB_NEN_ODN.clean_BI.ann.3mbChr5.vcf.gz --max-missing 0.8 --maf 0.05 --weir-fst-pop LAB_population.txt --weir-fst-pop NEN_population.txt --fst-window-size 2500 --fst-window-step 2500

#sort by highest fst value
sort -nrk 5 windowsize2500_2 | > windowsize2500_sort

#extract top 20 rows in file
sed -n 1,20p windowsize2500_sort | > 2500_top20

#used to view the data
less 2500_top20

###################################################################

##1000##

#run code for window-size = 1000 and window-step = 1000
vcftools --gzvcf LAB_NEN_ODN.clean_BI.ann.3mbChr5.vcf.gz --max-missing 0.8 --maf 0.05 --weir-fst-pop LAB_population.txt --weir-fst-pop NEN_population.txt --fst-window-size 1000 --fst-window-step 1000

#sort by highest fst value
sort -nrk 5 windowsize1000_2 | > windowsize1000_sort

#extract top 20 rows in file
sed -n 1,20p windowsize1000_sort | > 1000_top20

#used to view the data
less 1000_top20

###################################################################

##500##

#run code for window-size = 500 and window-step = 500
vcftools --gzvcf LAB_NEN_ODN.clean_BI.ann.3mbChr5.vcf.gz --max-missing 0.8 --maf 0.05 --weir-fst-pop LAB_population.txt --weir-fst-pop NEN_population.txt --fst-window-size 500 --fst-window-step 500

#sort by highest fst value
sort -nrk 5 windowsize500_2 | > windowsize500_sort

#extract top 20 rows in file
sed -n 1,20p windowsize500_sort | > 500_top20

#used to view the data
less 500_top20

###################################################################

#decided on 500 as the window value#

###################################################################

##bedtools##

#create env and install bedtools into it
conda create -y --name bedtools bioconda::bedtools

#activate bedtools environment
conda activate bedtools

#use the bedtools intersect function
bedtools intersect -wao -a 500_top20 -b C_excelsa_V5_braker2_wRseq.gff3 | > 500bedtools

#make ID column easier to read by isolating it with the scaf column
awk '{print $1, $15}' 500bedtools | > 500bedtoolsnew








