#!/bin/sh
# Please change #Spacer sequence without PAM!
read -p 'Input gRNA sequence (DO NOT INCLUDE PAM):' Spacer

if [ ${#Spacer} = 20 ];then
  echo "
  This is OfferOT v1.0 created by RN
  Spacer:$Spacer is imported successfully
  "
else
  echo "
  Error: Input needs to be 20 nt!
  "
  exit
fi

read -p 'Choose the length of seed sequence (8 or 12 nt):' VAR

# Importng
SandP=`echo $Spacer | sed -e 's/$/NGG/'`
SP=`echo $SandP`

if [ "$VAR" = 12 ]; then
  Seed=`echo $SP | cut -c 9-23`
elif [ "$VAR" = 8 ]; then
  Seed=`echo $SP | cut -c 13-23`
else
  echo "
  Error: Input needs to be 8 or 12!
  "
  exit
fi

#plus strand
# Data download
wget -q https://gggenome.dbcls.jp/ja/hg38/3/+/nogap/$SP.csv -O plus_1.csv
wget -q https://gggenome.dbcls.jp/ja/hg38/1/+/nogap/$Seed.csv -O plus_2.csv

# File import
OT1=`ls plus_1.csv`
OT2=`ls plus_2.csv`
echo "$OT1 was imported from GGGenome.
$OT2 was imported from GGGenome.
File processing started
"

# File process (OT1)
cut -f 9 -d "," ${OT1} > OT1_plus.csv
sed '1,5d' OT1_plus.csv | sed 's/"//g' > OT1_processed_plus.csv

# File process (OT1)
cut -f 9 -d "," ${OT2} > OT2_plus.csv
sed '1,5d' OT2_plus.csv | sed 's/"//g' > OT2_processed_plus.csv

grep -f OT2_processed_plus.csv OT1_processed_plus.csv > OT_alignment_plus.csv

grep "AG$" OT_alignment_plus.csv > AG.csv
grep "GG$" OT_alignment_plus.csv > GG.csv
cat AG.csv GG.csv > plus.csv

#plus strand
# Data download
wget -q https://gggenome.dbcls.jp/ja/hg38/3/-/nogap/$SP.csv -O minus_1.csv
wget -q https://gggenome.dbcls.jp/ja/hg38/1/-/nogap/$Seed.csv -O minus_2.csv

# File import
OT1m=`ls minus_1.csv`
OT2m=`ls minus_2.csv`
echo "$OT1m was imported from GGGenome.
$OT2m was imported from GGGenome.
"

# File process (OT1)
cut -f 9 -d "," ${OT1m} > OT1_minus.csv
sed '1,5d' OT1_minus.csv | sed 's/"//g' > OT1_processed_minus.csv

# File process (OT1)
cut -f 9 -d "," ${OT2m} > OT2_minus.csv
sed '1,5d' OT2_minus.csv | sed 's/"//g' > OT2_processed_minus.csv

grep -f OT2_processed_minus.csv OT1_processed_minus.csv > OT_alignment_minus.csv

grep "^CC" OT_alignment_minus.csv > CC.csv
grep "^CT" OT_alignment_minus.csv > CT.csv
cat CC.csv CT.csv > minus.csv

# Comparison step
echo "The candidates of OT was imported successfully."
listp=`cat plus.csv`
listm=`cat minus.csv`

echo "
OTs were detected successfully.
==============
$listp
$listm
=============="

OT_list=`cat plus.csv minus.csv`
grep "$OT_list" $OT1 > OT_list_p.csv
grep "$OT_list" $OT1m > OT_list_m.csv
cat OT_list_p.csv OT_list_m.csv > OT_list_final.csv

awk -F ',' '{print $9 "," $1 ":" $6 "-" $7}' OT_list_final.csv | sed -e 's/"//g' > UCSC_list_final.csv
awk -F ',' '{print $1"\t"$3"\t"$4}' OT_list_final.csv | sed -e 's/"//g' > OT_candidate.bed

mkdir -p analysis/OT
mkdir -p analysis/intermediate

mv OT_list_final.csv analysis/OT/
mv UCSC_list_final.csv analysis/OT/
mv *.csv analysis/intermediate
mv OT_candidate.bed analysis/OT/
