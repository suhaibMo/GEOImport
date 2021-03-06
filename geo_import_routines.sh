#!/bin/bash

geo_experiment() {
   expAcc=$1
   echo $expAcc | grep 'GSE' > /dev/null
   return $?
}

get_ena_study_id(){
 expAcc=$1
 geo_ena_mapping=$2 
geo_experiment $expAcc
 if [ $? -ne 0 ]; then
  echo "ERROR: Not GEO experiment type - $expAcc" >&2
 fi
 
 # mapping index
 idx=`cat $geo_ena_mapping | grep -n "$expAcc" | head -n1 | awk -F":" '{print $1}'`
 cat $geo_ena_mapping | awk 'FNR == '$idx' {print $2}'
 }

elibility_error_code(){
expAcc=$1
geo_experiment $expAcc
 if [ $? -ne 0 ]; then
  echo "ERROR: Not GEO experiment type - $expAcc" >&2
 fi
GEO_IMPORT_FOLDER=/nfs/production3/ma/home/atlas3-production/GEO_import/GEOImportDownload 
cat $GEO_IMPORT_FOLDER/${expAcc}_output/${expAcc}_atlas_eligibility.out | grep "Atlas eligibility fail codes:" | awk -F":" '{print $2}' | sed 's/-//g'
}

geo_fixable(){
expAcc=$1
error_code=`elibility_error_code $expAcc`
ERROR_CODES=/nfs/production3/ma/home/atlas3-production/GEO_import/geo_import_supporting_files/error_codes.txt 
code=`echo $error_code | sed 's/,/\t/g' | cut -f 1 | sed -e 's/^\s*//'`
comment=`cat $ERROR_CODES | grep -w "$code" | awk -F"\t" '{print $3}'`
echo -e $comment
}

exp_meta_info(){ 
expAcc=$1
error_code=`elibility_error_code $expAcc`
GEO_IMPORT_FOLDER=/nfs/production3/ma/home/atlas3-production/GEO_import/GEOImportDownload/${expAcc}_output
    if [ -s $GEO_IMPORT_FOLDER/${expAcc}_family.soft ]; then
   title=`cat $GEO_IMPORT_FOLDER/${expAcc}_family.soft | grep "Series_title" | awk -F"=" '{print $2}'`
   organism=`cat $GEO_IMPORT_FOLDER/${expAcc}_family.soft | grep -P "Platform_organism" | awk -F"=" '{print $2}' | sort -u`
   no_of_samples=`cat $GEO_IMPORT_FOLDER/${expAcc}_family.soft | grep -c "Series_sample_id"`
   else 
  echo "ERROR: soft file doesn't exist - $expAcc" 
 fi  
 if [ -s $GEO_IMPORT_FOLDER/${expAcc}-sdrf.txt ]; then
      no_of_replicates=`cat $GEO_IMPORT_FOLDER/${expAcc}-sdrf.txt | cut -d$'\t' -f 1 | tail -n+2 |  uniq -c | awk -F" " '{print $1}' | sort -u`
      factor_value=`cat $GEO_IMPORT_FOLDER/${expAcc}-sdrf.txt | head -n1 | tr "\t" "\n" | grep -P "FactorValue" | awk -F" " '{print $2}' | sed -r 's/(\[|\])//g' | tr "\n" " "`
    else 
  echo "ERROR: sdrf file doesn't exist - $expAcc" 
 fi 
 
 echo -e $title'\t'$organism'\t'$no_of_samples'\t'$no_of_replicates'\t'$factor_value
}

rename_magetab(){
expAcc=$1
if [ -e ${expAcc}.idf.txt ]; then
  mv ${expAcc}.idf.txt $(echo ${expAcc}.idf.txt | sed -e "s/idf/merged.idf/")
fi
}

atlas_loaded_experiments(){
dbConnection=$1
# experiment that are already loaded
expLoaded=`echo -e "select accession from experiment where type like 'RNASEQ%';" \
     | psql $dbConnection | tail -n +3 | head -n -2 | grep -v "\E-MTAB-" | grep -v "\E-ERAD" | sed 's/E-GEOD-/GSE/g' | sed 's/^ //g' | sort -u`

# experiments ongoing in Atlas
inAtlas=`echo -e "select jobobject from atlas_jobs;" \
     | psql $dbConnection | tail -n +3 | head -n -2 | grep -v "\E-MTAB-" | grep -v "\E-ERAD" | grep -v "any" | sed 's/E-GEOD-/GSE/g' | sed 's/^ //g' | sed '/^$/d' | sort -u`

echo -e "$expLoaded $inAtlas" |  tr " " "\n"
}

exp_loading_check(){
expAcc=$1
geo_ena_mapping=$2
dbConnection=$3
count=`echo "select count(*) from atlas_eligibility where geo_acc='$expAcc';" | psql $dbConnection | tail -n +3 | head -n1 | sed 's/ //g'`
count_in_atlas=`atlas_loaded_experiments $dbConnection | grep -c "$expAcc"`

if [ $count_in_atlas -ne 0 ]; then
      echo "$expAcc - previously loaded in atlas" >&1
       return 1
elif [ $count -ne 0 ]; then
      echo "$expAcc - already exist in atlas eligbility table" >&1
       return 1
else 
      echo "loading $expAcc in db"
      load_eligibility_to_db $expAcc $geo_ena_mapping $dbConnection
       return 0
 fi
}
 

load_eligibility_to_db(){
  expAcc=$1
  geo_ena_mapping=$2
  dbConnection=$3
  ena_id=`get_ena_study_id $expAcc $geo_ena_mapping`
   if [ $? -ne 0 ]; then
  echo "ENA-ID doesn't match for $expAcc" >&2
   fi
  atlas_id=`echo -e $expAcc | sed 's/GSE/E-GEOD-/g'`
  error_code=`elibility_error_code $expAcc`
  comment=`geo_fixable $expAcc`
  title=`exp_meta_info $expAcc | awk -F "\t" '{print $1}'`
  organism=`exp_meta_info $expAcc | awk -F "\t" '{print $2}'`
  no_of_samples=`exp_meta_info $expAcc | awk -F "\t" '{print $3}'`
  no_of_replicate=`exp_meta_info $expAcc | awk -F "\t" '{print $4}'`
  factor_value=`exp_meta_info $expAcc | awk -F "\t" '{print $5}'`
  
  echo "insert into atlas_eligibility values (current_timestamp(0),'$atlas_id','$ena_id','$expAcc','$error_code', \
      '$comment','$title','$organism','$no_of_samples','$no_of_replicate','$factor_value');" | psql $dbConnection
}


# query experiments for species that are having reference genome in ISL
# extract species and subspecies from RNA-seqer API call for ensembl, plants metazoa and wbps. 
query_atlas_eligibility() {
IFS="
"
dbConnection=$1

species_in_ensembl=`curl -s "https://www.ebi.ac.uk/fg/rnaseq/api/tsv/0/getOrganisms/ensembl" | cut -f 1 | tail -n +2 | awk '{ gsub("_", " ") ; print $0 }' | sed -e 's/^./\U&/'`
species_in_plants=`curl -s "https://www.ebi.ac.uk/fg/rnaseq/api/tsv/0/getOrganisms/plants" | cut -f 1 | tail -n +2 | awk '{ gsub("_", " ") ; print $0 }' | sed -e 's/^./\U&/'`
species_in_metazoa=`curl -s "https://www.ebi.ac.uk/fg/rnaseq/api/tsv/0/getOrganisms/metazoa" | cut -f 1 | tail -n +2 | awk '{ gsub("_", " ") ; print $0 }' | sed -e 's/^./\U&/'`
species_in_wbps=`curl -s "https://www.ebi.ac.uk/fg/rnaseq/api/tsv/0/getOrganisms/wbps" | cut -f 1 | tail -n +2 | awk '{ gsub("_", " ") ; print $0 }' | sed -e 's/^./\U&/'`

species_in_isl=`printf "%s\n" $species_in_ensembl $species_in_plants $species_in_metazoa $species_in_wbps`
output=""
echo -e "ae2_acc\tena_study_id\tgeo_acc\torganism\tstatus\texp_type"
for species in $species_in_isl; do
     var=`echo "select DISTINCT ae2_acc, ena_study_id, geo_acc, organism,status, exp_type from atlas_eligibility where organism like'%$species';" | psql $dbConnection | tail -n +3 | head -n -2`
     if [ $? -ne 0 ]; then
      echo "query didn't succeed for the $species"
     fi 
     output+="$var\n"
 done

echo -e "$output" | sed -e '/^ *$/d'
}