#!/bin/bash
index=0
path=${PWD}
DIRNAME=$(basename "$(dirname "$path")")/$(basename "$path")
for filename in silicon.*.griff; do
    ess_silicon_ana_bic $filename

    input="Elastic.txt"
    VAR=""
    while read line; do
	      VAR+="$line "
    done < $input
    #echo $VAR
    ess_griffanautils_extractevts $filename Elastic.$index.griff $VAR
    ess_silicon_sim -n15000000 -j4 event_gen="GriffGen" griff_file="Elastic.$index.griff" -l "ESS_SiliconArticle" -o="../../../ScreenedNuclearRecoil/$DIRNAME/Elastic/ScreenedNuclearRecoil.$index.griff"

    input="Inelastic.txt"
    VAR=""
    while read line; do
	      VAR+="$line "
    done < $input
    #echo $VAR
    ess_griffanautils_extractevts $filename Inelastic.$index.griff $VAR
    ess_silicon_sim -n1500000 -j4 event_gen="GriffGen" griff_file="Inelastic.$index.griff" -l "ESS_SiliconArticle" -o="../../../ScreenedNuclearRecoil/$DIRNAME/Inelastic/ScreenedNuclearRecoil.$index.griff"

    input="Coulomb.txt"
    VAR=""
    while read line; do
	      VAR+="$line "
    done < $input
    # echo $VAR
    ess_griffanautils_extractevts $filename Coulomb.$index.griff $VAR
    ess_silicon_sim -n1500000 -j4 event_gen="GriffGen" griff_file="Coulomb.$index.griff" -l "ESS_SiliconArticle" -o="../../../ScreenedNuclearRecoil/$DIRNAME/Coulomb/ScreenedNuclearRecoil.$index.griff"

    let index+=1

done


