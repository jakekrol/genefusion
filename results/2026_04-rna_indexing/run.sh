#!/usr/bin/env bash

PROCS=8

tissues=(blood kidney liver ovary pancreas)

### split tumor normal
for t in "${tissues[@]}"; do
    echo "# splitting tumor normal for ${t}"
    $GENEFUSION/scripts/python/split_tumor_normal.py \
        --indir /data/jake/stix-pcawg-rna/${t}_all_excord_bed \
        --outdir_normal ${t}_normal \
        --outdir_tumor ${t}_tumor 
done
# many tissues have no normal samples
# check if normal dir is empty, if so remove it
for t in "${tissues[@]}"; do
    if [ -d "${t}_normal" ] && [ -z "$(ls -A ${t}_normal)" ]; then
        echo "Directory ${t}_normal is empty. Removing it."
        rm -ir "${t}_normal"
    fi
done

### shard data
for t in "${tissues[@]}"; do
    dnormal="${t}_normal"
    dtumor="${t}_tumor"
    if [ -d $dnormal ]; then
        echo "# sharding $dnormal data"
        $GENEFUSION/scripts/python/giggle_shard.py \
            --input $dnormal \
            --output $(pwd) \
            --prefix ${t}_normal.data. \
            --pattern ".*\.bed\.gz$" \
            --n_shards 4 \
            --shard_meta ${t}_normal_shard_metadata.yaml
    else
        echo "# $dnormal does not exist skipping"
    fi
    if [ -d $dtumor ]; then
        echo "# sharding $dtumor data"
        $GENEFUSION/scripts/python/giggle_shard.py \
            --input $dtumor \
            --output $(pwd) \
            --prefix ${t}_tumor.data. \
            --pattern ".*\.bed\.gz$" \
            --n_shards 4 \
            --shard_meta ${t}_tumor_shard_metadata.yaml
    fi
done


### giggle index

echo "# giggle indexing shards"
# get data dirs
ls | grep data | grep -E "tumor|normal" | grep -v yaml > shard_data_dirs.txt
mapfile -t shards < shard_data_dirs.txt
rm index.input || echo "index.input does not exist, creating new one"
for line in "${shards[@]}"; do
    # get tissue and specimen
    prefix_tissue_specimen="$(cut -f 1 -d "." <(echo $line))"
    # # get shard number
    suffix_shard_num="$(cut -f 3 -d "." <(echo $line))"
    tissue="$(cut -f 1 -d "_" <(echo $prefix_tissue_specimen))"
    specimen="$(cut -f 2 -d "_" <(echo $prefix_tissue_specimen))"
    # need quotes for expansion of bed files to be treated as single arg
    printf "${line}\t${prefix_tissue_specimen}\t${suffix_shard_num}\n" >> index.input
done

cat index.input | \
    gargs -v --log=index.gargs.log -p "${PROCS}" \
    "giggle index -s -i \"{0}/*.gz\" -o {1}.index.{2}" 2>&1 | tee index_all.log


echo "# creating ped files for shards"
mapfile -t shards < shard_data_dirs.txt
for shard in "${shards[@]}"; do
    prefix="$(echo $shard | cut -f 1 -d '.')"
    shard_num="$(echo $shard | cut -f 3 -d '.')"
    echo "# creating ped for shard ${shard}"
    $GENEFUSION/scripts/python/data2ped.py \
        --indir "${shard}" \
        --outfile "${prefix}.${shard_num}.ped"
done

echo "# stix indexing shards"
# get giggle indices
mapfile -t indices < <(ls | grep index | grep -v "log$" | grep -v "input$")
for idx in "${indices[@]}"; do
   shard_name="$(echo $idx | cut -f 1 -d '.')"
   shard_num="$(echo $idx | cut -f 3 -d '.')"
   ped="${shard_name}.${shard_num}.ped"
   ped_db="${shard_name}.${shard_num}.ped.db"
   echo "# stix indexing $idx , $ped , to $ped_db"
   stix -i "${idx}" -p "${ped}" -d "${ped_db}" -c 5
done

echo "# making shardfiles"
./make_shardfile.py \
    --indir . \
    --outdir shardfiles \
    --tissues blood,kidney,liver,ovary,pancreas \
    --specimens normal,tumor \
    --index_pattern "*.index.*" \
    --db_pattern "*.ped.db"
    