# simple way to giggle index many shards

do_giggle() {
    giggle index -i "beds/*.gz" -o index -f -s
}

export -f do_giggle

nohup bash -c 'do_giggle' > giggle_indexer.log 2>&1 &

#stix -i index -p tumour.ped -d stix.ped.db
