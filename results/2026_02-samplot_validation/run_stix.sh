#!/usr/bin/env bash
INDEX=/data/jake/stix-pcawg-dna
DIR_SHARDFILES="$INDEX/shardfiles"
PROCS=30
OUTDIR="$(pwd)/stix_results"
OUTDIR_LOG="$(pwd)/stix_gargs_logs"
mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR_LOG"
for tissue in kidney liver ovary; do
    # find shardfile for tissue
    shardfile=$(ls "${DIR_SHARDFILES}" | grep -i "${tissue}" | grep tumor)
    shardfile="${DIR_SHARDFILES}/${shardfile}"
    shardfile=$(realpath "${shardfile}")
    if [ -z "${shardfile}" ]; then
        echo "Error: shardfile not found for tissue: ${tissue}"
        exit 1
    fi
    # queries
    stix_input=$(ls | grep stix_coords.tsv$ | grep ${tissue})
    stix_input=$(realpath "${stix_input}")
    if [ -z "${stix_input}" ]; then
        echo "Error: stix input file not found for tissue: ${tissue}"
        exit 1
    fi

    # out del
    outfile_del="${OUTDIR}/${tissue}_del.stix"
    outfile_del=$(realpath -m "${outfile_del}")

    # out inv
    outfile_inv="${OUTDIR}/${tissue}_inversions.stix"
    outfile_inv=$(realpath -m "${outfile_inv}")

    # query both del and inv
    for outfile in "${outfile_del}" "${outfile_inv}"; do
        logfile="${OUTDIR_LOG}/gargs_${tissue}_$(basename "${outfile}").log"
        if [[ "${outfile}" == *"_del.stix" ]]; then
            type="DEL"
        else
            type="INV"
        fi
        (
            cd $INDEX || { echo "Error: failed to cd to $INDEX"; exit 1; }
            gargs --log=$logfile -p $PROCS \
                "stix -t $type -l {1}:{2}-{3} -r {5}:{6}-{7} -B ${shardfile} > ${outfile}.{0}--{4}.out" \
                < "${stix_input}"
        )
    done
done
        
