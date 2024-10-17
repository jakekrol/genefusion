#!/usr/bin/env bash

echo "test"
echo exit

k=${1:-1000}

complain() {
    echo $1
    exit 1
}

# get query bed
bed=$(ls beds/ | head -n 1)
matching=false
# find a line with matching chr
for ((i=1; i<=k; i++)); do
    data=$(zcat beds/$bed | head -n $k | sed -n "${i}p")
    arr=($data)
    l_chr=${arr[0]}
    l_strand=${arr[3]}
    r_chr=${arr[4]}
    r_strand=${arr[7]}
    if [ "$l_chr" == "$r_chr" ] && [ "$l_strand" == "$r_strand" ]; then
        matching=true
        break
    fi
done
if [ "$matching" == false ];
then
    complain "No matching chromosome&strand interval found in the first $k lines of $bed"
fi
left="${l_chr}:${arr[1]}-${arr[2]}"
right="${r_chr}:${arr[5]}-${arr[6]}"

stix -i index -d stix.ped.db -t DEL -l $left -r $right -s 500 > out.stix
