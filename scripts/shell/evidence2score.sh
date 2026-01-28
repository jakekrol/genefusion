#!/usr/bin/env bash
# defaults
combine_modal="/data/jake/genefusion/scripts/shell/combine_modal.sh"
score_script="/data/jake/genefusion/scripts/python/score_fusions.py"
score_yaml="/data/jake/genefusion/config/score_dna1.0_t0.5_r0.5_u50.yaml"
while [[ $# -gt 0 ]]; do
  case $1 in
    -t|--tumor) tumor=$2; shift;; # read and sample tsv for tumor
    -n|--normal) normal=$2; shift;; # read and sample tsv for normal
    -c|--combine_modal) combine_modal=$2; shift;;
    -s|--score_script) score_script=$2; shift;; 
    -y|--score_yaml) score_yaml=$2; shift;; # yaml template for scoring
    -nt|--ntumor) ntumor=$2; shift;; # num. tumor samples
    -nn|--nnormal) nnormal=$2; shift;; # num. normal samples
    -o|--outdir) outdir=$2; shift;; # output directory
  esac
  shift
done

# validation
if [ -z "${tumor:-}" ] || [ -z "${normal:-}" ] || [ -z "${ntumor:-}" ] || [ -z "${nnormal:-}" ] || [ -z "${outdir:-}" ]; then
  echo "Usage: $0 -t <tumor read and sample tsv> -n <normal read and sample tsv> -nt <num tumor samples> -nn <num normal samples> -o <output dir> [-c <combine_modal.sh>] [-s <score_fusions.py>] [-y <score yaml>]"
  exit 1
fi
if [ ! -f "${tumor}" ]; then
  echo "Tumor read and sample tsv not found: ${tumor}"
  exit 1
fi
if [ ! -f "${normal}" ]; then
    echo "Normal read and sample tsv not found: ${normal}"
    exit 1
fi
if [ ! -f "${combine_modal}" ]; then
    echo "Combine modal script not found: ${combine_modal}"
    exit 1
fi
if [ ! -f "${score_script}" ]; then
    echo "Score script not found: ${score_script}"
    exit 1
fi
if [ ! -f "${score_yaml}" ]; then
    echo "Score yaml not found: ${score_yaml}"
    exit 1
fi

mkdir -p $outdir && cd $outdir

printf "${tumor}\n${normal}\n" > combine_modal.input
# symlink to score runner
ln -s ${score_script} ./score_fusions.py
ln -s ${combine_modal} ./combine_modal.sh
# get combine all modalities and onekg
./combine_modal.sh \
    -i combine_modal.input \
    -o $(pwd) \
    -t /data/jake/tmp
cp ${score_yaml} ./score_dna1.0_t0.5_r0.5_u50.yaml 
sed -i "s|^pop_size_dna_normal.*|pop_size_dna_normal: ${nnormal}|" score_dna1.0_t0.5_r0.5_u50.yaml
sed -i "s|^pop_size_dna_tumor.*|pop_size_dna_tumor: ${ntumor}|" score_dna1.0_t0.5_r0.5_u50.yaml
# score
./score_fusions.py --input modal_onekg_combined_read_sample_comb_fill_no_dups.tsv \
    --output scored_dna1.0_t0.5_r0.5_u50.tsv \
    --score_yaml score_dna1.0_t0.5_r0.5_u50.yaml 