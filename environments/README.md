# building giggle and stix from polymerization environment

cd $GENEFUSION/environments

conda env create -f environment-polymerization.yml

conda activate polymerization

# build giggle
export HTS_INC=$CONDA_PREFIX/include
export HTS_LIB=$CONDA_PREFIX/lib

cd /data/jake/giggle
make clean
make

# build stix
cd /data/jake/stix-fork
make clean
make