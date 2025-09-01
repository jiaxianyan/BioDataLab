# Nucleotide Sequence Databases, ASMdb database
conda install -y -c bioconda fastp
conda install -y -c bioconda samtools
conda install -y -c bioconda sra-tools
conda install -y -c bioconda BatMeth2

pip install GEOparse

# # install zlib
# ./configure --prefix=/disk1/glli/tools/zlib-1.2.11/
# make
# make install
# export C_INCLUDE_PATH=$C_INCLUDE_PATH:/disk1/glli/tools/zlib-1.2.11/include
# export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/disk1/glli/tools/zlib-1.2.11/include
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH::/disk1/glli/tools/zlib-1.2.11/lib
# export LIBRARY_PATH=$LIBRARY_PATH::/disk1/glli/tools/zlib-1.2.11/lib

# # install gsl
# ./configure --prefix=/disk1/glli/tools/gsl-2.4/
# make
# make install
# export C_INCLUDE_PATH=$C_INCLUDE_PATH:~/software/gsl-2.4/include
# export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:~/software/gsl-2.4/include
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/software/gsl-2.4/lib
# export LIBRARY_PATH=$LIBRARY_PATH:~/software/gsl-2.4/lib

# git clone https://github.com/GuoliangLi-HZAU/BatMeth2.git
# cd BatMeth2
# ./configure
# make
# make install