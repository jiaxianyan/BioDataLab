#!/usr/bin/env bash
set -euo pipefail

# 当前目录
CUR_DIR=$(pwd)
wget http://annovar.openbio.org/annovar_20210630.tar.gz
tar -zxvf annovar_20210630.tar.gz
cd annovar_20210630

source annovar/table_annovar.pl -buildannovar -update -dir ./humandb