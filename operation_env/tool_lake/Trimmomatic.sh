#!/usr/bin/env bash
set -euo pipefail

# Trimmomatic 版本
VERSION="0.39"

# 当前目录安装位置
CUR_DIR=$(pwd)
INSTALL_DIR="${CUR_DIR}/Trimmomatic"

# 下载链接（SourceForge）
URL="https://github.com/usadellab/Trimmomatic/releases/download/v${VERSION}/Trimmomatic-${VERSION}.zip"

echo "[INFO] Downloading Trimmomatic v${VERSION}..."
wget -nc "${URL}" -O "Trimmomatic-${VERSION}.zip"

echo "[INFO] Extracting..."
unzip -o "Trimmomatic-${VERSION}.zip" -d "${CUR_DIR}"

# 重命名目录为统一入口 Trimmomatic/
mv -f "${CUR_DIR}/Trimmomatic-${VERSION}" "${INSTALL_DIR}"

echo "[INFO] Installed Trimmomatic to ${INSTALL_DIR}"
echo "[SUCCESS] Example usage:"
echo "  java -jar ${INSTALL_DIR}/trimmomatic.jar PE -threads 4 input_R1.fq input_R2.fq out_R1.fq out_R1_unpaired.fq out_R2.fq out_R2_unpaired.fq ILLUMINACLIP:${INSTALL_DIR}/adapters/TruSeq3-PE.fa:2:30:10"
