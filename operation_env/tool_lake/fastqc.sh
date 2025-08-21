#!/usr/bin/env bash
set -euo pipefail

# FastQC 版本（可改最新）
VERSION="0.12.1"

# 当前目录
CUR_DIR=$(pwd)
INSTALL_DIR="${CUR_DIR}/FastQC"

# 下载地址（SourceForge 官方）
URL="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${VERSION}.zip"

echo "[INFO] Downloading FastQC v${VERSION}..."
wget -nc "${URL}" -O "fastqc_v${VERSION}.zip"

echo "[INFO] Extracting..."
unzip -o "fastqc_v${VERSION}.zip" -d "${CUR_DIR}"

# 给 fastqc 主程序加执行权限
chmod +x "${INSTALL_DIR}/fastqc"

echo "[INFO] Installed to ${INSTALL_DIR}"
echo "[SUCCESS] Example usage:"
echo "  ${INSTALL_DIR}/fastqc test.fastq.gz"
echo "Results will be written in the current directory (HTML + .zip)."
