#!/usr/bin/env bash
set -euo pipefail

# 版本（可改最新）
VERSION="2.5.4"

# 当前目录
CUR_DIR=$(pwd)
INSTALL_DIR="${CUR_DIR}/bowtie2-${VERSION}"

# 下载链接（SourceForge 官方）
URL="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${VERSION}/bowtie2-${VERSION}-linux-x86_64.zip/download"

echo "[INFO] Downloading Bowtie2 v${VERSION}..."
wget -nc "${URL}" -O "bowtie2-${VERSION}.zip"

echo "[INFO] Extracting..."
unzip -o "bowtie2-${VERSION}.zip" -d "${CUR_DIR}"

# 添加执行权限
# chmod +x "${INSTALL_DIR}/bowtie2" "${INSTALL_DIR}/bowtie2-*"

echo "[INFO] Installed to ${INSTALL_DIR}"
echo "[SUCCESS] Example usage:"
echo "  ${INSTALL_DIR}/bowtie2 --help"
echo "  ${INSTALL_DIR}/bowtie2-build ref.fasta ref_index"
echo "  ${INSTALL_DIR}/bowtie2 -x ref_index -U reads.fq -S output.sam"
