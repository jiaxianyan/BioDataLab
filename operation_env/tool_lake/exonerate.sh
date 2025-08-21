#!/usr/bin/env bash
set -euo pipefail

# Exonerate 版本（稳定版）
VERSION="2.4.0"

# 当前目录
CUR_DIR=$(pwd)
INSTALL_DIR="${CUR_DIR}/exonerate-${VERSION}"

# 下载地址（EBI FTP）
URL="https://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-${VERSION}-x86_64.tar.gz"

echo "[INFO] Downloading Exonerate v${VERSION}..."
wget -nc "${URL}" -O "exonerate-${VERSION}.tar.gz"

echo "[INFO] Extracting..."
tar -xvzf "exonerate-${VERSION}.tar.gz" -C "${CUR_DIR}"

# 添加执行权限
chmod +x "${INSTALL_DIR}/bin/exonerate"

echo "[INFO] Installed to ${INSTALL_DIR}"
echo "[SUCCESS] Example usage:"
echo "  ${INSTALL_DIR}/bin/exonerate --help"
echo "  ${INSTALL_DIR}/bin/exonerate query.fasta target.fasta"
