#!/usr/bin/env bash
set -euo pipefail

# InterProScan 版本（可在 EBI FTP 查看最新）
VERSION="5.72-104.0"

# 当前目录
CUR_DIR=$(pwd)
INSTALL_DIR="${CUR_DIR}/interproscan-${VERSION}"

# 下载地址（EBI FTP）
URL="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${VERSION}/interproscan-${VERSION}-64-bit.tar.gz"

echo "[INFO] Downloading InterProScan v${VERSION}..."
wget -nc "${URL}" -O "interproscan-${VERSION}.tar.gz"

echo "[INFO] Extracting..."
tar -xvzf "interproscan-${VERSION}.tar.gz" -C "${CUR_DIR}"

echo "[INFO] Installed to ${INSTALL_DIR}"
echo "[SUCCESS] Example usage:"
echo "  ${INSTALL_DIR}/interproscan.sh -i test_proteins.fasta -f tsv -o results.tsv"
