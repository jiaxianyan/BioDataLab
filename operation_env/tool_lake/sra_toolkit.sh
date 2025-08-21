#!/usr/bin/env bash
set -euo pipefail

# SRA Toolkit 版本（可改成最新）
SRA_VERSION=3.1.1

# 当前目录
CUR_DIR=$(pwd)
INSTALL_DIR="${CUR_DIR}/sra_toolkit"

# 自动检测系统架构
OS=$(uname -s)
ARCH=$(uname -m)

if [[ "$OS" == "Linux" && "$ARCH" == "x86_64" ]]; then
    FILE="sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz"
elif [[ "$OS" == "Linux" && "$ARCH" == "aarch64" ]]; then
    FILE="sratoolkit.${SRA_VERSION}-centos_linux64-aarch64.tar.gz"
else
    echo "[ERROR] Unsupported OS/ARCH: ${OS}/${ARCH}"
    exit 1
fi

URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/${FILE}"

echo "[INFO] Downloading SRA Toolkit ${SRA_VERSION}..."
wget -nc "${URL}" -O "${FILE}"

echo "[INFO] Extracting..."
tar -xzf "${FILE}"

# 重命名到 sra_toolkit（统一入口）
EXTRACTED_DIR=$(tar -tzf "${FILE}" | head -1 | cut -f1 -d"/")
mv -f "${EXTRACTED_DIR}" "${INSTALL_DIR}"

echo "[INFO] Installed SRA Toolkit to ${INSTALL_DIR}"
echo "[SUCCESS] Example usage:"
echo "  ${INSTALL_DIR}/bin/fasterq-dump --help"
