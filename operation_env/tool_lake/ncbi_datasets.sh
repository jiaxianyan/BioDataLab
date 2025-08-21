#!/usr/bin/env bash
set -euo pipefail

# 版本（可修改为最新）
VERSION="16.26.1"

# 当前目录
CUR_DIR=$(pwd)
INSTALL_DIR="${CUR_DIR}/ncbi_datasets"

# 自动检测系统架构
OS=$(uname -s)
ARCH=$(uname -m)

if [[ "$OS" == "Linux" && "$ARCH" == "x86_64" ]]; then
    FILE="datasets-linux-amd64.zip"
elif [[ "$OS" == "Linux" && "$ARCH" == "aarch64" ]]; then
    FILE="datasets-linux-aarch64.zip"
elif [[ "$OS" == "Darwin" && "$ARCH" == "x86_64" ]]; then
    FILE="datasets-macos-amd64.zip"
elif [[ "$OS" == "Darwin" && "$ARCH" == "arm64" ]]; then
    FILE="datasets-macos-arm64.zip"
else
    echo "[ERROR] Unsupported OS/ARCH: ${OS}/${ARCH}"
    exit 1
fi

URL="https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/${VERSION}/${FILE}"

echo "[INFO] Downloading NCBI datasets CLI ${VERSION}..."
wget -nc "${URL}" -O "${FILE}"

echo "[INFO] Extracting..."
unzip -o "${FILE}" -d "${INSTALL_DIR}"

echo "[INFO] Installed to ${INSTALL_DIR}"
echo "[SUCCESS] Example usage:"
echo "  ${INSTALL_DIR}/bin/datasets --help"
echo "  ${INSTALL_DIR}/bin/dataformat --help"
