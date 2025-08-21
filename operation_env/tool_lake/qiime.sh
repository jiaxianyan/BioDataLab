#!/usr/bin/env bash
set -euo pipefail

# QIIME 2 版本
QIIME_VERSION=2024.10
PYTHON_VERSION=3.9

# 安装位置（当前目录下 qiime_env）
INSTALL_DIR="$(pwd)/qiime_env"

echo "[INFO] Step 1: Create conda env at ${INSTALL_DIR}"
conda create -y -p "${INSTALL_DIR}" python=${PYTHON_VERSION}

echo "[INFO] Step 2: Download QIIME 2 YAML"
wget -nc https://data.qiime2.org/distro/core/qiime2-${QIIME_VERSION}-py${PYTHON_VERSION//./}-linux-conda.yml -O qiime2.yml

echo "[INFO] Step 3: Install QIIME 2 into ${INSTALL_DIR}"
conda env update -p "${INSTALL_DIR}" --file qiime2.yml

echo "[INFO] Step 4: Test QIIME 2"
"${INSTALL_DIR}/bin/qiime" --help

echo "[SUCCESS] Installed QIIME 2 ${QIIME_VERSION}"
echo "Run it with:"
echo "  ${INSTALL_DIR}/bin/qiime --help"
