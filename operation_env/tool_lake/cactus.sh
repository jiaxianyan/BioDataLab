#!/usr/bin/env bash
set -euo pipefail

# Cactus 版本
VERSION="2.7.0"

# 当前目录安装位置
CUR_DIR=$(pwd)
INSTALL_DIR="${CUR_DIR}/cactus_env"

echo "[INFO] Step 1: Create conda environment at ${INSTALL_DIR}"
conda create -y -p "${INSTALL_DIR}" -c conda-forge -c bioconda cactus=${VERSION}

echo "[INFO] Step 2: Test installation"
"${INSTALL_DIR}/bin/cactus-graphmap" --help || true
"${INSTALL_DIR}/bin/cactus-align" --help || true

echo "[SUCCESS] Cactus v${VERSION} installed in ${INSTALL_DIR}"
echo "Run it with:"
echo "  ${INSTALL_DIR}/bin/cactus jobStore input.txt output.hal"
