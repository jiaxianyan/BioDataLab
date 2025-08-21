#!/usr/bin/env bash
set -euo pipefail

# 版本号（可以在官网或 GitHub Releases 查看最新）
VERSION="2.0"

# 当前目录
CUR_DIR=$(pwd)
INSTALL_DIR="${CUR_DIR}/arts-${VERSION}"

# 下载地址（GitHub Release，假设有 jar 包或 zip 包）
URL="https://github.com/miguelinux314/arts/releases/download/v${VERSION}/arts-${VERSION}.zip"

echo "[INFO] Downloading ARTS v${VERSION}..."
wget -nc "${URL}" -O "arts-${VERSION}.zip"

echo "[INFO] Extracting..."
unzip -o "arts-${VERSION}.zip" -d "${CUR_DIR}"

# 赋予执行权限（假设有启动脚本或 jar）
chmod +x "${INSTALL_DIR}/arts.sh" || true

echo "[INFO] Installed ARTS to ${INSTALL_DIR}"
echo "[SUCCESS] Example usage:"
echo "  java -jar ${INSTALL_DIR}/arts.jar -i input.fasta -o results/"
