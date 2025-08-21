#!/usr/bin/env bash
set -euo pipefail

# 版本（可改成最新发布版本）
VERSION="2.3"

# 当前目录安装位置
CUR_DIR=$(pwd)
INSTALL_DIR="${CUR_DIR}/mash"

# 自动检测系统架构
OS=$(uname -s)
ARCH=$(uname -m)

if [[ "$OS" == "Linux" && "$ARCH" == "x86_64" ]]; then
    FILE="mash-Linux64-v${VERSION}.tar"
elif [[ "$OS" == "Linux" && "$ARCH" == "aarch64" ]]; then
    echo "[ERROR] Mash 官方没有预编译 aarch64 版本，请源码编译"
    exit 1
else
    echo "[ERROR] Unsupported OS/ARCH: ${OS}/${ARCH}"
    exit 1
fi

URL="https://github.com/marbl/Mash/releases/download/v${VERSION}/${FILE}"

echo "[INFO] Downloading Mash v${VERSION}..."
wget -nc "${URL}" -O "${FILE}"

echo "[INFO] Extracting..."
tar -xf "${FILE}"

# 解压出来的就是 bin，可重命名为 mash/
EXTRACTED_DIR="mash-Linux64-v${VERSION}"
mv -f "${EXTRACTED_DIR}" "${INSTALL_DIR}"

echo "[INFO] Installed Mash to ${INSTALL_DIR}"
echo "[SUCCESS] Example usage:"
echo "  ${INSTALL_DIR}/mash sketch -o ref.msh ref.fna"
echo "  ${INSTALL_DIR}/mash dist ref.msh query.fna"
