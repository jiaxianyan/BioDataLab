#!/usr/bin/env bash
set -euo pipefail

# 下载链接
TARBALL_URL="https://www.niehs.nih.gov/sites/default/files/2024-02/artbinmountrainier2016.06.05linux64.tgz"
TARBALL_NAME=$(basename "$TARBALL_URL")

echo "[INFO] Installation directory: $(pwd)"

# 下载（如果文件不存在）
if [ ! -f "$TARBALL_NAME" ]; then
    echo "[INFO] Downloading $TARBALL_NAME ..."
    curl -L -o "$TARBALL_NAME" "$TARBALL_URL"
else
    echo "[INFO] $TARBALL_NAME already exists, skipping download."
fi

# 解压
echo "[INFO] Extracting $TARBALL_NAME ..."
tar xzvf "$TARBALL_NAME"

echo "[INFO] Installation complete."
echo "Executables are available under: $(pwd)/art_bin_MountRainier"