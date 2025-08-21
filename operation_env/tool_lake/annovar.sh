#!/usr/bin/env bash
set -euo pipefail

# 当前目录
CUR_DIR=$(pwd)
TAR_FILE="annovar.latest.tar.gz"
INSTALL_DIR="${CUR_DIR}/annovar"

if [[ ! -f "${TAR_FILE}" ]]; then
    echo "[ERROR] ${TAR_FILE} not found!"
    echo "Please download it manually from:"
    echo "  https://www.openbioinformatics.org/annovar/annovar_download_form.php"
    exit 1
fi

echo "[INFO] Extracting ANNOVAR..."
mkdir -p "${INSTALL_DIR}"
tar -xvzf "${TAR_FILE}" -C "${INSTALL_DIR}"

echo "[INFO] Installed ANNOVAR to ${INSTALL_DIR}"
echo "[SUCCESS] Example usage:"
echo "  perl ${INSTALL_DIR}/table_annovar.pl input.vcf ${INSTALL_DIR}/humandb/ -buildver hg19 -out output -remove -protocol refGene -operation g -nastring . -vcfinput"
