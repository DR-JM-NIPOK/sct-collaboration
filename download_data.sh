#!/bin/bash
# download_data.sh - Download public data products

set -e

DATA_DIR="data"
mkdir -p $DATA_DIR

echo "Downloading DESI-DR2 BAO data..."
wget -O $DATA_DIR/desi_dr2_bao.txt https://data.desi.lbl.gov/public/dr2/bao/bao_measurements.txt

echo "Downloading DES-Y6 3x2pt data..."
wget -O $DATA_DIR/des_y6_3x2pt.txt https://des.ncsa.illinois.edu/releases/y6/3x2pt/des_y6_3x2pt.txt

echo "Downloading Planck PR4 likelihood data..."
# Planck data is available via the Planck Legacy Archive
echo "Planck data: https://pla.esac.esa.int"

echo "Downloading HSC-Y3 data..."
wget -O $DATA_DIR/hsc_y3_data.txt https://hsc-release.mtk.nao.ac.jp/dataservice/y3/weak_lensing/hsc_y3.txt

echo "Downloading KiDS-DR5 data..."
wget -O $DATA_DIR/kids_dr5_data.txt https://kids.strw.leidenuniv.nl/DR5/data/kids_dr5.txt

echo "All data downloaded to $DATA_DIR/"