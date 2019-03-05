#!/bin/bash

NLDAS_RAW_LOCATION="NLDASRAW"
MODIS_RAW_LOCATION="MODISRAW"

NLDAS_TEXT_LOCATION="NLDASTEXT"
MODIS_TEXT_LOCATION="MODISTEXT"


while [ "$1" != "" ]; do
  case $1 in 
    -s | --start )        shift
                          START=$1
                          ;;
    -e | --end )          END=$2
                          ;;
  esac
  shift
done

touch filelistNLDAS.txt
touch filelistMODIS.txt

R CMD BATCH --vanilla "--args start=\"$START\" end=\"$END\"" createfilelist.R

rm createfilelist.Rout

wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies -nc -i filelistNLDAS.txt -P "$NLDAS_RAW_LOCATION"
wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies -nc -i filelistMODIS.txt -P "$MODIS_RAW_LOCATION"

rm filelistNLDAS.txt
rm filelistMODIS.txt

for filename in $NLDAS_RAW_LOCATION/*.grb; do
  wgrib -s "$filename" | grep ":TMP:" | wgrib -i -text "$filename" -o "$NLDAS_TEXT_LOCATION/$(basename "$filename" .grb).txt"
done

for filename in $MODIS_RAW_LOCATION/*.hdf; do
  gdal_translate -of XYZ HDF4_EOS:EOS_GRID:\""$filename"\":MODIS_Grid_8Day_1km_LST:LST_Day_1km "$MODIS_TEXT_LOCATION/$(basename "$filename" .hdf).txt"
done
