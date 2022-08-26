#!/bin/bash
DATABASE_SCHEMA=sera_gis
OUTPUT_FOLDER=<OUTPUT_LOCATION>
USER=<USER>
PASSWORD=<PASSWORD>

listvar="CensusCounties_2010_Polygons_HighResolution CensusDivisions_2010_Polygons_HighResolution CensusStates_2010_Polygons_HighResolution CensusUrbanAreaCodes_2010_Polygons_HighResolution NERC_2014 NERC_And_California"

for f in $listvar
do
    pgsql2shp -f $OUTPUT_FOLDER$f -h 1lp11dnpdb1.nrel.gov -u $USER -P $PASSWORD sera sera_gis.$f
done