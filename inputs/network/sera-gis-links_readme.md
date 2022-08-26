SERA GIS shapefiles
===================
Date: July 2017

See <https://github.nrel.gov/SERA/sera-gis> repo for raw data, readme files, and scripts used to process the geographic tables stored in the sera-gis schema.


See [ProcessedShapefiles](<https://github.nrel.gov/SERA/sera-gis/tree/master/ProcessedShapefiles>) folder in sera-gis repo for final shapefiles that are of direct interest in the sera model.

Projection: NAD 83/Conus Albers

EPSG: 5070

*   Census Counties - [CensusCounties_2010_Polygons_HighResolution](<https://github.nrel.gov/SERA/sera-gis/tree/master/ProcessedShapefiles/CensusCounties_2010_Polygons_HighResolution.shp>)
*   Census Divisions - [CensusDivisions_2010_Polygons_HighResolution](<https://github.nrel.gov/SERA/sera-gis/tree/master/ProcessedShapefiles/CensusDivisions_2010_Polygons_HighResolution.shp>)
*   Census States - [CensusStates_2010_Polygons_HighResolution](<https://github.nrel.gov/SERA/sera-gis/tree/master/ProcessedShapefiles/CensusStates_2010_Polygons_HighResolution.shp>)
*   Census Urban Areas - [CensusUrbanAreaCodes_2010_Polygons_HighResolution](<https://github.nrel.gov/SERA/sera-gis/tree/master/ProcessedShapefiles/CensusUrbanAreaCodes_2010_Polygons_HighResolution.shp>)
*   Nerc Regions (2014) - [NERC_2014](<https://github.nrel.gov/SERA/sera-gis/tree/master/ProcessedShapefiles/NERC_2014.shp>)
*   Nerc Regions and California - [NERC_And_California](<https://github.nrel.gov/SERA/sera-gis/tree/master/ProcessedShapefiles/NERC_And_California.shp>)
*   Balancing Areas from REEDS - [BAs_2014_Polygons](<https://github.nrel.gov/SERA/sera-gis/tree/master/ProcessedShapefiles/BAs_2014_Polygons.shp>)	
	

The sera-gis schema contains additional geographic tables not listed in the ProcessedShapefiles folder.  These tables can be exported in bulk to shapefile format using the [ExportViewsToSHP](<https://github.nrel.gov/SERA/sera-gis/tree/master/Scripts/ExportViewsToSHP.sh>) bash script.
	
Data Maintaniers:

Emily Hueni – Emily.Hueni@nrel.gov

Dana Stright – Dana.Stright@nrel.gov