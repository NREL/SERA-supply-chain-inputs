create function pg_temp.getfraction(_tbl_1 text, _tbl_1_col text, _tbl_1_geom text, _tbl_2 text, _tbl_2_col text, _tbl_2_geom text)
returns table ("ID" character varying(75), "Zone" character varying(75), "Overlap" boolean, "Fraction" double precision)
as $func$
BEGIN
Return QUERY EXECUTE
'
SELECT ID, Zone, Overlap, final."Fraction" FROM 
(SElect *, CASE 
	WHEN ST_Area(st_intersection(geom, ClipToGeom))=0 then 0
	ELSE (ST_Area(ST_Intersection(st_intersection(geom, ClipToGeom), geom_2))/ST_Area(st_intersection(geom, ClipToGeom))) 
	END as "Fraction"
 FROM
(SELECT a.'||_tbl_1_col||' as ID, a.'||_tbl_1_geom||' as geom, b.'||_tbl_2_col||' as Zone, b.'||_tbl_2_geom||' as geom_2, st_overlaps(a.'||_tbl_1_geom||', b.'||_tbl_2_geom||') as Overlap,
case
when st_overlaps(a.'||_tbl_1_geom||', b.'||_tbl_2_geom||') OR st_contains(a.'||_tbl_1_geom||', b.'||_tbl_2_geom||') then 
	(select st_union (b.'||_tbl_2_geom||') from sera_gis.'||_tbl_2||' b Where st_intersects(a.'||_tbl_1_geom||', b.'||_tbl_2_geom||') group by a.'||_tbl_1_col||')
else b.'||_tbl_2_geom||'
end as ClipToGeom 
FROM sera_gis.'||_tbl_1||' a, sera_gis.'||_tbl_2||' b
where st_intersects(a.'||_tbl_1_geom||', b.'||_tbl_2_geom||'))intermediate) 
Final
';
END;
$func$
Language plpgsql;

Create Temporary Table a as (SELECT * from pg_temp.getfraction('"Census Urban Area Codes 2010"','"Census Urban Area Code"::character varying', '"Census Urban Area Geometry 500k"', '"Census Counties 2010"', '"County Name"', '"Census County Geometry 20m"'))

--SELECT * from pg_temp.getfraction('"BAs 2014"', '"BA Name"', '"BAs Geometry"','"Census Counties 2010"', '"County Name"', '"Census County Geometry 20m"')

--SELECT * from pg_temp.getfraction('"Census Counties 2010"', '"County Name"', '"Census County Geometry 20m"', '"BAs 2014"', '"BA Name"', '"BAs Geometry"')