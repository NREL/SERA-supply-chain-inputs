-- Important Note: Use nearest neighbor join and not intersect.  Along the coast the geometries do not always match  up  perfectly and so an st_intersect may misss some of the data (for example, the water stress ares)

-- Materialized View: supply_chain_outputs.geometry_current

-- DROP MATERIALIZED VIEW supply_chain_outputs.geometry_current;

CREATE MATERIALIZED VIEW supply_chain_outputs.geometry_current AS 
select DISTINCT ON("Network ID") a.*, b."Baseline Water Stress Category" From (SELECT a."Network ID",
    a."Position",
    a."X",
    a."Y",
    a.time_stamp,
    a.scenario,
    a.sub_scenario,
    a.geom
   FROM supply_chain_outputs.geometry a
  WHERE a.time_stamp = (( SELECT max(geometry.time_stamp) AS max
           FROM supply_chain_outputs.geometry))) a, sera_gis."Water Stress Areas" b where st_dwithin(st_transform(a.geom, 5070), b."Water Stress Geometry", 4000) order by "Network ID", st_distance(st_transform(a.geom, 5070), b."Water Stress Geometry")
WITH DATA;

ALTER TABLE supply_chain_outputs.geometry_current
  OWNER TO "sera-rw";
GRANT ALL ON TABLE supply_chain_outputs.geometry_current TO "sera-rw";
GRANT SELECT ON TABLE supply_chain_outputs.geometry_current TO public;
GRANT SELECT ON TABLE supply_chain_outputs.geometry_current TO "sera-ro";

-- Index: supply_chain_outputs.index_geometry_current_sub_scenario

-- DROP INDEX supply_chain_outputs.index_geometry_current_sub_scenario;

CREATE INDEX index_geometry_current_sub_scenario
  ON supply_chain_outputs.geometry_current
  USING btree
  (sub_scenario COLLATE pg_catalog."default");

