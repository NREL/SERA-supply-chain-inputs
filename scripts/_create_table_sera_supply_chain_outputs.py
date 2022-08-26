# -*- coding: utf-8 -*-
"""
Created on Monday August 28, 2017

@author: ehueni

Operating System: cross-platform
"""
from geoalchemy2 import Geometry, WKTElement
from sqlalchemy import create_engine
import pandas as pd
import geopandas as gpd
import numpy as np
import os
import datetime
from shapely.geometry import Point
import glob


#---------------------------User Defined Information------------------------
user_name = 'sera_outputs_insert'
user_password = 'N3W_dAtA!'

#place python file in directory where tsv files or define pathway to parent directory tests
cwd = os.getcwd()
cwd = 'C:\\Users\\ehueni\\Documents\\GitHub_Final\\supply-chain-inputs\\tests'
os.chdir(cwd)
print (cwd)
#---------------------------End User Defined Information------------------------


#Create time stamp for mountain time (minus 6 hours from UTC)
timestamp = pd.to_datetime('now') - datetime.timedelta(hours=6)

output_tables = ['cash.tsv', 'construction.tsv', 'flow.tsv', 'geometry.tsv', 'impact.tsv', 'sale.tsv']

# engine to write dataframe to postgres
engine = create_engine('postgresql://' + user_name + ':' + user_password + '@1lp11dnpdb1.nrel.gov/sera')


for dir in os.listdir(cwd):
    scenario_name = dir
    print ('---------------------------------', dir, '-------------------')
    for dirpath, dirnames, filenames in os.walk(dir):
        for file in [f for f in filenames if f.endswith(".tsv")]:
            if file in output_tables:
                print (file)
                sub_scenario_name = dirpath

                filepath = os.path.join(dirpath, file)
                table = file[:-4]

                # Read File
                df = pd.read_csv(filepath, sep='\t', lineterminator='\n')

                # clean column names to remove extra space at end of names
                df.rename(columns=lambda x: x.strip(), inplace=True)

                # create primary key string with quotes around column name
                column_names = ['"' + l + '"' for l in list(df)]
                print(column_names)

                pk = ''
                if table == 'construction':
                    pk = column_names[0]
                elif table == 'impact':
                    pk = ', '.join(column_names[:4])
                elif table == 'cash':
                    pk = ', '.join(column_names[:3])
                else:
                    pk = ', '.join(column_names[:2])

                # add new columns to df
                df['time_stamp'] = timestamp
                df['scenario'] = dir
                df['sub_scenario'] = sub_scenario_name

                # spatially enable geometry table and append data
                if table == 'geometry':
                    the_geometry = [Point(xy) for xy in zip(df.X, df.Y)]
                    geodataframe = gpd.GeoDataFrame(df, crs=4326, geometry=the_geometry)
                    geodataframe['geom'] = geodataframe['geometry'].apply(lambda x: WKTElement(x, srid=4326))
                    geodataframe.drop('geometry', 1, inplace=True)
                    geodataframe.to_sql(table, engine, schema="supply_chain_outputs", if_exists='append', index=False,
                                        dtype={'geom': Geometry('POINT', srid=4326)})
                else:
                    df.to_sql(table, engine, schema="supply_chain_outputs", if_exists='append', index=False)

                # Add primary key, views, and permissions - will only work if this is creating a new table
                try:
                    sql_alter_table = '''ALTER TABLE supply_chain_outputs.%s OWNER TO "sera-rw";
                    GRANT ALL ON TABLE supply_chain_outputs.%s TO "sera-rw";
                    GRANT SELECT ON TABLE supply_chain_outputs.%s TO public;
                    GRANT SELECT ON TABLE supply_chain_outputs.%s TO "sera-ro";''' % (table, table, table, table)
                    engine.execute(sql_alter_table)

                    sql_pk = '''ALTER TABLE supply_chain_outputs.%s ADD PRIMARY KEY(%s, "scenario", "sub_scenario", "time_stamp");''' % (table, pk)
                    engine.execute(sql_pk)

                    
					#Important Note: When creating new views please take  into account that the geometry_current may need to undergo additional joins to other geometry tables. Please see geometry_current.sql for latest version of the geometry SQL
					sql_create_views = '''CREATE MATERIALIZED VIEW supply_chain_outputs.%s_current AS SELECT * FROM supply_chain_outputs.%s a WHERE "time_stamp" = (SELECT MAX("time_stamp")FROM supply_chain_outputs.%s)''' % (table, table, table)
                    engine.execute(sql_create_views)

                    sql_create_views_alter_table = '''ALTER TABLE supply_chain_outputs.%s_current OWNER TO "sera-rw";
                    GRANT ALL ON TABLE supply_chain_outputs.%s_current TO "sera-rw";
                    GRANT SELECT ON TABLE supply_chain_outputs.%s_current TO public;
                    GRANT SELECT ON TABLE supply_chain_outputs.%s_current TO "sera-ro";''' % (table, table, table, table)
                    engine.execute(sql_create_views_alter_table)
                except:
                    pass


#Refresh Materialized Views
for j in output_tables:
    x = j[:-4]
    print(x)
    sql_refresh_mat_view = '''REFRESH MATERIALIZED VIEW supply_chain_outputs.%s_current;''' % (x)
    print(sql_refresh_mat_view)
    engine.execute(sql_refresh_mat_view)




print(timestamp)