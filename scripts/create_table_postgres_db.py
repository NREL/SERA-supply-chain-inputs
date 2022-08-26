# -*- coding: utf-8 -*-
"""
Created on Thursay July 20, 2017

@author: ehueni

Purpose: to read excel table into dataframe and upload to postgres table
"""

from sqlalchemy import create_engine
import pandas as pd
import psycopg2
import numpy as np


UserName = 'ehueni'
Password =  '*******'
Server = '1lp11dnpdb1.nrel.gov'
Database = 'sera'
Schema = 'supply_chain_outputs'

#Read File and list columns
df = pd.read_excel("C:\\Users\\ehueni\\Documents\\GitHub_Final\\CECB\\data_2017\\Vision Vehicles Data_Updated_August2017_v2.xlsx", sheetname='HDV')
NewTableName = ''


#print(df)

#Creates new table in database with all of the input data
connection_info = 'postgresql://%s:%s@%s/%s' % (UserName, Password, Server, Database)
engine = create_engine(connection_info)
df.to_sql(NewTableName, engine, schema=Schema)

print (connection_info)