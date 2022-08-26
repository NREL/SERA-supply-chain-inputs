# -*- coding: utf-8 -*-
"""
Created on October, 2017

@author: ehueni

Purpose: To reformat biomass spreadsheet sent by EIA in October, 2017 to a format and layout for sera biomass prices..
"""

import pandas as pd



file = "C:\\Users\\ehueni\\Documents\\BiomassData\\AEO2017_BiomassPrices_20171010 (3).xlsx"
sheet = 'AEO2017'
output_name = 'energy-prices-eia-biomass.tsv'
Output = "C:\\Users\\ehueni\\Documents\\GitHub_Final\\supply-chain-inputs\\inputs\\prices\\"+output_name
eia_file_name = 'energy-prices-reference-all-eia-real-dollars.tsv'
rows_of_interest_start = 2
rows_of_interest_stop = 18
ConversionPriceFactor = 17


#---------------------------------------------
#eia_file = "C:\\Users\\ehueni\\Documents\\GitHub_Final\\supply-chain-inputs\\inputs\\prices\\"+eia_file_name


Headers = ['Material', 'Zone', 'Year', 'Price [$/unit]', 'Billable?']

df = pd.read_excel(file, sheet, skiprows=7, parse_cols='B:AN')


df_sub = df.iloc[rows_of_interest_start:rows_of_interest_stop]

print(list(df_sub))

columns_to_pivot =[2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029, 2030, 2031, 2032, 2033, 2034, 2035, 2036, 2037, 2038, 2039, 2040, 2041, 2042, 2043, 2044, 2045, 2046, 2047, 2048, 2049, 2050]
columns_to_keep = ['Type and Region']

df_2 = pd.melt(df_sub, id_vars=columns_to_keep, value_vars=columns_to_pivot, var_name="Year", value_name="Unconverted_Price")

df_2['Price [$/unit]'] = df_2['Unconverted_Price'] * ConversionPriceFactor

df_2['Zone'] = df_2['Type and Region'].str.strip()

#No Averaging---------------------------------------
'''
df_2['Billable?'] = 'True'
df_2['Material'] = 'Biomass'

df_2.to_csv(Output, sep='\t', encoding='utf-8', columns=Headers, index=False, line_terminator='\n')
'''


#Averaging by division ---------------------------------------

for i, j in enumerate(df_2['Zone']):
    if 'East North Cent' in j:
        df_2['Zone'][i] = 'East North Central'
    elif 'South Atlantic' in j:
        df_2['Zone'][i] = 'South Atlantic'
    elif 'East South' in j:
        df_2['Zone'][i] = 'East South Central'
    elif 'West North' in j:
        df_2['Zone'][i] = 'West North Central'
    elif 'Mountain' in j:
        df_2['Zone'][i] = 'Mountain'

df_4 = df_2.groupby(['Zone', 'Year'])['Price [$/unit]'].mean()
df_4 = df_4.to_frame().reset_index()

print(df_4)

df_4['Billable?'] = 'True'
df_4['Material'] = 'Biomass [ton]'

df_4.to_csv(Output, sep='\t', encoding='utf-8', columns=Headers, index=False, line_terminator='\n')
