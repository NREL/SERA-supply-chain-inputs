# -*- coding: utf-8 -*-
"""
Created on Monday September 26, 2017

@author: ehueni

Purpose: Example of reading an excel spreadsheet and outputing a tsv file that is formatted for sera (encoding, not quotes, correct line  endings, ect.)
"""

import pandas as pd


#AdjustableVariables
WorkingDir = "C:\\Users\\ehueni\\Documents\\CreateNewLinksAndNodes\\"


Workbook = WorkingDir+"CanadaData_Working.xlsx" #"CanadaData_Working.xlsx", "NewLinks.xlsx", "links_old_2.xlsx"
Sheet_Name = 'Nodes_Working' #Nodes_Working, Existing_Working, 'NewLinks'
Output = WorkingDir+Sheet_Name+'_Final.tsv'

df = pd.read_excel(Workbook, sheetname=Sheet_Name)
df.drop('City', axis=1, inplace = True)
#df['Capacity [kg/yr]'] = df['Capacity [kg/yr]'].astype(int)

print(df)

df.to_csv(Output, sep='\t', encoding='utf-8', index=False, line_terminator='\n')