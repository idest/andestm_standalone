import pandas as pd
from pandas import read_excel, ExcelWriter

#eq1 = read_excel("CSN_2000-2012_CAT.xlsx", sheet_name="Hoja1")
#eq2 = read_excel("CSN_2013-2018_CAT.xlsx", sheet_name="Hoja1")
#
#mask1_lon = eq1['Lon.']<-1000
#mask1_lat = eq1['Lat.']<-1000
#eq1.loc[mask1_lon, 'Lon.'] = eq1.loc[mask1_lon, 'Lon.']/1000.
#eq1.loc[mask1_lat, 'Lat.'] = eq1.loc[mask1_lat, 'Lat.']/1000.
#
#mask2_lon = eq2['Lon.']<-1000
#mask2_lat = eq2['Lat.']<-1000
#eq2.loc[mask2_lon, 'Lon.'] = eq2.loc[mask2_lon, 'Lon.']/1000.
#eq2.loc[mask2_lat, 'Lat.'] = eq2.loc[mask2_lat, 'Lat.']/1000.

#writer = ExcelWriter('CSN.xlsx')
#eq1.to_excel(writer,'Sheet1')
#eq2.to_excel(writer,'Sheet2')
#writer.save()

eq1 = read_excel('data/earthquakes/CSN.xlsx', sheet_name='Sheet1')
eq2 = read_excel('data/earthquakes/CSN.xlsx', sheet_name='Sheet2')

eqs = pd.concat([eq1, eq2], ignore_index=True)

writer = ExcelWriter('data/earthquakes/CSN_2000_2018.xlsx')
eqs.to_excel(writer, 'Sheet1')
writer.save()

"""
earthquakes = eq1
latitude = -18.
eq = earthquakes[
    (earthquakes['Lat.'] >= latitude-0.1) &
    (earthquakes['Lat.'] < latitude+0.1)]
print(eq1.head())
print(eq2.head())
print(eq)
"""
