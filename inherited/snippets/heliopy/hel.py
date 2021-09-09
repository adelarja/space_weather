"""
El dataframe 'data' tiene las siguientes columnas:

['BF1', 'BF1LOG', 'BRMSF1', 'BGSM_0', 'BGSM_1', 'BGSM_2', 'BRMSGSM_0', 'BRMSGSM_1', 'BRMSGSM_2', 'BGSE_0', 'BGSE_1', 'BGSE_2', 'BGSEa_0', 'BGSEa_1', 'BGSEa_2', 'BRMSGSE_0', 'BRMSGSE_1', 'BRMSGSE_2', 'DIST', 'PGSM_0', 'PGSM_1', 'PGSM_2', 'PGSE_0', 'PGSE_1', 'PGSE_2', 'DISTV', 'PGSMV_0', 'PGSMV_1', 'PGSMV_2', 'PGSEV_0', 'PGSEV_1', 'PGSEV_2']

Se exporta un csv con Time, BGSM_0, BGSM_1, y BGSM_1 separados por coma.

No estoy seguro de si lo que necesitamos es BGSM:
Magnetic field vector in GSM cartesian coordinates (1 min) [BGSM]
"""

import heliopy.data.wind as wind
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

starttime = datetime(2016, 1, 1, 0, 0, 0)
endtime = starttime + timedelta(hours=2)

data = wind.mfi_h0(starttime, endtime)

print(data.columns)
data = data.to_dataframe()

vector = data.loc[:,['BGSM_0', 'BGSM_1', 'BGSM_2']]

print(vector)

with open('file.csv', 'w') as f:
    f.write(vector.to_csv())