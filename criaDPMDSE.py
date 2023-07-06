"""
Arquivo de apoio para criar o DPMFSE de acordo com os branchs da redem contando apenas com PMUS
O arquivo é ordenado por barra (subestação)
Toda barra tem medidas nos correntes nos seus ramos adjacentes 
Toda barra tem medidas de tensão nos seus ramos adjacentes 
"""

#%%
import pandas as pd
import numpy as np



sys="IEEE6"

dfBran=pd.read_csv(sys+"/DBRAN.csv",header=None)
dfBran.columns=["id","tipo","de","para","r","x","bsh","tap"]


diretos=dfBran[["de","para"]].values
inversos=dfBran[["para","de"]].values

#%%
todos=np.vstack((diretos,inversos))

Imeas=np.ones(len(todos),dtype=int)*6
Imeas=np.column_stack([Imeas,todos])
Vmeas=np.ones(len(todos),dtype=int)*5
Vmeas=np.column_stack([Vmeas,todos])
DPMFSE=np.vstack([Imeas,Vmeas])

dfDMES=pd.DataFrame(data=DPMFSE,columns=["tipo","de","para"])
dfDMES.sort_values(by=["de","para"],inplace=True,ignore_index=True)
dfDMES["existente"]=0


dfDMES["id_MFS"]=0
#%%
i=0
for item in todos:
    mask=(dfDMES["de"]==item[0]) & (dfDMES["para"]==item[1]) 
    dfDMES.loc[mask,"id_MFS"]=i
    i=i+1


# %%
dfDMES.to_csv("teste.csv",header=None,index=None,columns = ["existente","tipo","de","para","id_MFS"])


# %%
