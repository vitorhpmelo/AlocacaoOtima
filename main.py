#%%
from classes import *
from readfiles import *
from networkstruc import *
from networkcalc import *
from inidiviuos import *
import numpy as np







#%%
tolpiv=1e-9


#%%

sys="IEEE6"

dfDBAR,dfDBRAN,dfDMED,dfDINVI = read_files(sys)


[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

graph=create_graph(bars,ram)


sort_order = {1: 0, 2: 1, 0: 2}

dfDINVI["sort_order"]=dfDINVI["instalado"].map(sort_order)
dfDINVI=dfDINVI.sort_values("sort_order")

lstInd=[]
lstInd.append(dfDINVI)

#%%



individuos=ciraindIviduos(lstInd,ind_i)



[HT,Oorder]=montaH(graph,individuos[0])


[Htri,obs,nMCs,MC]=fatoraH(HT,graph,Oorder,individuos[0])
   

# %%



