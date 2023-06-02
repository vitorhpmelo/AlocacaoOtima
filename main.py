"""
Código para alocação ótima de medidores no SIN, implemetantação em python
"""






#%%
from classes import *
from readfiles import *
from networkstruc import *
from networkcalc import *
from inidiviuos import *
import numpy as np

lambda_param = 2.5  # Mean of the Poisson distribution
size = 10  # Size of the vector


#%%




#%%

sys="IEEE6"

dfDBAR,dfDBRAN,dfDMED,dfDMEDS = read_files(sys)

#%%
[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)

#
[ram,nbran]=create_bran(dfDBRAN,ind_i)

graph=create_graph(bars,ram)

N=100
HT=montaH(graph,dfDMEDS)
## programa do vigliassi vai me fornecer isso daqui
fitas=[]
for i in range(N):
    fitas.append(np.random.binomial(n=1, p=0.3, size= len(dfDMEDS)))


fitas[0]=np.array([0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0])
#%%    
populacao=[]

## avalia o individuo inicial 


## cria individuos e analisa os individuos
for fita in fitas:
    indv=ciraindIviduos(dfDMEDS,fita)
    populacao.append(indv)# ordena o DMED


#%%
populacao_nc=[]

dTabela={}
i=0
for indv in populacao:
    #esta colocando individuos vazios
    populacao_nc.append(fatoraH(HT,graph,indv))


# %%
