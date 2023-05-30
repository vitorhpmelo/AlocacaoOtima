#%%
from classes import *
from readfiles import *
from networkstruc import *
import numpy as np
from scipy.sparse import csr_matrix
import scikits.umfpack as um
import inspect as inspect

#%%

umfpack = um.UmfpackContext()


umfpack.control[um.UMFPACK_STRATEGY]=um.UMFPACK_STRATEGY_SYMMETRIC
umfpack.control[um.UMFPACK_ORDERING]=um.UMFPACK_ORDERING_NONE


umfpack.control[um.UMFPACK_SINGLETONS]=1


umfpack.control[um.UMFPACK_FIXQ]=1
umfpack.control[um.UMFPACK_AGGRESSIVE]=0
umfpack.control[um.UMFPACK_PIVOT_TOLERANCE]=1
umfpack.control[um.UMFPACK_IRSTEP]=0
umfpack.control[um.UMFPACK_DROPTOL]=1



#%%

sys="IEEE6"

dfDBAR,dfDBRAN,dfDMED,dfDINVI = read_files(sys)


sort_order = {1: 0, 2: 1, 0: 2}



[bars,nbars,pv,pq,ind_i]=creat_bar(dfDBAR)
[ram,nbran]=create_bran(dfDBRAN,ind_i)

graph=create_graph(bars,ram)


dfDINVI["sort_order"]=dfDINVI["instalado"].map(sort_order)
dfDINVI=dfDINVI.sort_values("sort_order")

HT=np.zeros((len(graph),len(dfDINVI)))
Oorder=np.zeros(len(dfDINVI))
#%%
i=0
for idx, row in dfDINVI.iterrows():
    # if the measurement is a Flow
    if row["type"]==2:
        HT[ind_i[row["de"]]][i]=1
        HT[ind_i[row["para"]]][i]=-1
    # if the measurement is a injection 
    elif row["type"]==0:
        k=ind_i[row["de"]]
        m=0
        for item in graph[k].ladjk:
            HT[item][i]=-1
            m=m+1
        for item in graph[k].ladjm:
            HT[item][i]=-1
            m=m+1
        HT[k][i]=m
    i=i+1
    Oorder[i-1]=i

# %%

mtx=csr_matrix(HT.T)
# b=umfpack.symbolic( mtx )
# a=umfpack.numeric( mtx )
#%%
# L, U, P, Q, R, do_recip = umfpack.get_numeric(a)
L, U, P, Q, R, do_recip = umfpack.lu( mtx )

print(P)
# %%
