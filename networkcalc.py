from classes import *
from readfiles import *
from networkstruc import *
import numpy as np

tolpiv=1e-9




def montaH(graph,plano):
    i=0
    HT=np.zeros((len(graph),len(plano)))
    Oorder=np.zeros(len(plano))
    for med in plano:
        # if the measurement is a Flow
        if med.tipo==2:
            HT[med.de][i]=1
            HT[med.para][i]=-1
        # if the measurement is a injection 
        elif med.tipo==0:
            k=med.de
            m=0
            for item in graph[k].ladjk:
                HT[item][i]=-1
                m=m+1
            for item in graph[k].ladjm:
                HT[item][i]=-1
                m=m+1
            HT[k][i]=m
        Oorder[i]=i
        i=i+1
    return [HT,Oorder]



def fatoraH(HT,graph,Oorder,plano):
    obs=1
    for  i in range(len(graph)-1):
        if np.abs(HT[i][i])<tolpiv:
            for j in range(i+1,len(plano)):
                obs=0
                if np.abs(HT[i][j])>tolpiv:
                    HT[:,[j,i]]=HT[:,[i,j]]
                    aux=Oorder[i]
                    Oorder[i]=Oorder[j]
                    Oorder[j]=aux
                    obs=1
                    break
            if obs==0:
                break
        for j in range(i+1,len(graph)):
            if np.abs(HT[j][i])>tolpiv:
                HT[j,:]=HT[j,:]-(HT[j][i]/HT[i][i])*HT[i,:]

    i=0
    lstMC=[]
    nMCs=0
    if obs==0:
        return HT,obs,nMCs,lstMC
    for row in HT:
        if np.sum(np.abs(row)>tolpiv)==1:
            lstMC[nMCs]=Oorder[i]
            nMCs=nMCs+1
    i=i+1


    return HT,obs,nMCs,lstMC
