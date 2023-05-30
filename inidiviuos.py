from classes import *
from readfiles import *
from networkstruc import *
import numpy as np

def ciraindIviduos(lstIndv,ind_i):

    indiviuos=[]
    
    for dfDINVI in lstIndv:
        plano=[]
        for idx, row in dfDINVI.iterrows():
            med=medida(row["instalado"],row["type"],ind_i[row["de"]],ind_i[row["para"]])
            plano.append(med)
        indiviuos.append(plano)

    return indiviuos

