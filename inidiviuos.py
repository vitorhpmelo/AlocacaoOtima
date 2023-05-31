from classes import *
from readfiles import *
from networkstruc import *
import numpy as np

def ciraindIviduos(lstIndv,ind_i):
    """
    Função auxiliar para criar individuos na fase de testes
    dentro dela está a lista de individuos
    """



    indiviuos=[]
    
    for dfDINVI in lstIndv:
        plano=[]
        flagPMU=0
        for idx, row in dfDINVI.iterrows():
            med=medida(row["instalado"],row["type"],ind_i[row["de"]],ind_i[row["para"]])
            plano.append(med)
        ind=individuo(plano,flagPMU)
        indiviuos.append(ind)

    return indiviuos

