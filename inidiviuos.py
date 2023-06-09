from classes import *
from readfiles import *
from networkstruc import *
import numpy as np











def ciraindIviduos(dfDMEDS,fita):
    """
    Função auxiliar para criar individuos na fase de testes
    dentro dela está a lista de individuos
    """

    sort_order = {1: 0, 2: 1, 0: 2}

    DMEDS_indviduo=dfDMEDS.copy()#cria o DMED no individuo
    DMEDS_indviduo["fita"]=fita #insere a fita
    DMEDS_indviduo["instalado_fita"]=fita

    ## garante que os pares I e V ramos são instalados

    MFS_instaladas=list(set(DMEDS_indviduo[(DMEDS_indviduo["fita"]==1)].id_MFS))
    mask=DMEDS_indviduo["id_MFS"].isin(MFS_instaladas)
    DMEDS_indviduo.loc[mask,"instalado_fita"]=1


    # força as medidas existentes
    DMEDS_indviduo["instalado"]=(DMEDS_indviduo["instalado_fita"]) | (DMEDS_indviduo["existente"]) # força as medidas existentes coisas existentes 
    
    
    DMEDS_indviduo["instalado_candidatas"]=DMEDS_indviduo["instalado"] # cria a coluna instalado & candidatas
    barras_c_candidatas=list(set(DMEDS_indviduo[DMEDS_indviduo["instalado"]==1].i_de)) # cria lista de barras com medida ("PMUS")
    mask=(DMEDS_indviduo["i_de"].isin(barras_c_candidatas))&(DMEDS_indviduo["instalado"]==0)
    DMEDS_indviduo.loc[mask,"instalado_candidatas"]=2 #Insere as PMUs Candidadas
    DMEDS_indviduo["Oorder_DPFMSE"]=DMEDS_indviduo.index.to_list() # ordenação original do dataframe Oorder_DPFMSE
    indv=individuo(lista=fita,plano=DMEDS_indviduo) # cria o invidiuo
    indv.nPMUs_instaladas=len(barras_c_candidatas) # preenche as PMUs instaladas
    indv.nMFS_instaladas=sum(fita) # preenche o numero de MFS instaladas
    indv.FlagPMUV=sum(DMEDS_indviduo[DMEDS_indviduo["type"]==5].instalado)>0 # verifica se tem PMU de tnesão para fatorar até a linha adicional
    DMEDS_indviduo["sort_order"]=DMEDS_indviduo["instalado_candidatas"].map(sort_order) # ordena o DMED 
    DMEDS_indviduo=DMEDS_indviduo.sort_values(["sort_order","Oorder_DPFMSE"], ignore_index=True)# ordena o DMED 
    indv.plano=DMEDS_indviduo # coloca o DMEDS_individuo
    indv.startPMUS=sum((DMEDS_indviduo["instalado_candidatas"]==2) | (DMEDS_indviduo["instalado_candidatas"]==1)) 
    indv.startCandidatas=sum(DMEDS_indviduo["instalado_candidatas"]==1)
    return indv

