from classes import *
from readfiles import *
from networkstruc import *
import numpy as np

tolpiv=1e-9




def montaH(graph,dfDMEDS):

    """
    Função para montar a matriz H
    @param: graph lista de objetos do tipo nó com todas as indormações sobre a rede
    @param: DMEDS dataframe com as informações de todas as medidas possiveis
    @return: HT Matriz Jacobiana transposta
    @return: Oorder vetor com a ordenação das medidas   
    """
    
    i=0
    
            
    HT=np.zeros((len(graph)+1,len(dfDMEDS)))
    
    for idx,med in dfDMEDS.iterrows():
        # if the measurement is a Flow
        if med.type==2:
            HT[med.i_de][i]=1
            HT[med.i_para][i]=-1
        # if the measurement is a injection 
        elif (med.type==0) | (med.type==6):
            k=med.i_de
            m=0
            for item in graph[k].ladjk:
                HT[item][i]=-1
                m=m+1
            for item in graph[k].ladjm:
                HT[item][i]=-1
                m=m+1
            HT[k][i]=m
        elif med.type==5:
            HT[med.i_de][i]=1
            HT[med.i_de][i]=-1
        i=i+1
        
    return HT



def fatoraH(HT,graph,ind):
    """
    Função para fatorar a matriz H
    @param: matriz jacobiana HT array by dimensional de float do numpy
    @param: graph lista de objetos do tipo nó com todas as indormações sobre a rede
    @param: Oorder vetor com a ordenação das medidas
    @param: ind objeto do tipo inividuo com o plano de medições e os indicadores de fit
    @return: H_fat matriz fatorada
    @return: Numero de medidas criticas
    @return: N medidas adicionadas 
    """


    nMFScandidatas=0
    nPMUscandidatas=0

    plano=ind.plano
    startPMUS=ind.startPMUS

    obs=1

    Hfat=HT.copy()
    flag=ind.FlagPMUV
    Oorder=list(plano.index) # guarda os index originais do data frame
    plano["instalado_mod"]=plano["instalado"]
    # 
    for  i in range(len(graph)-1+flag):
     # Verifica se o pivo é nulo
        if np.abs(Hfat[i][i])<tolpiv:
            # se for ele procura a próxima medida com pivo não nulo
            for j in range(i+1,len(plano)):
                if np.abs(Hfat[i][j])>tolpiv:
                    # quando ele contra a primeira medida com pivo não nulo e faz a troca i por j
                    Hfat[:,[j,i]]=Hfat[:,[i,j]]
                    aux=Oorder[i] # index do DF plano
                    Oorder[i]=Oorder[j] # index final do DF plano
                    Oorder[j]=aux # index final do DF plano
                    # altera plano, tem q reordenar o plano pra corresponder a jacobiana e a ao Oorder
                    if plano.iloc[Oorder[j]].instalado_candidatas==2:# se precisou de 1 MFS candidata para a observabilidade só incrementa
                        nMFScandidatas=nMFScandidatas+1# se precisou de 1 MFS candidata para a observabilidade só incrementa                       
                        plano.iloc[Oorder[j]].loc["instalado_mod"]=1
                    if plano.iloc[Oorder[j]].instalado_candidatas==0: # se precisou de 1 PMU candidata 
                        nPMUscandidatas=nPMUscandidatas+1 # incrementa 
                        plano.iloc[Oorder[j]].loc["instalado_mod"]=1
                        for k in range(startPMUS,len(plano)): # atualiza todas as medidas daquela PMU para MFS
                            if plano.iloc[Oorder[k]].i_de==plano[j].de: # se a medida está na mesma barra
                               plano.iloc[Oorder[k]].loc["instalado_candidatas"]=2 # muda o status para candidata
                               aux=Oorder[startPMUS]  
                               Oorder[startPMUS]=Oorder[k] # coloca ela no fim do arquivo de candidatas
                               Oorder[k]=aux # coloca a primeira PMU candidata pro lugar dela
                               Hfat[:,[k,startPMUS]]=Hfat[:,[startPMUS,k]]# altera a ordem da Jacobiana
                               startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS
                            #muda o status de todas as outras medidas na mesma barra para candidata
                            #reordena a Ht e o plano para elas virem primeiro                                     
                    break
        for j in range(i+1,len(graph)):
            if np.abs(Hfat[j][i])>tolpiv:
                Hfat[j,:]=Hfat[j,:]-(Hfat[j][i]/Hfat[i][i])*Hfat[i,:]
    #%%pensar como tirar os dados
    
    
    i=0
    lstMC=[]
    nMCs=0
    if obs==0:
        return Hfat,obs,nMCs,lstMC
    for row in Hfat:
        if np.sum(np.abs(row)>tolpiv)==1:
            lstMC[nMCs]=Oorder[i]
            nMCs=nMCs+1
    i=i+1

    return Hfat,obs,nMCs,lstMC
