from classes import *
from readfiles import *
from networkstruc import *
import numpy as np

tolpiv=1e-9




def montaH(graph,ind):

    """
    Função para montar a matriz H
    @param: graph lista de objetos do tipo nó com todas as indormações sobre a rede
    @param: ind objeto do tipo inividuo com o plano de medições e os indicadores de fit
    @return: HT Matriz Jacobiana transposta
    @return: Oorder vetor com a ordenação das medidas   
    """
    
    i=0
    plano=ind.plano

    nvar=len(graph)+ind.FlagPMUV
        
    HT=np.zeros((nvar,len(plano)))

    Oorder=np.zeros(len(plano))
    for med in plano:
        # if the measurement is a Flow
        if med.tipo==2:
            HT[med.de][i]=1
            HT[med.para][i]=-1
        # if the measurement is a injection 
        elif (med.tipo==0) | (med.tipo==6):
            k=med.de
            m=0
            for item in graph[k].ladjk:
                HT[item][i]=-1
                m=m+1
            for item in graph[k].ladjm:
                HT[item][i]=-1
                m=m+1
            HT[k][i]=m
        elif med.tipo==5:
            HT[med.de][i]=1
   
        Oorder[i]=i
        i=i+1
        
    return [HT,Oorder]



def fatoraH(HT,graph,Oorder,ind):
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


    # 
    for  i in range(len(graph)-1):
     # Verifica se o pivo é nulo
        if np.abs(HT[i][i])<tolpiv:
            # se for ele procura a próxima medida com pivo não nulo
            for j in range(i+1,len(plano)):
                if np.abs(HT[i][j])>tolpiv:
                    # quando ele contra a primeira medida com pivo não nulo e faz a troca
                    HT[:,[j,i]]=HT[:,[i,j]]
                    aux=Oorder[i]
                    Oorder[i]=Oorder[j]
                    Oorder[j]=aux
                    
                    # altera plano, tem q reordenar o plano pra corresponder a jacobiana e a ao Oorder
                    if plano[j].instalado==2:# se precisou de 1 MFS candidata para a observabilidade só incrementa
                        nMFScandidatas=nMFScandidatas+1# se precisou de 1 MFS candidata para a observabilidade só incrementa                       
                    if plano[j].instalado==0: # se precisou de 1 PMU candidata 
                        nPMUscandidatas=nPMUscandidatas+1 # incrementa 
                        startPMUS
                        for k in range(startPMUS,len(plano)): # atualiza todas as medidas daquela PMU para MFS
                            if plano[k].de==plano[j].de: # se a medida está na mesma barra
                               plano[k].instalado==2 # muda o status para candidata
                               aux=plano[startPMUS]  
                               plano[startPMUS]=plano[k] # coloca ela no fim do arquivo de candidatas
                               plano[k]=aux # coloca a primeira PMU candidata pro lugar dela
                               HT[:,[k,startPMUS]]=HT[:,[startPMUS,k]]# altera a ordem da Jacobiana
                               startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS
                            #muda o status de todas as outras medidas na mesma barra para candidata
                            #reordena a Ht e o plano para elas virem primeiro                                     
                    aux=plano[i] # mantem o plano na mesma ordem da HT
                    plano[i]=plano[j]  # mantem o plano na mesma ordem da HT
                    plano[j]=aux  # mantem o plano na mesma ordem da HT
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
