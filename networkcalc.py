from classes import *
from readfiles import *
from networkstruc import *
import numpy as np

tolpiv=1e-9


def reorder_columns(matrix, permutation_vector):
    # Convert matrix and permutation vector to numpy arrays
    matrix = np.array(matrix)
    permutation_vector = np.array(permutation_vector)
    
    # Get the number of columns in the matrix
    num_columns = matrix.shape[1]
    
    # Check if permutation vector is valid
    if len(permutation_vector) != num_columns:
        raise ValueError("Permutation vector length does not match the number of columns in the matrix.")
    
    # Reorder the columns of the matrix
    reordered_matrix = matrix[:, permutation_vector]
    
    return reordered_matrix

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
    melhorar, na permutacao pode mandar uma medida instalada lá para trás 
    """


    nMFScandidatas=0
    nPMUscandidatas=0

    plano=ind.plano
    startPMUS=ind.startPMUS

    obs=1

    Hfat=HT.copy()

    Hfat=reorder_columns(Hfat,plano.Oorder.array)

    flag=ind.FlagPMUV
    Permutacao=list(plano.index) # vetor de permutacoes tem a mesma ordem da Htriang final e em suas poiscoes tem a posicao da medida no plano
    plano["instalado_mod"]=plano["instalado"]
    # 
    for  i in range(len(graph)-1+flag):
     # Verifica se o pivo é nulo
        if np.abs(Hfat[i][i])<tolpiv:
            # se for ele procura a próxima medida com pivo não nulo
            for j in range(i+1,len(plano)):
                if np.abs(Hfat[i][j])>tolpiv:
                    Hfat[:,[j,i]]=Hfat[:,[i,j]] # quando ele contra a primeira medida com pivo não nulo e faz a troca i por j
                    aux=Permutacao[i] # coloca permutacao no vetor de permutacoes
                    Permutacao[i]=Permutacao[j] # coloca permutacao no vetor de permutacoes
                    Permutacao[j]=aux # coloca permutacao no vetor de permutacoes
                    ## Depois daqui o i guarda a posicao da medida permutada e o j da medida que estava na posicao anteriormente
                    if plano.iloc[Permutacao[i]].instalado_candidatas==2: #verifica se a medida permutada era candidata
                        nMFScandidatas=nMFScandidatas+1# se precisou de 1 MFS candidata para a observabilidade só incrementa                       
                        plano.at[Permutacao[i],"instalado_mod"]=1 # muda o status dela para mod
                    elif plano.iloc[Permutacao[i]].instalado_candidatas==0: # se precisou de 1 PMU candidata 
                        nPMUscandidatas=nPMUscandidatas+1 # incrementa 
                        plano.at[Permutacao[i],"instalado_mod"]=1 #muda o status das medidas
                        offset=startPMUS
                        for k in range(offset,len(plano)): # atualiza todas as medidas daquela PMU para MFS
                            if plano.iloc[Permutacao[k]].i_de==plano.iloc[Permutacao[i]].i_de: # se a medida está na mesma barra
                               plano.at[Permutacao[k],"instalado_candidatas"]=2 # muda o status para candidata
                               Hfat[:,[k,startPMUS]]=Hfat[:,[startPMUS,k]]# altera a ordem da Jacobiana
                               aux=Permutacao[startPMUS]  
                               Permutacao[startPMUS]=Permutacao[k] # coloca ela no fim do arquivo de candidatas
                               Permutacao[k]=aux # coloca a primeira PMU candidata pro lugar dela
                               startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS
                            #muda o status de todas as outras medidas na mesma barra para candidata
                            #reordena a Ht e o plano para elas virem primeiro
                    if plano.iloc[Permutacao[j]].instalado_candidatas==1:
                        #se a medida for permutada eu tenho q passar ela pro topo das candidatas
                        pass

                    break
        for j in range(i+1,len(graph)):
            if np.abs(Hfat[j][i])>tolpiv:
                Hfat[j,:]=Hfat[j,:]-(Hfat[j][i]/Hfat[i][i])*Hfat[i,:]  


    ## ordenar a partir das medidas básixas para colocar as instaldas primeiro
    ## Verificar MCs a partir das medidas básicas
    ## ordenar só no dataframe



    i=0
    lstMC=[]
    nMCs=0
    # ordernar a base da matriz pelas instaladas 
    for row in Hfat:
        if np.sum(np.abs(row)>tolpiv)==1:
            lstMC[nMCs]=Permutacao[i]
            nMCs=nMCs+1
        i=i+1
    ind.nMCs=nMCs
    ind.nMedidas_adicionadas=nMFScandidatas
    ind.nMedidas_adicionadas=nPMUscandidatas

    ## falta fazer isso daqui verificar a partir da matriz fatorada, quais as medidas que devem ser instaladas para não ter mais medidas críticas    
    # if ind.nMCs>0:


    return Hfat,obs,nMCs,lstMC
