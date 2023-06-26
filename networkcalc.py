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
        # if the measurement is a Flow or a current fwo
        if med.type==2| (med.type==6):
            HT[med.i_de][i]=1
            HT[med.i_para][i]=-1
        # if the measurement is a injection 
        elif (med.type==0):
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
            HT[-1][i]=-1
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
    nMFSadicioanadas=0
    nPMUscandidatas=0
    nPMUSadicionadas=0

    plano=ind.plano
    startPMUS=ind.startPMUS#util na hora da permutacao
    startCandidatas=ind.startCandidatas#util na hora da permutacao
    tollfill=1e-9
    obs=1

    Hfat=HT.copy()

    Hfat=reorder_columns(Hfat,plano.Oorder.array)

    flag=ind.FlagPMUV
    Permutacao=list(plano.index) # vetor de permutacoes tem a mesma ordem da Htriang final e em suas poiscoes tem a posicao da medida no plano
    plano["instalado_mod"]=plano["instalado"]
    # 
    for  i in range(len(graph)-1+flag):
     # Verifica se o pivo é nulo

        if np.abs(Hfat[i][i]) >tolpiv:
            if plano.iloc[Permutacao[i]].instalado_candidatas!=1:
                if plano.iloc[Permutacao[i]].instalado_candidatas==2:
                    plano.at[Permutacao[i],"instalado_candidatas"]=1
                    plano.at[Permutacao[i],"instalado_mod"]=1
                    if plano.iloc[Permutacao[i]].type==5 & flag==0:
                        flag=1
                if plano.iloc[Permutacao[i]].instalado_candidatas==0:
                    plano.at[Permutacao[i],"instalado_candidatas"]=1
                    plano.at[Permutacao[i],"instalado_mod"]=1
                    if (plano.iloc[Permutacao[i]].type==5) & (flag==0):
                        flag=1
                    offset=startPMUS
                    for k in range(offset,len(plano)): # atualiza todas as medidas daquela PMU para MFS
                        if plano.iloc[Permutacao[k]].i_de==plano.iloc[Permutacao[i]].i_de: # se a medida está na mesma barra
                            plano.at[Permutacao[k],"instalado_candidatas"]=2 # muda o status para candidata
                            Hfat[:,[k,startPMUS]]=Hfat[:,[startPMUS,k]]# altera a ordem da Jacobiana
                            aux=Permutacao[startPMUS]  
                            Permutacao[startPMUS]=Permutacao[k] # coloca ela no fim do arquivo de candidatas
                            Permutacao[k]=aux # coloca a primeira PMU candidata pro lugar dela
                        startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS
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
                        nMFSadicioanadas=nMFSadicioanadas+1               
                        plano.at[Permutacao[i],"instalado_mod"]=1 # muda o status dela para mod
                        plano.at[Permutacao[i],"instalado_candidatas"]=1
                        if (plano.iloc[Permutacao[i]].type==5) & (flag==0):
                            flag=1
                        if plano.iloc[Permutacao[j]].instalado_candidatas==1: #coloca a medida que foi retirada para o final das instaladas
                            Hfat[:,[j,startCandidatas]]=Hfat[:,[startCandidatas,j]]
                            aux=Permutacao[startCandidatas]  
                            Permutacao[startCandidatas]=Permutacao[j] # coloca ela no fim do arquivo de candidatas
                            Permutacao[j]=aux # coloca a primeira PMU candidata pro lugar dela    
                        startCandidatas=startCandidatas+1 # incrementa o identificador de onde começam as MFS_candidatas
                    elif plano.iloc[Permutacao[i]].instalado_candidatas==0: # se precisou de 1 PMU candidata 
                        nPMUscandidatas=nPMUscandidatas+1 # incrementa 
                        nMFSadicioanadas=nMFSadicioanadas+1
                        nPMUSadicionadas=nPMUSadicionadas+1
                        plano.at[Permutacao[i],"instalado_mod"]=1 #muda o status das medidas
                        plano.at[Permutacao[i],"instalado_candidatas"]=1
                        if (plano.iloc[Permutacao[i]].type==5) & (flag==0):
                            flag=1
                        if plano.iloc[Permutacao[j]].instalado_candidatas==1:
                            Hfat[:,[j,startCandidatas]]=Hfat[:,[startCandidatas,j]] # volta a medida para o topo das candidatas
                            aux=Permutacao[startCandidatas]  
                            Permutacao[startCandidatas]=Permutacao[j] # coloca ela no fim do arquivo de candidatas
                            Permutacao[j]=aux # coloca a primeira PMU candidata pro lugar dela
                            startCandidatas=startCandidatas+1 # incrementa o identificador de onde começam as PMUS
                            Hfat[:,[j,startPMUS]]=Hfat[:,[startPMUS,j]] # volta a medida do topo das candidatas para as candidatas
                            aux=Permutacao[startPMUS]  
                            Permutacao[startPMUS]=Permutacao[j] # coloca ela no fim do arquivo de candidatas
                            Permutacao[j]=aux # coloca a primeira PMU candidata pro lugar dela
                            startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS
                        startCandidatas=startCandidatas+1 # incrementa o identificador de onde começam as MFS_candidatas
                        startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS
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
                    break
        Hfat[i,:]=Hfat[i,:]/Hfat[i][i]
        for j in range(i+1,len(graph)):
                if np.abs(Hfat[j][i])>tollfill:
                    Hfat[j,:]=Hfat[j,:]-(Hfat[j][i])*Hfat[i,:]  
                else:  
                    Hfat[j][i]=0  


    for i in range(1,(len(graph)-1+flag)):
        for j in range(0,i):
            if np.abs(Hfat[j][i])>tollfill:
                Hfat[j,:] = Hfat[j,:]-Hfat[j,i]*Hfat[i,:]
            else:
                Hfat[j][i]=0




    plano["order_after_observa"]=0
    #desinstala medidas



    ## verifica PMUs instaladas no conjunto básico
    PMUsInstaladas=list(set(plano.iloc[Permutacao[0:len(graph)-1+flag]].i_de.tolist()))

    #
    mask1=plano.iloc[Permutacao[len(graph)-1+flag:]].i_de.isin(PMUsInstaladas) 
    mask1=(~(plano.de !=-999) )| mask1

    plano.loc[mask1,"instalado_candidatas"]=2
    plano.loc[mask1,"instalado_mod"]=0 #para todas
   
    ##desinstala todas as medidas depois da observabilidade

    mask2=~plano.iloc[Permutacao[len(graph)-1+flag:]].i_de.isin(PMUsInstaladas) 
    mask2=(~(plano.de !=-999) )| mask2

    plano.loc[mask2,"instalado_candidatas"]=0
    plano.loc[mask2,"instalado_mod"]=0 #para todas

    #ordem depois da observabilidade
    plano.loc[Permutacao,"order_after_observa"]=list(range(len(Permutacao)))

    plano=plano.sort_values("order_after_observa")

    sort_order = {1: 0, 2: 1, 0: 2}
    plano["sort_order_observabilidade"]=plano["instalado_candidatas"].map(sort_order) # ordena o DMED 
    plano=plano.sort_values("sort_order_observabilidade",kind="stable")# ordena o DMED 

    # aqui está errado o vetor permutacao depois da observabilidade 
    Permutacao_depois_Observabilidade=plano["order_after_observa"].tolist()

    Hfat2=reorder_columns(Hfat,Permutacao_depois_Observabilidade)
    Permutacao=np.array(Permutacao)
    Permutacao_depois_Observabilidade=Permutacao[Permutacao_depois_Observabilidade]








    

    ncadidatas=len(plano[plano["instalado_candidatas"]==2])
    nPMUscandidatas=len(plano[plano["instalado_candidatas"]==0])
    
    startCandidatas=len(plano)-ncadidatas-nPMUscandidatas
    startPMUS=len(plano)-ncadidatas



    i=0
    lstMC=[]
    varMC=[]
    nMCs=0
    # ordernar a base da matriz pelas instaladas 
    for i in range(0,len(graph)-1+flag):
        if np.sum(np.abs(Hfat2[i][0:startCandidatas])>tolpiv)<=1:
            lstMC.append(plano.iloc[Permutacao_depois_Observabilidade[i]]["Oorder"])
            varMC.append(i)
            nMCs=nMCs+1




    ind.nMCs=nMCs
    ind.nPMUs_instaladas=ind.nPMUs_instaladas+nPMUSadicionadas
    ind.nMFS_instaladas=ind.nMFS_instaladas+nMFSadicioanadas
    ind.calcula_custo()

    ind.lista_observavel=plano.iloc[plano["Oorder"].sort_values().index]["instalado_mod"].array

    #rotina para gerar plano não critico
    if nMCs>0:
        plano_nc=plano.copy()
        inv_nao_critico=individuo(plano_nc,ind.lista_observavel)
        inv_nao_critico.nMFS_instaladas=ind.nMFS_instaladas
        inv_nao_critico.nPMUs_instaladas=ind.nPMUs_instaladas
    else:
        return



    Permutacao_nc=Permutacao_depois_Observabilidade.copy()


    Hinstaladas_n=np.zeros((len(graph)+1,1))




    
    m=0
    for MC in varMC:
        i=MC    
        if sum(np.abs(Hinstaladas_n[i][:]))>0:
            pass
        else:
            for j in range(startCandidatas,len(plano_nc)):
                if (Hfat2[i][j]>0) & (j<startPMUS):
                    plano_nc.at[Permutacao_nc[j],"instalado_candidatas"]=1
                    plano_nc.at[Permutacao_nc[j],"instalado_mod"]=1
                    inv_nao_critico.nMFS_instaladas=inv_nao_critico.nMFS_instaladas+1
                    Hinstaladas_n=np.concatenate((Hinstaladas_n,np.array(Hfat[:,j]).reshape(-1,1)),axis=1)
                    break
                elif (Hfat2[i][j]>0):
                    plano_nc.at[Permutacao_nc[j],"instalado_candidatas"]=1
                    plano_nc.at[Permutacao_nc[j],"instalado_mod"]=1
                    inv_nao_critico.nPMUs_instaladas=inv_nao_critico.nPMUs_instaladas+1
                    inv_nao_critico.nMFS_instaladas=inv_nao_critico.nMFS_instaladas+1
                    Hinstaladas_n=np.concatenate((Hinstaladas_n,np.array(Hfat2[:,j]).reshape(-1,1)),axis=1)
                    aux=Permutacao_depois_Observabilidade[startPMUS]  
                    Permutacao_depois_Observabilidade[startPMUS]=Permutacao_depois_Observabilidade[j] # coloca ela no fim do arquivo de candidatas
                    Permutacao_depois_Observabilidade[j]=aux # coloca a primeira PMU candidata pro lugar dela
                    Hfat2[:,[j,startPMUS]]=Hfat2[:,[startPMUS,j]] # volta a medida do topo das candidatas para as candidatas
                    startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS
                    offset=startPMUS
                    for k in range(offset,len(plano_nc)): # atualiza todas as medidas daquela PMU para MFS
                        if plano_nc.iloc[Permutacao[k]].i_de==plano_nc.iloc[Permutacao[i]].i_de: # se a medida está na mesma barra
                            plano_nc.at[Permutacao[k],"instalado_candidatas"]=2 # muda o status para candidata
                            Hfat2[:,[k,startPMUS]]=Hfat2[:,[startPMUS,k]]# altera a ordem da Jacobiana
                            aux=Permutacao_depois_Observabilidade[startPMUS]  
                            Permutacao_depois_Observabilidade[startPMUS]=Permutacao_depois_Observabilidade[k] # coloca ela no fim do arquivo de candidatas
                            Permutacao_depois_Observabilidade[k]=aux # coloca a primeira PMU candidata pro lugar dela
                            startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS
                    break
                m=m+1
    if nMCs>0:
        inv_nao_critico.lista_observavel=plano_nc.iloc[plano_nc["Oorder"].sort_values().index]["instalado_mod"].array
        inv_nao_critico.calcula_custo()
    else:
        inv_nao_critico=[]

    return inv_nao_critico

