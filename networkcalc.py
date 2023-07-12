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

    #%%-------------------------------------------inicializa variáveis-----------------------%%#
    nMFSadicioanadas=0
    nPMUscandidatas=0
    nPMUSadicionadas=0

    plano=ind.plano
    startPMUS=ind.startPMUS#util na hora da permutacao
    startCandidatas=ind.startCandidatas#util na hora da permutacao
    tollfill=1e-9

    Hfat=HT.copy()



    # plano.Oorder.array tem a ordem original das medidas colocando instaldas primeiro , primeiro as medidas instaladas, depois as medidas candidatas e depois as PMUs candidatas 
    # o plano já está pré-ordenado só a matriz q está na ordem original

    Hfat=reorder_columns(Hfat,plano.Oorder_DPFMSE) #reordena a matriz jacobiana de acrodo com um vetor permutacao de colunas

    flag=ind.FlagPMUV #se existe PMUs
    Permutacao=np.array(list(plano.index),dtype=int) # vetor de permutacoes tem a mesma ordem da Htriang final e em suas poiscoes tem a posicao da medida no plano
    plano["instalado_mod"]=plano["instalado"]
    # O vetor Permutacao permite acessar de acordo com a coluna i, a medida referente no plano original
    # ele traduz a ordem das colunas com a ordem do plano data frame
    for  i in range(len(graph)-1+flag):
     # Verifica se o pivo é nulo
        if np.abs(Hfat[i][i]) >tolpiv: #checar depois
            # se o pivo é nulo verifica se a medida é instalada
            if plano.iloc[Permutacao[i]].instalado_candidatas!=1:
                #Verifica se o pivo é de uma medida candidata
                if plano.iloc[Permutacao[i]].instalado_candidatas==2: # senão for ele instala ela se for candidata
                    #Se sim ele deve instalar a medida
                    plano.at[Permutacao[i],"instalado_candidatas"]=1 # senão for ele instala ela se for candidata
                    plano.at[Permutacao[i],"instalado_mod"]=1  # senão for ele instala ela se for candidata
                    if plano.iloc[Permutacao[i]].type==5 & flag==0: # se instalou PMU e antes não tinha ele muda a flag pra 1
                        flag=1
                    startCandidatas=startCandidatas+1
                    ## e em seguida instlar a medida par dela
                    mask=(plano.id_MFS==plano.iloc[Permutacao[i]].id_MFS)&(plano.instalado_mod==0) ## encontra posicao do par
                    instalar=int(np.where(Permutacao==plano[mask].index[0])[0])  # medeida que faz parte do mesmo par                        
                    permutaMedida(Hfat,Permutacao,instalar,startCandidatas)
                    plano.at[Permutacao[startCandidatas],"instalado_mod"]=1 # muda o status dela para mod
                    plano.at[Permutacao[startCandidatas],"instalado_candidatas"]=1 #instala a medida
                    startCandidatas=startCandidatas+1
                if plano.iloc[Permutacao[i]].instalado_candidatas==0: # se ela for uma PMU candidata
                    # se a medida do pivo é uma PMU candidata ele instala a medida
                    plano.at[Permutacao[i],"instalado_candidatas"]=1 # se ela for uma PMU candidata ele instala
                    plano.at[Permutacao[i],"instalado_mod"]=1# se ela for uma PMU candidata ele instala
                    if (plano.iloc[Permutacao[i]].type==5) & (flag==0):  # se instalou PMU e antes não tinha ele muda a flag pra 1
                        flag=1
                    startPMUS=startPMUS+1
                    startCandidatas=startCandidatas+1
                    # Em seguida ele instala o par dela
                    mask=(plano.id_MFS==plano.iloc[Permutacao[i]].id_MFS)&(plano.instalado_mod==0) ## encontra posicao do par
                    instalar=int(np.where(Permutacao==plano[mask].index[0])[0])  # medeida que faz parte do mesmo par                        
                    permutaMedida(Hfat,Permutacao,instalar,startCandidatas)
                    plano.at[Permutacao[startCandidatas],"instalado_mod"]=1 # muda o status dela para mod
                    plano.at[Permutacao[startCandidatas],"instalado_candidatas"]=1 #instala a medida
                    startPMUS=startPMUS+1
                    startCandidatas=startCandidatas+1
                    # Depois ele passa as outras na mesma barra para candidatas 
                    for idx,row in plano[mask3].iterrows():
                        plano.at[idx,"instalado_candidatas"]=2
                        mask=(plano.id_MFS==plano.iloc[Permutacao[i]].id_MFS)&(plano.instalado==0) ## encontra posicao do par
                        instalar=int(np.where(Permutacao==idx)[0])
                        permutaMedida(Hfat,Permutacao,instalar,startPMUS)
                        startPMUS=startPMUS+1 
        if np.abs(Hfat[i][i])<tolpiv:
            # se for ele procura a próxima medida com pivo não nulo
            for j in range(i+1,len(plano)):
                if np.abs(Hfat[i][j])>tolpiv: # se o pivo é não nulo
                    permutaMedida(Hfat,Permutacao,i,j) # permuta a medida i com j, j dá a informcao
                    ## Depois daqui o i guarda a posicao da medida permutada e o j da medida que estava na posicao anteriormente
                    if plano.iloc[Permutacao[i]].instalado_candidatas==2: #verifica se a medida permutada era candidata
                        plano.at[Permutacao[i],"instalado_mod"]=1 # muda o status dela para mod
                        plano.at[Permutacao[i],"instalado_candidatas"]=1 #instala a medida
                        ## preciso do contrário 
                        if (plano.iloc[Permutacao[i]].type==5) & (flag==0): #verifica se ela há pmus no plano e se ela for pmus ele coloca a flag para 1
                            flag=1
                        if plano.iloc[Permutacao[j]].instalado_candidatas==1: #coloca a medida que foi retirada para o final das instaladas
                            permutaMedida(Hfat,Permutacao,j,startCandidatas)
                        #instala a medida no topo
                        startCandidatas=startCandidatas+1 # incrementa o identificador de onde começam as MFS_candidatas
                        mask=(plano.id_MFS==plano.iloc[Permutacao[i]].id_MFS)&(plano.instalado_mod==0) ## encontra posicao do par
                        instalar=int(np.where(Permutacao==plano[mask].index[0])[0]) # medeida que faz parte do mesmo par                        
                        permutaMedida(Hfat,Permutacao,instalar,startCandidatas)
                        plano.at[Permutacao[startCandidatas],"instalado_mod"]=1 # muda o status dela para mod
                        plano.at[Permutacao[startCandidatas],"instalado_candidatas"]=1 #instala a medida
                        startCandidatas=startCandidatas+1 # incrementa o identificador de onde começam as MFS_candidatas
                    elif plano.iloc[Permutacao[i]].instalado_candidatas==0: # se precisou de 1 PMU candidata 
                        plano.at[Permutacao[i],"instalado_mod"]=1 #muda o status das medidas
                        plano.at[Permutacao[i],"instalado_candidatas"]=1
                        if (plano.iloc[Permutacao[i]].type==5) & (flag==0):
                            flag=1
                        if (plano.iloc[Permutacao[j]].instalado_candidatas==1) & (startCandidatas < startPMUS): # se a medida é instalada ela tem q ir pro topo    
                            permutaMedida(Hfat,Permutacao,j,startCandidatas) #manda pro topo das candidatas e o topo das candidatas pras PMUs
                            startCandidatas=startCandidatas+1 # incrementa o identificador de onde começam as PMUS
                            permutaMedida(Hfat,Permutacao,j,startPMUS) # manda a medida que está no start candidatas pro start PMUS
                            startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS                        
                        elif plano.iloc[Permutacao[j]].instalado_candidatas==2 & (startCandidatas < startPMUS): # se a medida é instalada ela tem q ir pro topo
                            permutaMedida(Hfat,Permutacao,j,startPMUS)
                            startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS  
                        elif (plano.iloc[Permutacao[j]].instalado_candidatas==1)& (startCandidatas == startPMUS):
                            permutaMedida(Hfat,Permutacao,j,startCandidatas)
                            startCandidatas=startCandidatas+1
                            startPMUS=startPMUS+1

                        #pega medida par e coloca no conjunto das instaladas
                        mask=(plano.id_MFS==plano.iloc[Permutacao[i]].id_MFS)&(plano.instalado_mod==0)
                        instalar=int(np.where(Permutacao==plano[mask].index[0])[0])
                        permutaMedida(Hfat,Permutacao,instalar,startCandidatas)
                        plano.at[Permutacao[startCandidatas],"instalado_mod"]=1 # muda o status dela para mod
                        plano.at[Permutacao[startCandidatas],"instalado_candidatas"]=1 #instala a medida
                        if startCandidatas < startPMUS:
                            startCandidatas=startCandidatas+1 # incrementa o identificador de onde começam as MFS_candidatas
                            permutaMedida(Hfat,Permutacao,instalar,startPMUS)
                            startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS
                        else:
                            startCandidatas=startCandidatas+1 
                            startPMUS=startPMUS+1

                        mask3=(plano["instalado_candidatas"]==0) & (plano.i_de==plano.iloc[Permutacao[i]].i_de)
                        #novo modo de instalar candidatas procurando apenas no pandas, deve dar agilidade ao código
                        for idx,row in plano[mask3].iterrows():
                            plano.at[idx,"instalado_candidatas"]=2
                            mask=(plano.id_MFS==plano.iloc[Permutacao[i]].id_MFS)&(plano.instalado_mod==0) ## encontra posicao do par
                            instalar=int(np.where(Permutacao==idx)[0]) # medeida que faz parte do mesmo par                        
                            permutaMedida(Hfat,Permutacao,instalar,startPMUS)
                            startPMUS=startPMUS+1 
                    break
        Hfat[i,:]=Hfat[i,:]/Hfat[i][i]
        for j in range(i+1,len(graph)):
                if np.abs(Hfat[j][i])>tollfill:
                    Hfat[j,:]=Hfat[j,:]-(Hfat[j][i])*Hfat[i,:]  
                else:  
                    Hfat[j][i]=0  

    #termina de fazer a Hdelta
    for i in range(1,(len(graph)-1+flag)):
        for j in range(0,i):
            if np.abs(Hfat[j][i])>tollfill:
                Hfat[j,:] = Hfat[j,:]-Hfat[j,i]*Hfat[i,:]
            else:
                Hfat[j][i]=0




    plano["order_after_observa"]=0
    
    #%%---------------------Desinstala medidas--------------------------%%#



    ## verifica PMUs instaladas no conjunto básico
    PMUsInstaladas=list(set(plano.iloc[Permutacao[0:len(graph)-1+flag]].i_de.tolist()))
    ## Ids das medidas basicas

    IDF_MFSbasicos=list(set(plano.iloc[Permutacao[0:len(graph)-1+flag]].id_MFS.tolist()))
    
    ## Medida existente




    plano.iloc[Permutacao[len(graph)-1+flag:]].existente.isin(IDF_MFSbasicos)
    
    #
    maskPMUs=plano.iloc[Permutacao[len(graph)-1+flag:]].i_de.isin(PMUsInstaladas) #Se a medida não está no conjunto básico mas está em uma barra com PMU
    maskIDs=plano.iloc[Permutacao[len(graph)-1+flag:]].id_MFS.isin(IDF_MFSbasicos) # Se a medida é par de uma que foi instalada
    maskExisente=plano.existente==1 # se a medida é existente

    mask1= (maskPMUs) & (~maskIDs)& (~maskExisente) # se ela pertence a uma PMU, não faz par com uma básica e não é existente
    mask1=(~(plano.de !=-999) )| mask1 # ela é considera candidata e não instalada

    plano.loc[mask1,"instalado_candidatas"]=2 # todas essas medidas são candidatas
    plano.loc[mask1,"instalado_mod"]=0 #e não estão instaladas

    ind.nPMUs_instaladas=len(PMUsInstaladas) #novo número de PMUs instaladas
    
   
    

    maskPMUs=~plano.iloc[Permutacao[len(graph)-1+flag:]].i_de.isin(PMUsInstaladas) #SE a medida não está no conjunto básico nem a barra dela está com PMUs
    maskIDs=plano.iloc[Permutacao[len(graph)-1+flag:]].id_MFS.isin(IDF_MFSbasicos) #Se a medida é par de uma básica
    maskExisente=plano.existente==1 # se a medida é existente
    
    mask2= (maskPMUs) & ~(maskIDs) & (~maskExisente)  # se ela não pertence a uma PMU, não faz par com uma básica e não é existente 
    mask2=(~(plano.de !=-999) )| mask2 # ela será desinstalada

    plano.loc[mask2,"instalado_candidatas"]=0 # ela passa a ser PMU candidata
    plano.loc[mask2,"instalado_mod"]=0 #desisntala todas


    #reordena o plano 
    plano.loc[Permutacao,"order_after_observa"]=list(range(len(Permutacao))) #ordem da permutacao

    #reordena o plano de acordo com a permutacao
 
    sort_order = {1: 0, 2: 1, 0: 2} #atribui novamente o sort_order 
    plano["sort_order_observabilidade"]=plano["instalado_candidatas"].map(sort_order) # ordena o DMED 
    
    
    plano2=plano.sort_values(by=["sort_order_observabilidade","order_after_observa"],kind="stable",ignore_index=0)# retordena o plano de acordo com a prioridade de maneira "stable" ou seja, não altera a ordem se não precisar
    
    
    Permutacao_depois_Observabilidade=plano2.order_after_observa.to_list()

    del plano2

    # aqui está errado o vetor permutacao depois da observabilidade 
    # esse vetor pode ser utilizado para reordenar a matriz e 
    # o próprio vetor plano para fazer a referencia a medida original

    Hfat2=reorder_columns(Hfat,Permutacao_depois_Observabilidade) # reordena a matriz de acordo com o permutacao
    Permutacao=np.array(Permutacao) 
    Permutacao_depois_Observabilidade=np.array(Permutacao[Permutacao_depois_Observabilidade]) # reordena o permutacao para corresponder a ordem inicial do plano, ou seja o coneteúdo da posicao i faz referencia ao indice do elemento no plano

    # Permutacao_depois_Observabilidade tem a ordem dos indices originaais do plano



    ##----------------Desinstala medidas ok-----------------


    #preenche o número de medidas candidatas

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
            lstMC.append(plano.iloc[Permutacao_depois_Observabilidade[i]]["Oorder_DPFMSE"])
            varMC.append(i)
            nMCs=nMCs+1


    ## Resultados de saída das MCs
    ind.nMCs=nMCs
    ind.calcula_custo()

    
    
    ind.lista_observavel=plano.iloc[plano["Oorder_DPFMSE"].sort_values().index]["instalado_mod"].array

    #rotina para gerar plano não critico
    if nMCs>0:
        plano_nc=plano.copy()
        inv_nao_critico=individuo(plano_nc,ind.lista_observavel)
    else:
        return



    Permutacao_nc=Permutacao_depois_Observabilidade.copy()


    Hinstaladas_n=np.zeros((len(graph)+1,1))




    ###
    m=0
    for MC in varMC:
        i=MC    
        if sum(np.abs(Hinstaladas_n[i][:]))>0:
            pass
        else:
            for j in range(startCandidatas,len(plano_nc)):
                if (Hfat2[i][j]>0) & (j<startPMUS):
                    plano_nc.at[Permutacao_nc[j],"instalado_candidatas"]=1 # instala a medida no plano se a medida for candidata
                    plano_nc.at[Permutacao_nc[j],"instalado_mod"]=1 # instala a medida no plano se a medida for candidata
                    Hinstaladas_n=np.concatenate((Hinstaladas_n,np.array(Hfat[:,j]).reshape(-1,1)),axis=1)
                    ##instala o par dela
                    mask=(plano_nc.id_MFS==plano_nc.iloc[Permutacao_nc[j]].id_MFS)&(plano_nc.instalado_mod==0) ## encontra posicao do par
                    instalar=int(np.where(Permutacao_nc==plano_nc[mask].index[0])[0]) # medeida que faz parte do mesmo par                        
                    plano_nc.at[Permutacao_nc[instalar],"instalado_mod"]=1 # muda o status dela para mod
                    plano_nc.at[Permutacao_nc[instalar],"instalado_candidatas"]=1 #instala a medida
                    Hinstaladas_n=np.concatenate((Hinstaladas_n,np.array(Hfat[:,instalar]).reshape(-1,1)),axis=1)
                    break
                elif (Hfat2[i][j]>0):
                    plano_nc.at[Permutacao_nc[j],"instalado_candidatas"]=1 # se for PMU candidata ele instala
                    plano_nc.at[Permutacao_nc[j],"instalado_mod"]=1       

                    Hinstaladas_n=np.concatenate((Hinstaladas_n,np.array(Hfat2[:,j]).reshape(-1,1)),axis=1)
                    
                    mask=(plano_nc.id_MFS==plano_nc.iloc[Permutacao_nc[j]].id_MFS)&(plano_nc.instalado_mod==0) ## encontra posicao do par
                    permutaMedida(Hfat2,Permutacao_nc,j,startPMUS)
                    startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS
                    instalar=int(np.where(Permutacao_nc==plano_nc[mask].index[0])[0]) # medeida que faz parte do mesmo par                        
                    plano_nc.at[Permutacao_nc[instalar],"instalado_mod"]=1 # muda o status dela para mod
                    plano_nc.at[Permutacao_nc[instalar],"instalado_candidatas"]=1 #instala a medida
                    Hinstaladas_n=np.concatenate((Hinstaladas_n,np.array(Hfat2[:,instalar]).reshape(-1,1)),axis=1)
                    permutaMedida(Hfat2,Permutacao_nc,instalar,startPMUS)
                    startPMUS=startPMUS+1
                    #coloca nas medidas basicas, mas eu sempre verifico na matriz Hinstaladas, se uma medida instalada previamente já resolveu
                    offset=startPMUS
                    for k in range(offset,len(plano_nc)): # atualiza todas as medidas daquela PMU para MFS
                        if plano_nc.iloc[Permutacao_nc[k]].i_de==plano_nc.iloc[Permutacao_nc[i]].i_de: # se a medida está na mesma barra
                            plano_nc.at[Permutacao_nc[k],"instalado_candidatas"]=2 # muda o status para candidata
                            permutaMedida(Hfat2,Permutacao_nc,k,startPMUS)
                            startPMUS=startPMUS+1 # incrementa o identificador de onde começam as PMUS
                    break
                m=m+1
    if nMCs>0:
        inv_nao_critico.lista_observavel=plano_nc.iloc[plano_nc["Oorder_DPFMSE"].sort_values().index]["instalado_mod"].array
        inv_nao_critico.calcula_custo()
    else:
        inv_nao_critico=[]

    return inv_nao_critico



def permutaMedida(Hfat,Permutacao,x,y):
    """
    Troca posicao de x com y na matriz e no vetor 
    """
    
    Hfat[:,[x,y]]=Hfat[:,[y,x]]# altera a ordem da Jacobiana
    aux=Permutacao[y]  
    Permutacao[y]=Permutacao[x] # coloca ela no fim do arquivo de candidatas
    Permutacao[x]=aux # coloca a primeira PMU candidata pro lugar dela
