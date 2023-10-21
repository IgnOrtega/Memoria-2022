import random
from scipy.stats import bernoulli
import pickle 
from librerias  import alg_fabrica_grafo as alg_grafo
from datetime import datetime
import numpy as np
from numpy.random import choice
from librerias import proceso_infeccion_met_1_3_cython as pvp_cython
#from librerias import proceso_infeccion_met_1_3_random_cython as pvp_cython_random

def met_muestra_tip_1(G,porc_nodos_muestra): # 11 -10 -2022  Uniforme sin repetición
    nodos_G=list(G.nodes())
    cant_nodos=len(nodos_G)
    cant_nodos_muestra=int(porc_nodos_muestra * cant_nodos)
    nodos_muestra=random.sample(nodos_G, cant_nodos_muestra)
    return nodos_muestra


def met_muestra_tip_2(nodos_G,cant_nodos,porc_nodos_muestra,list_nei_G,dic_grados_G): # RW  # 11 -10 -2022 RW sin repetición
    cant_nodos_dist=int(porc_nodos_muestra * cant_nodos)
    cant_saltos_x=2*cant_nodos_dist+200
    
    # Agregar un nodo inicial adecuado
    # revisar si el nodo es aislado
    while True:
        start_nodo=random.sample(nodos_G, 1)
        start_nodo=start_nodo[0]
        x_i=start_nodo
        
        if dic_grados_G[x_i]!=0:
            break
    nodos_muestra=[]            
    # Agregar un nodo x1,x2,...,        
    nodos_muestra.append(x_i)
    
    while True:
        for i in range(cant_saltos_x):  # ya se agrego un nodo a la muestra
            # elegimos un vecino del nodo x_i que se convertirá en el nodo x_i
            neig_x_i=list_nei_G[x_i]        
            # si i+1 es menor que el largo de la muestra necesaria, el algoritmo sigue
            x_i=random.sample(neig_x_i, 1)
            x_i=x_i[0]
            nodos_muestra.append(x_i)
        set_nodos_muestra=set(nodos_muestra)
        delta=cant_nodos_dist-len(set_nodos_muestra)
        if delta>0:
            cant_saltos_x=2*delta + 200
        else:
            break
    nodos_muestra= list(set_nodos_muestra)
    return nodos_muestra[0:cant_nodos_dist] 


def proceso_infeccion_met_1_3(last_list_aceptados ,list_nei_G,dic_grados_G,prob_acept):
    # resetear lista de pendientes:
    last_list_aceptados_new=[]
    for x_i in last_list_aceptados:
        #calcular vecinos de ultimos nodos aceptados
        nei_xi=list_nei_G[x_i]
        d_xi=dic_grados_G[x_i]
        eva_list_pend=bernoulli.rvs(prob_acept,size= d_xi)
        for ii in range(d_xi):
            if eva_list_pend[ii]==1:
                #print(nei_xi[ii])
                last_list_aceptados_new.append(nei_xi[ii])
    
    return last_list_aceptados_new

 
def proceso_infeccion_met_1_3_random(last_list_aceptados ,list_nei_G,prob_acept):
    # resetear lista de pendientes:
    last_list_aceptados_new=[]
    for x_i in last_list_aceptados:
        #calcular vecinos de ultimos nodos aceptados
        nei_xi=list_nei_G[x_i]
        for u in nei_xi:
            if random.random() < prob_acept:
                last_list_aceptados_new.append(u)
    return last_list_aceptados_new    

def met_muestra_tip_3(nodos_G,cant_nodos,list_nei_G,dic_grados_G,prob_acept,porc_nodos):
    cant_nodos_pedida=int(cant_nodos*porc_nodos)
    lista=[]
    last_list_aceptados=[]
    
    while True:
        start_nodo=random.sample(nodos_G, 1)
        start_nodo=start_nodo[0]
        x_i=start_nodo
        if dic_grados_G[x_i]!=0:
            break
    delta_=+1        
    last_list_aceptados.append(x_i)
    lista.append(x_i)
    
    while 0 < delta_: #Si quedan elementos pendientes:
        
        # calcular lista de infectados producidos por x_i del tamaño requerido
        #last_list_aceptados_aux = proceso_infeccion_met_1_3(last_list_aceptados ,list_nei_G,dic_grados_G,prob_acept) 
        last_list_aceptados_aux = pvp_cython.proceso_infeccion_met_1_3(last_list_aceptados ,list_nei_G,dic_grados_G,prob_acept)  
        #last_list_aceptados_aux =proceso_infeccion_met_1_3_random(last_list_aceptados ,list_nei_G,prob_acept)
        #last_list_aceptados_aux =pvp_cython_random.proceso_infeccion_met_1_3_random(last_list_aceptados ,list_nei_G,prob_acept)
        
        if len(last_list_aceptados_aux)== 0: # caso en que todos los postulantantes sean rechazados.
            while True:
                start_nodo=random.sample(nodos_G, 1)
                start_nodo=start_nodo[0]
                x_i=start_nodo
                if dic_grados_G[x_i]!=0:
                    break
            last_list_aceptados_aux=[x_i]

        last_list_aceptados=random.sample(last_list_aceptados_aux,len(last_list_aceptados_aux))
        lista.extend(last_list_aceptados)

        set_3=set(lista) # set_3 es la muestra con la ultima oleada
        larg_lista_actual=len(set_3) 

        # calcular cantidad de nodos requeridos
        delta_=cant_nodos_pedida-larg_lista_actual
        
    lista=list(set_3)[0:cant_nodos_pedida]
    return lista



def proceso_infeccion_met_1_4(last_list_aceptados ,list_nei_G,dic_grados_G,cant_acept):
    # resetear lista de pendientes:
    last_list_aceptados_new=[]
    for x_i in last_list_aceptados:
        #calcular vecinos de ultimos nodos aceptados
        nei_i=list_nei_G[x_i]
        d_i=dic_grados_G[x_i]
        if d_i>cant_acept:
            nodos_elegidos=random.sample(nei_i,cant_acept)
        else:
            nodos_elegidos=nei_i
        last_list_aceptados_new.extend(nodos_elegidos)
    return last_list_aceptados_new




def met_muestra_tip_4(nodos_G,cant_nodos,list_nei_G,dic_grados_G,cant_acept,porc_nodos):
    cant_nodos_pedida=int(cant_nodos*porc_nodos)
    lista=[]
    last_list_aceptados=[]
    
    while True:
        start_nodo=random.sample(nodos_G, 1)
        start_nodo=start_nodo[0]
        x_i=start_nodo
        if dic_grados_G[x_i]!=0:
            break
    delta_=+1        
    last_list_aceptados.append(x_i)
    lista.append(x_i)
    suma=0
    while 0 < delta_: #Si quedan elementos pendientes:
        
        # calcular lista de infectados producidos por x_i del tamaño requerido
        last_list_aceptados_aux = proceso_infeccion_met_1_4(last_list_aceptados ,list_nei_G,dic_grados_G,cant_acept)
        last_list_aceptados=random.sample(last_list_aceptados_aux,len(last_list_aceptados_aux))
        lista.extend(last_list_aceptados)

        set_3=set(lista) # set_3 es la muestra con la ultima oleada
        larg_lista_actual=len(set_3) 

        # calcular cantidad de nodos requeridos
        delta_=cant_nodos_pedida-larg_lista_actual
    lista=list(set_3)[0:cant_nodos_pedida]
    return lista



def met_muestra_tip_5(G,porc_nodos_muestra): # 11 -10 -2022  Uniforme con  repeticion
    
    nodos_G=list(G.nodes())
    cant_nodos=len(nodos_G)
    cant_nodos_muestra=int(porc_nodos_muestra * cant_nodos)
    
    nodos_muestra=random.choices(nodos_G, k=cant_nodos_muestra)
    return nodos_muestra

def met_muestra_tip_6(nodos_G,cant_nodos,porc_nodos_muestra,list_nei_G,dic_grados_G): 
    cant_nodos_dist=int(porc_nodos_muestra * cant_nodos)
    cant_saltos_x=cant_nodos_dist-1
    while True:
        start_nodo=random.sample(nodos_G, 1)
        start_nodo=start_nodo[0]
        x_i=start_nodo
        if dic_grados_G[x_i]!=0:
            break
    nodos_muestra=[]            
    # Agregar un nodo x1,x2,...,        
    nodos_muestra.append(x_i)
    
    for i in range(cant_saltos_x):  # ya se agrego un nodo a la muestra
        # elegimos un vecino del nodo x_i que se convertirá en el nodo x_i
        neig_x_i=list_nei_G[x_i]        
        # si i+1 es menor que el largo de la muestra necesaria, el algoritmo sigue

        x_i=random.sample(neig_x_i, 1)
        x_i=x_i[0]
        nodos_muestra.append(x_i)
    return nodos_muestra

def proceso_infeccion_met_2_7(last_list_aceptados ,list_nei_G,dic_grados_G,prob_acept):
    # resetear lista de pendientes:
    last_list_aceptados_new=[]
    sample_arista_add=[]
    
    for x_i in last_list_aceptados:
        #calcular vecinos de ultimos nodos aceptados
        nei_xi=list_nei_G[x_i]
        d_xi=dic_grados_G[x_i]
        eva_list_pend=bernoulli.rvs(prob_acept,size= d_xi)
        for ii in range(d_xi):
            if eva_list_pend[ii]==1:
                u=nei_xi[ii]
                last_list_aceptados_new.append(u)
                sample_arista_add.append((x_i,u)) 
    
    return last_list_aceptados_new,sample_arista_add




# Nota:
# esto lista_a=random.sample(lista_a,len(lista_a)) es más lento que "lista_b=random.sample(lista_a,len(lista_a))", yo creo que sucede por que estamos operando 
# y sacando elementos de un lugar que van cambiando valores, entonces esto complica el reemplazo. 

def met_muestra_tip_7(nodos_G,cant_nodos,list_nei_G,dic_grados_G,prob_acept,porc_nodos):
    cant_nodos_pedida=int(cant_nodos*porc_nodos)
    lista=[]
    last_list_aceptados=[]
    
    while True:
        start_nodo=random.sample(nodos_G, 1)
        start_nodo=start_nodo[0]
        x_i=start_nodo
        if dic_grados_G[x_i]!=0:
            break
    delta_=+1        
    last_list_aceptados.append(x_i)
    lista.append(x_i)
    sample_arista=[]
    while 0 < delta_: #Si quedan elementos pendientes:
        
        # calcular lista de infectados producidos por x_i del tamaño requerido
        last_list_aceptados_aux,sample_arista_add = proceso_infeccion_met_2_7(last_list_aceptados ,list_nei_G,dic_grados_G,prob_acept)
        if len(last_list_aceptados_aux)== 0: # caso en que todos los postulantantes sean rechazados.
            while True:
                start_nodo=random.sample(nodos_G, 1)
                start_nodo=start_nodo[0]
                x_i=start_nodo
                if dic_grados_G[x_i]!=0:
                    break
            last_list_aceptados_aux=[x_i]
        else:
            delta_=cant_nodos_pedida-len(lista) - len(last_list_aceptados_aux)
            if delta_> 0:
                last_list_aceptados=last_list_aceptados_aux
                lista.extend(last_list_aceptados)
                sample_arista.extend(sample_arista_add)
                
            
    len_add=len(last_list_aceptados_aux)
    #print(cant_nodos_pedida,len_add,len(lista))
    random_selection=random.sample(list(range(0,len_add)),cant_nodos_pedida-len(lista))
    for ii in random_selection:
        lista.append(last_list_aceptados_aux[ii])
        sample_arista.append(sample_arista_add[ii])
        
    return lista,sample_arista

def proceso_infeccion_met_2_8(last_list_aceptados ,list_nei_G,dic_grados_G,cant_acept):
    # resetear lista de pendientes:
    last_list_aceptados_new=[]
    sample_arista_add=[]
    
    for x_i in last_list_aceptados:
        #calcular vecinos de ultimos nodos aceptados
        nei_i=list_nei_G[x_i]
        d_i=dic_grados_G[x_i]
        if d_i>cant_acept:
            nodos_elegidos=random.sample(nei_i,cant_acept)
        else:
            nodos_elegidos=nei_i
        sample_arista_add.extend(list(zip([x_i]*d_i,nodos_elegidos)))
        last_list_aceptados_new.extend(nodos_elegidos)
    return last_list_aceptados_new,sample_arista_add


def met_muestra_tip_8(nodos_G,cant_nodos,list_nei_G,dic_grados_G,cant_acept,porc_nodos):
    cant_nodos_pedida=int(cant_nodos*porc_nodos)
    lista=[]
    last_list_aceptados=[]
    
    while True:
        start_nodo=random.sample(nodos_G, 1)
        start_nodo=start_nodo[0]
        x_i=start_nodo
        if dic_grados_G[x_i]!=0:
            break
    delta_=+1        
    last_list_aceptados.append(x_i)
    lista.append(x_i)
    sample_arista=[]
    while 0 < delta_: #Si quedan elementos pendientes:
        # calcular lista de infectados producidos por x_i del tamaño requerido
        last_list_aceptados_aux,sample_arista_add = proceso_infeccion_met_2_8(last_list_aceptados ,list_nei_G,dic_grados_G,cant_acept)
        delta_=cant_nodos_pedida-len(lista) - len(last_list_aceptados_aux)
        
        if delta_> 0:
            last_list_aceptados=last_list_aceptados_aux
            lista.extend(last_list_aceptados)
            sample_arista.extend(sample_arista_add)
            
    len_add=len(last_list_aceptados_aux)
    random_selection=random.sample(list(range(0,len_add)),cant_nodos_pedida-len(lista))
    for ii in random_selection:
        lista.append(last_list_aceptados_aux[ii])
        sample_arista.append(sample_arista_add[ii])
        
    return lista,sample_arista



def generate_samples(densidad,tipo_grafo,n,parametro,inicio, termino,porc ):

    # input comple
    L_metodo_infec=["1","2","3","4","DP","IP"]
    L_porc_nodos_infec=[0.01,0.1,0.5, 0.9]

    L_tipo_sample_p_gral=["5","6","7","8"]
    L_porc_nodos_muestra=[0.01,0.1,0.5,0.9]

    direccion=alg_grafo.nombre_ubi_save(densidad, tipo_grafo,parametro)
    if tipo_grafo=="escala":
        alpha_1="m"
        alpha=alpha_1+str(int(parametro))        
        alfa=str(int(parametro))
    else:
        alpha_1="p"
        alpha=alpha_1+str(int(1/parametro))    
        alfa=str(int(1/parametro))


    for num_graph in range(inicio,termino):   
        date_1=datetime.now()    
        # cargar grafo
        name_graph ="densidad_"+densidad+"graph_n"+str(n)+alpha+"n°"+str(num_graph)
        file = open(direccion+name_graph, 'rb')
        G= pickle.load(file)
        file.close()
        # crear sample de nodos y aristas
        dic_lista_infec_y_sample,dic_arista_sample=f_dict_sample_nodes_edges(G,L_metodo_infec, L_porc_nodos_infec, L_tipo_sample_p_gral, L_porc_nodos_muestra,porc)

        # guardar sample de nodos
        directorio="datos_muestra"+direccion[12:]
        name_save_1="subset"+str(tipo_grafo)+"_num_graph"+str(num_graph)+"_alfa"+str(alfa)
        fichero_bin= open( directorio+name_save_1, "wb")
        pickle.dump(dic_lista_infec_y_sample, fichero_bin)
        fichero_bin.close()    

        # guardar sample de aristas
        directorio="datos_muestra"+direccion[12:]
        name_save_1="dic_arista_sample"+str(tipo_grafo)+"_num_graph"+str(num_graph)+"_alfa"+str(alfa)
        fichero_bin= open( directorio+name_save_1, "wb")
        pickle.dump(dic_arista_sample, fichero_bin)
        fichero_bin.close()    

        date_2=datetime.now()    
        print(num_graph,date_1,date_2, date_2-date_1)


def f_dict_sample_nodes_edges(G, L_metodo_infec, L_porc_nodos_infec, L_tipo_sample_p_gral, L_porc_nodos_muestra,porc):
    nodos_G=list(G.nodes())    
    cant_nodos=len(nodos_G)    
    dic_grados_G=dict(G.degree( list(G.nodes()) ))
    lista_degree=list(dic_grados_G.values())

    list_nei_G=[]
    for x_i in nodos_G:
        list_nei_G.append(list(G.neighbors(x_i)) )       
        
    list_nei_G.append(nodos_G)
    list_nei_G.append(cant_nodos)
    #  variable:
    grado_mean=np.mean(lista_degree)
    prob_acept_p_inf=0.5
    cant_acept_p_inf= int(prob_acept_p_inf*grado_mean)

    #cant_personas= ?
    lista_valores_c=[]
    for porc_nodos in L_porc_nodos_muestra:
        tamaño_sample=int(porc_nodos*cant_nodos)
        c=2
        cant_olas=6
        tamaño_sample_c=find_valor_c(c, cant_olas)
        while tamaño_sample_c< tamaño_sample:
            c+=1
            tamaño_sample_c=find_valor_c(c, cant_olas)
        #print(c,tamaño_sample_c)    
        lista_valores_c.append(c)
    dict_valores_c={L_porc_nodos_muestra[i]:lista_valores_c[i] for i in range(len(lista_valores_c))}
    
    # fijo 
    lista_infec_y_samplekey,lista_infec_y_sample=[],[]
    lista_arista_sample_key,lista_arista_sample=[],[]
    #Inicio  for 
    for metodo_infec in L_metodo_infec:
        for porc_nodos_infec in L_porc_nodos_infec:
            for tipo_muestreo in L_tipo_sample_p_gral:
                for porc_nodos in L_porc_nodos_muestra:
                    print("metodo_infec",metodo_infec,"porc_nodos_infec",porc_nodos_infec,"tipo_muestreo",tipo_muestreo,"porc_nodos",porc_nodos)
                    
                    cant_acept_p_inf=dict_valores_c[porc_nodos_infec]    
                    prob_acept_p_inf=cant_acept_p_inf/grado_mean      
                    if prob_acept_p_inf> 1:
                        prob_acept_p_inf=1                                  
                    ## Infectar nodos:
                    if metodo_infec =="1":
                            lista_infectados=met_muestra_tip_1(G,porc_nodos_infec) # muestra sin repeticion
                    elif metodo_infec =="2":# list_nei_G: lista de listas
                        lista_infectados=met_muestra_tip_2(nodos_G,cant_nodos,porc_nodos_infec,list_nei_G,dic_grados_G)
                    elif metodo_infec =="3":# list_nei_G: lista de listas
                        lista_infectados=met_muestra_tip_3(nodos_G,cant_nodos,list_nei_G,dic_grados_G,prob_acept_p_inf,porc_nodos_infec) # muestra sin repeticion
                    elif metodo_infec =="4": # list_nei_G: lista de listas
                        lista_infectados=met_muestra_tip_4(nodos_G,cant_nodos,list_nei_G,dic_grados_G,cant_acept_p_inf,porc_nodos_infec) # muestra sin repeticion
                    elif metodo_infec =="DP":   # caso: grado es directamente proporcional a la prob de salir en la muestra
                        number_of_items_to_pick= int(porc_nodos_infec*cant_nodos)
                        vec=np.array(lista_degree)
                        total=sum(vec)
                        probability_distribution=vec/total
                        lista_infectados = choice(nodos_G, number_of_items_to_pick, p=probability_distribution,replace=False)
                    elif metodo_infec =="IP":   # caso: 1/grado es directamente proporcional a la prob de salir en la muestra        
                                                # caso: grado es inversamente proporcional a la prob de salir en la muestra        
                        number_of_items_to_pick= int(porc_nodos_infec*cant_nodos)
                        vec=np.array(lista_degree)
                        vec=1/vec
                        total=sum(vec)
                        probability_distribution=vec/total
                        lista_infectados = choice(nodos_G, number_of_items_to_pick, p=probability_distribution,replace=False)            

                    key_N_H="tipo_subset"+"metodo_infec"+str(metodo_infec)+"porc_nodos_infec" + str(porc_nodos_infec)+"tipo_muestreo"+str(tipo_muestreo)+"porc_nodos" + str(porc_nodos)
                    lista_infec_y_samplekey.append(key_N_H)
                    lista_infec_y_sample.append(lista_infectados)                    
                    
                    # Muestreo de poblacion general

                    cant_acept=dict_valores_c[porc_nodos]    
                    prob_acept=cant_acept/grado_mean
                    if prob_acept> 1:
                        prob_acept=1
                    if tipo_muestreo =="5":
                        Muestra_general=met_muestra_tip_5(G,porc_nodos)
                    elif tipo_muestreo =="6":
                        Muestra_general=met_muestra_tip_6(nodos_G,cant_nodos,porc_nodos,list_nei_G,dic_grados_G) # list_nei_G: lista de listas
                    else:
                        if tipo_muestreo =="7":
                            Muestra_general,sample_arista=met_muestra_tip_7(nodos_G,cant_nodos,list_nei_G,dic_grados_G,prob_acept,porc_nodos)
                        else:
                            Muestra_general,sample_arista=met_muestra_tip_8(nodos_G,cant_nodos,list_nei_G,dic_grados_G,cant_acept,porc_nodos)
                        
                        key_sample_arista="tipo_sample"+"metodo_infec"+str(metodo_infec)+"porc_nodos_infec" + str(porc_nodos_infec)+"tipo_muestreo"+str(tipo_muestreo)+"porc_nodos" + str(porc_nodos)
                        lista_arista_sample_key.append(key_sample_arista)
                        lista_arista_sample.append(sample_arista)

                    key_sample="tipo_sample"+"metodo_infec"+str(metodo_infec)+"porc_nodos_infec" + str(porc_nodos_infec)+"tipo_muestreo"+str(tipo_muestreo)+"porc_nodos" + str(porc_nodos)
                    lista_infec_y_samplekey.append(key_sample)
                    lista_infec_y_sample.append(Muestra_general)
                    
    dic_lista_infec_y_sample={lista_infec_y_samplekey[i]:lista_infec_y_sample[i] for i in range(len(lista_infec_y_samplekey))}
    dic_arista_sample={lista_arista_sample_key[i]:lista_arista_sample[i] for i in range(len(lista_arista_sample_key))}    

    return dic_lista_infec_y_sample,dic_arista_sample    




def find_valor_c(c, cant_olas):
    # c personas, Recordar que la ola 1 se elimina
    cant_pos_x=c   # cant nodos ola x=2
    total_x=cant_pos_x # cant total hasta la ola x=2

    for num_ola in range(3,cant_olas+1):
        cant_pos_x*=c
        total_x+=cant_pos_x
    return total_x

