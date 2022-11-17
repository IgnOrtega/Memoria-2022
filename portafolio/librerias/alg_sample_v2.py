import random
import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import bernoulli
import pickle 


# met_muestra_tip_1
# compute_hop_recomendado_v2
# met_muestra_tip_2_v2
# add_modo_random_no_aislado
# proceso_infeccion_met_2
# met_muestra_tip_3
# proceso_infeccion_met_1
# met_muestra_tip_4
# met_muestra_tip_5
# met_muestra_tip_6

# inclusion_prob_manager
# Prob_inclusion_sample_tip_2
# Prob_inclusion_sample_tip_1
# sample_manager
# add_modo_random_no_aislado
# sample_infec_pobla

def met_muestra_tip_1(G,porc_nodos_muestra): # 11 -10 -2022  Uniforme sin repetición
    
    nodos_G=list(G.nodes())
    cant_nodos=len(nodos_G)
    cant_nodos_muestra=int(porc_nodos_muestra * cant_nodos)
    nodos_muestra=random.sample(nodos_G, cant_nodos_muestra)
    return nodos_muestra


def compute_hop_recomendado_v2(G,cant_nodos_dist,Muestra_general,dic_grados_G):  # 12 /10 /2022
    #L cantidad de saltos 
    
    # Estimacion de la distribución de grados Poblacional y n_k
    # -------------------------------------------- # --------------------------------------------     
    
    col1=[dic_grados_G[x_i] for x_i in Muestra_general]
    n_col1="Grados nodos G"
    dic_df={n_col1:col1}
    df_grados_G=pd.DataFrame(dic_df)
    vector_k=np.array( df_grados_G[n_col1].value_counts().keys().tolist()   ) # lista de grados en el mismo orden d_i
    vector_nk=np.array( df_grados_G[n_col1].value_counts().tolist()    )  # cant de repeticiones de  cierto grado= f(d_i) 
    max_nk=max(vector_nk)

    # calcular denominador
    dim=vector_nk.shape
    dim=dim[0]
    den=0
    suma=0
    for i in range(dim):
        x=vector_k[i]
        q_x=vector_nk[i]
        suma=suma+q_x/x

    func_grado_pob=[0]*dim
    for i in range(dim):
        q_x=vector_nk[i]    
        x=vector_k[i]    
        func_grado_pob[i]=q_x/x
    func_grado_pob=np.array(func_grado_pob)    
    func_grado_pob=func_grado_pob/suma
    
    cant_nodos=len(list(G.nodes))    
    vector_nk=func_grado_pob*cant_nodos
    
    # -------------------------------------------- # -------------------------------------------- 
    vector_p_k= vector_nk/cant_nodos
    barra_k=sum( np.multiply(vector_k,vector_p_k ) )
    vector_P_k_salto= ( 1/barra_k )* np.multiply(vector_k, vector_p_k )
    
    L=cant_nodos_dist -1
    # L=0  # cant_nodos_muestra=0
    if L>= 0:
        vec_v_l_k_c=vector_p_k  # 
        vec_v_l_k_a=vec_v_l_k_c
        cant_hop_recomendado=L
    if L>=1:    
        vec_v_l_k_c=vector_p_k + vector_P_k_salto
        vec_v_l_k_b=vec_v_l_k_c
        cant_hop_recomendado=L
    if L>=2:
        vec_alpha_k= (1 - 1 / barra_k ) * vector_P_k_salto    
        for l in range(2,L+1): # notar que que si saltamos la cant "cant_nodos_dist" entonces la cant de nodos distintos hasta aqui, es menor o igual a la variable "cant_nodos_dist"
            vec_beta_k =1- np.multiply(1/vector_nk , vec_v_l_k_a ) 
            # nuevos vectores
            vec_v_l_k_c =vec_v_l_k_b + np.multiply(vec_alpha_k,vec_beta_k)  # c=2
            # actualizar valores a,b
            vec_v_l_k_a = vec_v_l_k_b  # 1=a
            vec_v_l_k_b = vec_v_l_k_c  # b=2
        V_l=sum(vec_v_l_k_c)

        while V_l < cant_nodos_dist :
            L=L+1
            vec_beta_k =1- np.multiply(1/vector_nk , vec_v_l_k_a ) 
            # nuevos vectores
            vec_v_l_k_c =vec_v_l_k_b + np.multiply(vec_alpha_k,vec_beta_k)  # c=2
            # actualizar valores a,b
            vec_v_l_k_a = vec_v_l_k_b  # 1=a
            vec_v_l_k_b = vec_v_l_k_c  # b=2 
            V_l=sum(vec_v_l_k_c)
        cant_hop_recomendado=L
    return cant_hop_recomendado, barra_k 

def met_muestra_tip_2(G,cant_nodos_dist,list_nei_G,dic_grados_G): # RW  # 11 -10 -2022 RW sin repetición
    nodos_G=list(G.nodes())
    cant_nodos=len(nodos_G)
    cant_saltos=cant_nodos_dist-1
    
    # Agregar un nodo inicial adecuado
    # revisar si el nodo es aislado
    while True:
        
        start_nodo=random.sample(nodos_G, 1);
        start_nodo=start_nodo[0];
        x_i=start_nodo;
        
        if dic_grados_G[x_i]!=0:
            break
    nodos_muestra=[]            
    # Agregar un nodo x1,x2,...,        
    nodos_muestra.append(x_i)
    
    for i in range(cant_saltos):  # ya se agrego un nodo a la muestra
        # elegimos un vecino del nodo x_i que se convertirá en el nodo x_i
        neig_x_i=list_nei_G[x_i]        
        # si i+1 es menor que el largo de la muestra necesaria, el algoritmo sigue
        
        x_i=random.sample(neig_x_i, 1)
        x_i=x_i[0]
        nodos_muestra.append(x_i)

        
    cant_aprox_hop_need,barra_k=compute_hop_recomendado_v2(G,cant_nodos_dist,nodos_muestra,dic_grados_G)
    hop_faltantes=cant_aprox_hop_need-cant_saltos
    
    for i in range(hop_faltantes):
        # elegimos un vecino del nodo x_i que se convertirá en el nodo x_i
        neig_x_i=list_nei_G[x_i]
        # si i+1 es menor que el largo de la muestra necesaria, el algoritmo sigue
        
        x_i=random.sample(neig_x_i, 1)
        x_i=x_i[0]
        nodos_muestra.append(x_i)
        
    nodos_muestra=list(nodos_muestra)
    cant_aprox_hop_need_2,barra_k=compute_hop_recomendado_v2(G,cant_nodos_dist,nodos_muestra,dic_grados_G)
    est_2E= barra_k*cant_nodos
    muestra_sr=list(set(nodos_muestra))

    return nodos_muestra,muestra_sr, est_2E,cant_aprox_hop_need

def add_modo_random_no_aislado(G, nodos_G):
    # Agregar un nodo inicial adecuado
    # revisar si el nodo es aislado
    while True:
        
        start_nodo=random.sample(nodos_G, 1);
        start_nodo=start_nodo[0];
        x_i=start_nodo;
        
        if len(list(G.neighbors(x_i)))!=0:
            break
    return x_i

def proceso_infeccion_met_2(G,list_aceptados,last_list_aceptados,list_nei_G,prob_acept):
    from scipy.stats import bernoulli
    # resetear lista de pendientes:    
    last_list_aceptados_new=[]
    for x_i in last_list_aceptados:
        #calcular vecinos de ultimos nodos aceptados
        nei_xi=list_nei_G[x_i]
        d_xi=len(nei_xi)
        
        eva_list_pend=bernoulli.rvs(prob_acept,size= d_xi)
        for ii in range(d_xi):
            if eva_list_pend[ii]==1:
                last_list_aceptados_new.append(nei_xi[ii])
        
    # resetear last_list_aceptados
    last_list_aceptados=random.sample(last_list_aceptados_new,len(last_list_aceptados_new))
    #actualizar lista de nodos aceptados
    list_aceptados=list_aceptados+last_list_aceptados
    return list_aceptados, last_list_aceptados


def met_muestra_tip_3(G,porc_nodos,list_nei_G,prob_acept):  # 05 /11 /2022  # sin repeticion  RDS probabilistico sin repeticion
    nodos_G=list(G.nodes())
    cant_nodos=len(nodos_G)
    cant_nodos_pedida=int(porc_nodos * cant_nodos)
    list_nei_G_v2=[]
    for x_i in nodos_G:
        list_nei_G_v2.append(list(list_nei_G[x_i] ) )

    lista=[]
    larg_lista_actual=0
    last_list_aceptados=[]
    delta_=+1
    
    x_i=add_modo_random_no_aislado(G, nodos_G) # obtener nodos inicial aleatorio
    last_list_aceptados.append(x_i)
    lista.append(x_i)
    set_3=set(lista)
    
    while 0 < delta_:
        set_1=set_3 # set_1 es la muestra sin la ultima oleada
        #Si quedan elementos pendientes:
        # Utilizar la lista con los nodos pendientes es decir no hacer nada
        # calcular lista de infectados producidos por x_i del tamaño requerido
        lista, last_list_aceptados = proceso_infeccion_met_2(G,lista,last_list_aceptados ,list_nei_G_v2,prob_acept)  # hacer el procedimiento x50
        # ver si lista cumple con el tamaño pedido
        set_3=set(lista) # set_3 es la muestra con la ultima oleada
        larg_lista_actual=len(set_3) 
        
        # calcular cantidad de nodos requeridos
        delta_=cant_nodos_pedida-larg_lista_actual
        #print(cant_nodos_pedida,larg_lista_actual,delta_)
        
        # if # no se cumple largo pedido se crea una lista nueva
        # if # Si se cumple largo pedido se sale, corta la lista y se devuelve el resultado
        
    set_delta = set_3- set_1 # los ultimos nodos agregados
    # print("set_1",len(set_1),"set_3",len(set_3), "set_delta",len(set_delta),"delta_",delta_)
    delta_=cant_nodos_pedida-len(set_1)
    lista=list(set_1)+  random.sample(set_delta,delta_)
    return lista
def proceso_infeccion_met_1(G,list_aceptados,last_list_aceptados,list_nei_G,cant_acept): # 12 /10 /2022
    # resetear lista de pendientes:    
    last_pendientes=last_list_aceptados
    last_list_aceptados=[]
    for x_i in last_pendientes:
        #calcular vecinos de ultimos nodos aceptados
        
        last_list_aceptados=random.sample(list_nei_G[x_i],cant_acept)+last_list_aceptados
    last_list_aceptados=random.sample(last_list_aceptados,len(last_list_aceptados))
    #actualizar lista de nodos aceptados
    list_aceptados=list_aceptados+last_list_aceptados
    return list_aceptados, last_list_aceptados

def met_muestra_tip_4(G,porc_nodos,list_nei_G,cant_acept):  # 11 /10 /2022  RDS fijo sin repeticion
    nodos_G=list(G.nodes())
    cant_nodos=len(nodos_G)
    cant_nodos_pedida=int(porc_nodos * cant_nodos)
    lista=[]
    larg_lista_actual=0
    last_list_aceptados=[]
    delta_=+1
    
    x_i=add_modo_random_no_aislado(G, nodos_G) # obtener nodos inicial aleatorio
    last_list_aceptados.append(x_i)
    lista.append(x_i)
    set_3=set(lista)
    
    while 0 < delta_:
        set_1=set_3
        #Si quedan elementos pendientes:
        # Utilizar la lista con los nodos pendientes es decir no hacer nada
        # calcular lista de infectados producidos por x_i del tamaño requerido
        lista, last_list_aceptados = proceso_infeccion_met_1(G,lista,last_list_aceptados,list_nei_G,cant_acept)
        # ver si lista cumple con el tamaño pedido
        set_3=set(lista)
        larg_lista_actual=len(set_3) 
        
        # calcular cantidad de nodos requeridos
        delta_=cant_nodos_pedida-larg_lista_actual
        
        # if # no se cumple largo pedido se crea una lista nueva
        # if # Si se cumple largo pedido se sale, corta la lista y se devuelve el resultado
    set_delta = set_3- set_1 # los ultimos nodos agregados
    # print("set_1",len(set_1),"set_3",len(set_3), "set_delta",len(set_delta),"delta_",delta_)
    delta_=cant_nodos_pedida-len(set_1)
    lista=list(set_1)+  random.sample(set_delta,delta_)
    return lista

def met_muestra_tip_7(G,porc_nodos,list_nei_G,prob_acept):  # 11 /10 /2022  RDS probabilistico con repeticion
    nodos_G=list(G.nodes())
    cant_nodos=len(nodos_G)
    cant_nodos_pedida=int(porc_nodos * cant_nodos)

    list_nei_G_v2=[]
    for x_i in nodos_G:
        list_nei_G_v2.append(list(list_nei_G[x_i] ) )
        

    lista=[]
    larg_lista_actual=0
    last_list_aceptados=[]
    delta_=+1
    
    x_i=add_modo_random_no_aislado(G, nodos_G) # obtener nodos inicial aleatorio
    last_list_aceptados.append(x_i)
    lista.append(x_i)
    
    while 0 < delta_:
        #Si quedan elementos pendientes:
        # Utilizar la lista con los nodos pendientes es decir no hacer nada
        # calcular lista de infectados producidos por x_i del tamaño requerido
        lista, last_list_aceptados = proceso_infeccion_met_2(G,lista,last_list_aceptados ,list_nei_G_v2,prob_acept)  # hacer el procedimiento x50
        # ver si lista cumple con el tamaño pedido
        larg_lista_actual=len(lista) 
        
        # calcular cantidad de nodos requeridos
        delta_=cant_nodos_pedida-larg_lista_actual
        
        # if # no se cumple largo pedido se crea una lista nueva
        # if # Si se cumple largo pedido se sale, corta la lista y se devuelve el resultado

    return lista[0:cant_nodos_pedida]

def met_muestra_tip_8(G,porc_nodos,list_nei_G,cant_acept):  # 11 /10 /2022  RDS fijo con repeticion
    nodos_G=list(G.nodes())
    cant_nodos=len(nodos_G)
    cant_nodos_pedida=int(porc_nodos * cant_nodos)
    lista=[]
    larg_lista_actual=0
    last_list_aceptados=[]
    delta_=+1
    
    x_i=add_modo_random_no_aislado(G, nodos_G) # obtener nodos inicial aleatorio
    last_list_aceptados.append(x_i)
    lista.append(x_i)
    
    while 0 < delta_:
        #Si quedan elementos pendientes:
        # Utilizar la lista con los nodos pendientes es decir no hacer nada
        # calcular lista de infectados producidos por x_i del tamaño requerido
        
        lista, last_list_aceptados = proceso_infeccion_met_1(G,lista,last_list_aceptados,list_nei_G,cant_acept)
        # ver si lista cumple con el tamaño pedido
        larg_lista_actual=len(lista) 
        
        # calcular cantidad de nodos requeridos
        delta_=cant_nodos_pedida-larg_lista_actual
        
        # if # no se cumple largo pedido se crea una lista nueva
        # if # Si se cumple largo pedido se sale, corta la lista y se devuelve el resultado

    return lista[0:cant_nodos_pedida]

def met_muestra_tip_5(G,porc_nodos_muestra): # 11 -10 -2022  Uniforme con  repeticion
    
    nodos_G=list(G.nodes())
    cant_nodos=len(nodos_G)
    cant_nodos_muestra=int(porc_nodos_muestra * cant_nodos)
    
    nodos_muestra=random.choices(nodos_G, k=cant_nodos_muestra)
    return nodos_muestra

def sample_manager(input_sm):     # 17 / 10 /2022
    
    G=input_sm[0]
    tipo_muestreo=input_sm[1]
    porc_nodos=input_sm[2]
    list_nei_G=input_sm[3]
    dic_grados_G=input_sm[4]

    ## Inicio algoritmo
    cant_nodos_G=len(list(G.nodes()))
    if tipo_muestreo =="1":
        output_sample_manager=met_muestra_tip_1(G,porc_nodos)
    elif tipo_muestreo =="2":
        cant_nodos_dist=int(porc_nodos*cant_nodos_G)
        output_sample_manager=met_muestra_tip_2(G,cant_nodos_dist,list_nei_G,dic_grados_G)
        output_sample_manager=[output_sample_manager[1],output_sample_manager[2],output_sample_manager[3]]        
    elif tipo_muestreo =="3":
        prob_acept=input_sm[5]        
        output_sample_manager=met_muestra_tip_3(G,porc_nodos,list_nei_G,prob_acept)
    elif tipo_muestreo =="4":
        cant_acept=input_sm[6]                
        output_sample_manager=met_muestra_tip_4(G,porc_nodos,list_nei_G,cant_acept)
    elif tipo_muestreo =="5":
        output_sample_manager=met_muestra_tip_5(G,porc_nodos)
    elif tipo_muestreo =="6":
        cant_nodos_dist=int(porc_nodos*cant_nodos_G)
        output_sample_manager=met_muestra_tip_2(G,cant_nodos_dist,list_nei_G,dic_grados_G)
        output_sample_manager=[output_sample_manager[0],output_sample_manager[2],output_sample_manager[3]]
    elif tipo_muestreo =="7":
        prob_acept=input_sm[5]                
        output_sample_manager=met_muestra_tip_7(G,porc_nodos,list_nei_G,prob_acept)
    elif tipo_muestreo =="8":
        cant_acept=input_sm[6]        
        output_sample_manager=met_muestra_tip_8(G,porc_nodos,list_nei_G,cant_acept)
    return output_sample_manager


# Muestreo de poblacion infectada o no infectada
def sample_infec(input_base, input_muestra):
    set_nodos_G=input_base[0]
    set_infectados=input_base[1]
    porc_nodos_inf=input_base[2]    
    
    tipo_muestreo=input_muestra[0]
    porc_nodos=input_muestra[1]
    prob_acept=input_muestra[2]

    if porc_nodos_inf <= 0.5:

        set_noinfec= set_nodos_G - set_infectados
        H = G.subgraph(set_noinfec)
    else:
        H = G.subgraph(set_infectados)
    dict_nei_H={x_i:set(H.neighbors(x_i))  for x_i in list(H.nodes())}    

    nodes_H=list(H.nodes()) 
    cant_nodos_H=len(nodes_H)
    if tipo_muestreo =="3":
        output_sample_manager=met_muestra_tip_3(H,porc_nodos,dict_nei_H,prob_acept)
    elif tipo_muestreo =="7":
        output_sample_manager=met_muestra_tip_7(H,porc_nodos,dict_nei_H,prob_acept)    
    Muestra_general=output_sample_manager    
    return Muestra_general


