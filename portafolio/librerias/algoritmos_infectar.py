
#import networkx as nx

from librerias  import alg_sample_v2 as alg_sample
from numpy.random import choice
import numpy as np

def manager_proceso_infec(input_sm):     # 12 / 10 /2022
    G=input_sm[0]
    tipo_muestreo=input_sm[1]
    porc_nodos=input_sm[2]
    list_nei_G=input_sm[3]
    dic_grados_G=input_sm[4]
    prob_acept=input_sm[5]
    cant_acept=input_sm[6]
    
    ## Inicio algoritmo
    cant_nodos_G=len(list(G.nodes()))
    if tipo_muestreo =="1":
        nodos_inf=alg_sample.met_muestra_tip_1(G,porc_nodos) # muestra sin repeticion
    elif tipo_muestreo =="2":
        cant_nodos_dist=int(porc_nodos*cant_nodos_G)
        nodos_muestra,muestra_sr, est_2E, cant_hop=alg_sample.met_muestra_tip_2(G,cant_nodos_dist,list_nei_G,dic_grados_G)
        nodos_inf=muestra_sr                                                    
    elif tipo_muestreo =="3":
        nodos_inf=alg_sample.met_muestra_tip_3(G,porc_nodos,list_nei_G,prob_acept) # muestra sin repeticion
    elif tipo_muestreo =="4":
        nodos_inf=alg_sample.met_muestra_tip_4(G,porc_nodos,list_nei_G,100) # muestra sin repeticion
        nodos_inf=list(set(G.nodes())-set(nodos_inf))
        
    elif tipo_muestreo =="5":   # caso: grado es directamente proporcional a la prob de salir en la muestra
        list_of_candidates=list(G.nodes())
        number_of_items_to_pick= int(porc_nodos*len(list(G.nodes())))
        vec=np.array(list(dic_grados_G.values()))
        total=sum(vec)
        probability_distribution=vec/total
        nodos_inf = choice(list_of_candidates, number_of_items_to_pick, p=probability_distribution,replace=False)

    elif tipo_muestreo =="6":   # caso: 1/grado es directamente proporcional a la prob de salir en la muestra        
                                # caso: grado es inversamente proporcional a la prob de salir en la muestra        
        list_of_candidates=list(G.nodes())
        number_of_items_to_pick= int(porc_nodos*len(list(G.nodes())))
        
        vec=np.array(list(dic_grados_G.values()))
        vec=1/vec
        total=sum(vec)
        probability_distribution=vec/total
        nodos_inf = choice(list_of_candidates, number_of_items_to_pick, p=probability_distribution,replace=False)            
    return nodos_inf



    
    