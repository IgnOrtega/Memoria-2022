def Prob_inclusion_sample_tip_1(Total_nodos, tamaño_sample):  # 17/10/2022
    # Muestreo sin Reemplazo
    #E_i: Es el Evento en donde el nodo i no esta en la muestra
    prob_Ei= ( Total_nodos - tamaño_sample) / Total_nodos
    return 1- prob_Ei

def Prob_inclusion_sample_tip_2(cant_nodos , dic_grados_G , Muestra_general,est_2aristas,cant_aprox_hop_need): # 17/10/2022
    M=cant_aprox_hop_need+1
    cte= (cant_nodos-1)/cant_nodos
    dic_prob_incl={x_i: 1- cte*( (est_2aristas- dic_grados_G[x_i])/est_2aristas )**(M-1) for x_i in Muestra_general }
    return dic_prob_incl

# Muestra_general,muestra_sr, est_2E, cant_nodos=len(list(G.nodes()))

def Prob_inclusion_sample_tip_5(Total_nodos, tamaño_sample):    # 17/10/2022
    # Muestreo sin Reemplazo
    #E_i: Es el Evento en donde el nodo i no esta en la muestra
    prob_Ei= (( Total_nodos - 1) / Total_nodos )**tamaño_sample 
    return 1- prob_Ei

def Prob_sel_sample_tip_5(cant_nodos):   # 17/10/2022
    return 1/cant_nodos

def Prob_sel_sample_tip_678(cant_nodos , dic_grados_G , est_2aristas, Muestra_general):   # 17/10/2022
    M=len(Muestra_general)
    cte= 1/cant_nodos
    dic_prob_sel={x_i: dic_grados_G[x_i]/est_2aristas  for x_i in Muestra_general }
    return cte ,dic_prob_sel