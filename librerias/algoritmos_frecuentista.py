import numpy as np
import random
import pickle 
import pandas as pd
from librerias  import alg_fabrica_grafo as alg_grafo
#import numba as nb
#from librerias import compute_c_ab_sample_7_8_cython  as cython


def N_u_PIMPLE_zip(input_base):
    [N,M,suma_1]=input_base[1:4]
    est_N_u_PIMPLE=N/M *suma_1
    return est_N_u_PIMPLE

def N_u_EMV_zip(input_base):
    N=input_base[1]
    [suma_2,suma_3]=input_base[4:6]
    est_2=N*suma_2/suma_3
    return est_2

def gnsum_zip(input_base): # listo 14/05
    tipo_muestreo=input_base[0]
    
    if tipo_muestreo=="5":
        [N,M,suma_1,suma_2,suma_3,suma_4,suma_5]=input_base[1:-3]
        if suma_4==0:
            return 0
        est_3=N/M*suma_2/suma_5*suma_4

    elif tipo_muestreo == "6" or tipo_muestreo == "7" or tipo_muestreo == "8":
        [tipo_muestreo,N,M]=input_base[0:3]
        [suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8,suma_9,suma_10]=input_base[3:-2]
        if suma_4==0:
            return 0
        
        num=N*suma_1/suma_9
        den=suma_4/suma_10
        est_3=num/den
        
    return est_3

def estimator_RDS_I_zip(input_base,est_C_ab,est_C_ba): # listo 14/05  # solo falta definir est_C_ab y est_C_ba para el segundo if 6,7,8
    [tipo_muestreo]=input_base[0]
    if tipo_muestreo=="5":
        [tipo_muestreo,N,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]=input_base

        if est_C_ab ==0 and est_C_ba==0:
            return 0
        if suma_4==0 or suma_6==0:
            return 0    
        est_D_A=suma_5/suma_4
        est_D_B=suma_7/suma_6
        N_A=N* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)

    elif tipo_muestreo == "6" or tipo_muestreo == "7" or tipo_muestreo == "8":
        [tipo_muestreo,N]=input_base[0:2]
        [suma_4,suma_5,suma_6,suma_7,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[6:]
        if est_C_ab ==0 and est_C_ba==0:
            return 0
        if suma_4==0 or suma_6==0:
            return 0            
        est_D_A= suma_4/suma_10
        est_D_B= suma_6/suma_12

        N_A=N* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)
    return N_A


def estimador_rds2_zip(input_base):  # listo 14/05
    tipo_muestreo=input_base[0]
    if tipo_muestreo=="5":
        [tipo_muestreo,N,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]=input_base
        RDS_2=N/M*suma_4
    elif tipo_muestreo == "6" or tipo_muestreo == "7" or tipo_muestreo == "8":
        [N,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8,suma_9,suma_10]=input_base[1:-2]
        RDS_2=N* suma_10/suma_9
    return RDS_2

def alg_estimador_manager_zip(n_metodo_est,input_base,est_C_ab,est_C_ba):
    L_name_est=["PIMPLE","EMV","GNSUM","RDS_I","RDS_II"]
    metodo_est=L_name_est[n_metodo_est - 1]

    if metodo_est=="PIMPLE":
        est_N_u=N_u_PIMPLE_zip(input_base)
    elif metodo_est=="EMV":       
        est_N_u=N_u_EMV_zip(input_base)
    elif metodo_est=="GNSUM":
        est_N_u=gnsum_zip(input_base)
    elif metodo_est=="RDS_I":
        est_N_u=estimator_RDS_I_zip(input_base,est_C_ab,est_C_ba)
    elif metodo_est=="RDS_II":      
        est_N_u=estimador_rds2_zip(input_base)
    return est_N_u,metodo_est




#@nb.jit(nopython=True)   # no ejecutar nopython de numba (es muy restrictivo)
#@nb.jit(nopython=False)   # no ejecutar nopython de numba (es muy restrictivo)
def calculos_basicos(tipo_muestreo,Muestra_general,dict_nei_G,dic_grados_G,set_infectados,set_noinfectados):
    
    if tipo_muestreo == "5":
        suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8=[0]*8
        num_c_ab,den_c_ab=0,0
        num_c_ba,den_c_ba=0,0
        
        for x_i in Muestra_general:
            neig_i=dict_nei_G[x_i]
            d_i=dic_grados_G[x_i]
            yih= len(neig_i & set_infectados)   #
            yiJ= d_i-yih
            
            func_directriz=len({x_i} & set_infectados) #
            g_directriz=-func_directriz+1

            suma_1=suma_1+yih/d_i
            suma_2=suma_2+yih
            suma_3=suma_3+d_i
            suma_4=suma_4+func_directriz
            suma_5=suma_5+func_directriz*d_i
            suma_6=suma_6+g_directriz
            suma_7=suma_7+g_directriz*d_i
            suma_8=suma_8+yiJ
            
            
            l_u=random.sample(list(neig_i), 1)
            u=l_u[0]
            u_func_directriz=len({u} & set_infectados)
            u_g_directriz=-u_func_directriz+1
            num_c_ab=num_c_ab+d_i*func_directriz*u_g_directriz
            den_c_ab=den_c_ab+d_i*func_directriz
            num_c_ba=num_c_ba+d_i*g_directriz*u_func_directriz
            den_c_ba=den_c_ba+d_i*g_directriz

        if den_c_ab== 0 or den_c_ba==0:
            est_cab,est_cba=0,0
        else:
            est_cab=num_c_ab/den_c_ab
            est_cba=num_c_ba/den_c_ba
        input_base=[suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]
        return input_base,est_cab,est_cba
    
    else: # muestreo tipo 6,7,8 para calcular las sumas.
        
        suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8,suma_9,suma_10,suma_11,suma_12=[0]*12

        for x_i in Muestra_general[1:]:
            neig_i=dict_nei_G[x_i]
            yih= len(neig_i & set_infectados)  
            yiJ= len(neig_i & set_noinfectados)  
            d_i=dic_grados_G[x_i]
            func_directriz=len({x_i} & set_infectados)
            g_directriz=-func_directriz+1

            suma_1=suma_1+yih/d_i
            suma_2=suma_2+yih
            suma_3=suma_3+d_i
            suma_4=suma_4+func_directriz
            suma_5=suma_5+func_directriz*d_i
            suma_6=suma_6+g_directriz
            suma_7=suma_7+g_directriz*d_i
            suma_8=suma_8+yiJ
            suma_9=suma_9+1/d_i
            suma_10=suma_10+func_directriz/d_i
            suma_11=suma_11+yiJ/d_i        
            suma_12=suma_12+g_directriz/d_i
        input_base=[suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8,suma_9,suma_10,suma_11,suma_12]
        
    # calcular C_AB con la muestra de aristas
    if tipo_muestreo == "6":
        num_c_ab,den_c_ab=0,0
        num_c_ba,den_c_ba=0,0
        for i in range(len(Muestra_general[1:-1])):
            x_i=Muestra_general[i]
            func_directriz=len({x_i} & set_infectados)
            g_directriz=-func_directriz+1

            u=Muestra_general[i+1]
            u_func_directriz=len({u} & set_infectados)
            u_g_directriz=-u_func_directriz+1

            num_c_ab=num_c_ab+func_directriz*u_g_directriz
            den_c_ab=den_c_ab+func_directriz
            num_c_ba=num_c_ba+g_directriz*u_func_directriz
            den_c_ba=den_c_ba+g_directriz            

        if den_c_ab== 0 or den_c_ba==0:
            est_cab,est_cba=0,0
        else:
            est_cab=num_c_ab/den_c_ab
            est_cba=num_c_ba/den_c_ba
        
        return input_base,est_cab,est_cba
    else:            
        return input_base
    

def f_estimaciones_grupo(G,dic_inf_y_muestreo,dic_arista_sample,name_df,num_graph,tipo_graph, L_metodo_infec, L_porc_nodos_infec, L_tipo_sample_p_gral, L_porc_nodos_muestra, L_n_metodo_est):
    
    set_nodos=set(G.nodes())
    N=len(set_nodos)
    dict_nei_G={x_i:set(G.neighbors(x_i))  for x_i in list(G.nodes())}    
    dic_grados_G=dict(G.degree( list(G.nodes()) ))

    col_num_graph=[]
    col_tipo_graph=[]
    col_metodo_infec=[]
    col_porc_nodos_infec=[]
    col_tipo_sample_p_gral=[]
    col_porc_nodos_muestra=[]
    col_n_metodo_est=[]

    col_Valor_Estimado=[]
    col_valor_error=[]
    col_valor_error_rel=[]    
    # Opciones variables
    #Inicio  for 
    iteracion=0
    
    for metodo_infec in L_metodo_infec:
        for porc_nodos_infec in L_porc_nodos_infec:
            for tipo_muestreo in L_tipo_sample_p_gral:
                for porc_nodos in L_porc_nodos_muestra:
                    # cargar grupo a estimar H #print(key_N_H) #lista_keys=list(dic_inf_y_muestreo.keys()) #print(lista_keys[0])
                    key_N_H="tipo_subset"+"metodo_infec"+str(metodo_infec)+"porc_nodos_infec" + str(porc_nodos_infec)+"tipo_muestreo"+str(tipo_muestreo)+"porc_nodos" + str(porc_nodos)
                    lista_infectados=dic_inf_y_muestreo[key_N_H]
                    set_infectados=set(lista_infectados)
                    set_noinfectados=set(G.nodes())- set_infectados
                    N_H=len(lista_infectados)
                    
                    # cargar sample nodos
                    key_sample="tipo_sample"+"metodo_infec"+str(metodo_infec)+"porc_nodos_infec" + str(porc_nodos_infec)+"tipo_muestreo"+str(tipo_muestreo)+"porc_nodos" + str(porc_nodos)
                    Muestra_general=dic_inf_y_muestreo[key_sample]
                    M=len(Muestra_general)
                    
                    # cargar sample aristas
                    if tipo_muestreo=="7" or tipo_muestreo=="8":
                        sample_aristas=dic_arista_sample[key_sample]
                    
                    if tipo_muestreo=="5" or tipo_muestreo=="6":
                        input_sumas,est_cab,est_cba=calculos_basicos(tipo_muestreo,Muestra_general,dict_nei_G,dic_grados_G,set_infectados,set_noinfectados) 
                    else:
                        input_sumas=calculos_basicos(tipo_muestreo,Muestra_general,dict_nei_G,dic_grados_G,set_infectados,set_noinfectados) 
                        est_cab,est_cba=compute_c_ab_sample_7_8(sample_aristas,set_infectados)
                        #est_cab,est_cba=cython.compute_c_ab_sample_7_8_cython(sample_aristas,set_infectados)

                    input_base=[tipo_muestreo,N,M]
                    input_base.extend(input_sumas)
                    
                    for n_metodo_est in L_n_metodo_est:  
                                iteracion=iteracion+1
                            
                                Estimacion_alg,metodo_est=alg_estimador_manager_zip(n_metodo_est,input_base,est_cab,est_cba)
                                print(iteracion,metodo_infec,porc_nodos_infec,tipo_muestreo, porc_nodos, n_metodo_est,Estimacion_alg,metodo_est)                                    
                                # N_H = int(N*porc_nodos_infec)
                                error=N_H-Estimacion_alg

                                # guardar informacion de la iteracion en listas correspondientes
                                col_num_graph.append(num_graph)
                                col_tipo_graph.append(tipo_graph)
                                
                                col_metodo_infec.append(metodo_infec)
                                col_porc_nodos_infec.append(porc_nodos_infec)
                                
                                col_tipo_sample_p_gral.append(tipo_muestreo)
                                col_porc_nodos_muestra.append(porc_nodos)
                                
                                col_n_metodo_est.append(metodo_est)
                                col_Valor_Estimado.append(Estimacion_alg)
                                col_valor_error.append(error)
                                col_valor_error_rel.append( error/ N_H)
    # Termino for # Generar DF
    n_col1="Id_grafo"
    n_col2="Tipo_grafo"
    n_col3="Metodo_infeccion"
    n_col4="Porc_infectados"
    n_col5="Tipo_sample_p_gral"
    n_col6="Porc_nodos_muestra"
    n_col9="Num_metodo_est"  
    n_col10="Valor_Estimado"
    n_col11="Valor_error(Exac - Est)"
    n_col12="Valor_error_rel(Exac - Est)"
    
    col1=col_num_graph
    col2=col_tipo_graph
    col3=col_metodo_infec
    col4=col_porc_nodos_infec
    col5=col_tipo_sample_p_gral
    col6=col_porc_nodos_muestra
    col9=col_n_metodo_est
    col10=col_Valor_Estimado
    col11=col_valor_error
    col12=col_valor_error_rel
     
    dic_df={n_col1:col1 ,n_col2:col2 ,n_col3:col3 ,n_col4:col4 ,n_col5:col5 ,n_col6:col6,n_col9:col9,n_col10:col10,n_col11:col11,n_col12:col12}
    Df_est_m1=pd.DataFrame(dic_df)

    # creamos fichero binario 
    fichero_bin= open( name_df, "wb")
    pickle.dump(Df_est_m1, fichero_bin)
    fichero_bin.close()
    return Df_est_m1            

def generate_calculos(densidad,tipo_grafo,n,parametro,inicio, termino):
    # input configuracion
    L_metodo_infec=["1","2","3","4","DP","IP"]
    L_porc_nodos_infec=[0.01,0.1,0.5, 0.9]
    L_tipo_sample_p_gral=["5","6","7","8"]
    L_porc_nodos_muestra=[0.01,0.1,0.5,0.9]    
    L_n_metodo_est=[1,2,3,4,5]
    
    # cargar info de tipo de grafo
    direccion=alg_grafo.nombre_ubi_save(densidad, tipo_grafo,parametro)
    if tipo_grafo=="escala":
        tipo_grafo="algoritmo_escala"
        tipo_grafo_sample="escala"
        alpha="m"+str(int(parametro))        
    else:
        tipo_grafo="algoritmo_aleatorio"
        tipo_grafo_sample="aleatorio"    
        alpha="p"+str(int(1/parametro))    
    
    for num_graph in range(inicio,termino):   

        # cargar grafo
        name_graph ="densidad_"+densidad+"graph_n"+str(n)+alpha+"n°"+str(num_graph)
        file = open(direccion+name_graph, 'rb')
        G= pickle.load(file)
        file.close()
        
        # cargar diccionario de sample de nodos
        if tipo_grafo=="algoritmo_escala":
            alpha="m"+str(int(parametro))        
            alfa=str(int(parametro))
        else:
            alpha="p"+str(int(1/parametro))    
            alfa=str(int(1/parametro))
        directorio="datos_muestra"+direccion[12:]
        name_save_1="subset"+str(tipo_grafo_sample)+"_num_graph"+str(num_graph)+"_alfa"+str(alfa)
        file = open(directorio+name_save_1, 'rb')
        dic_inf_y_muestreo = pickle.load(file) 
        file.close()
        
        # cargar diccionario de sample de aristas
        directorio="datos_muestra"+direccion[12:]
        name_save_1="dic_arista_sample"+str(tipo_grafo_sample)+"_num_graph"+str(num_graph)+"_alfa"+str(alfa)
        file = open(directorio+name_save_1, 'rb')
        dic_arista_sample = pickle.load(file) #close the file
        file.close()

        name_df="data frame/caso_sin_falsos/"+"s_7_8_"+name_graph
        Df_est_m1=f_estimaciones_grupo(G,dic_inf_y_muestreo,dic_arista_sample,name_df,num_graph,tipo_grafo, L_metodo_infec, L_porc_nodos_infec, L_tipo_sample_p_gral, L_porc_nodos_muestra, L_n_metodo_est)
        print(Df_est_m1.head(10))
        #print(direccion+name_graph)

#@nb.jit(nopython=False)   # no ejecutar nopython de numba (es muy restrictivo)
def compute_c_ab_sample_7_8(sample_aristas,set_infectados):
    # calcular valores c_ab
    suma_num_c_ab,suma_den_c_ab=0,0
    suma_num_c_ba,suma_den_c_ba=0,0
    for v,u in sample_aristas:
        v_f_directriz=len({v} & set_infectados)
        v_g_directriz=-v_f_directriz+1

        u_f_directriz=len({u} & set_infectados)
        u_g_directriz=-u_f_directriz+1

        suma_num_c_ab+=v_f_directriz*u_g_directriz
        suma_den_c_ab+=v_f_directriz    
        suma_num_c_ba+=v_g_directriz*u_f_directriz
        suma_den_c_ba+=v_g_directriz  
    if suma_den_c_ab==0 or suma_den_c_ba==0:
        est_cab=0
        est_cba=0
    else:
        est_cab=suma_num_c_ab/suma_den_c_ab
        est_cba=suma_num_c_ba/suma_den_c_ba

    return est_cab,est_cba                    
        




