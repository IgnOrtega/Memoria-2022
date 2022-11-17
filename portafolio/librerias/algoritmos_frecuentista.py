import numpy as np
from librerias  import alg_sample_v2 as alg_sample
from librerias.alg_prob_incl_y_sel  import *
import random
import pickle 


def gnsum_seed_zip(input_base):
    
    [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
    [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
    [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
    
    est_2aristas=(M-1) * t/(suma_9)
    if suma_4==0:
        est_3=0
        return est_3
    est_3=1/(M-1)*  est_2aristas*suma_1*suma_10/suma_4
        
    return est_3


def gnsum_adapt_seed_zip(input_base):
    [porc_nodos_inf,tipo_muestreo]=input_base[:2]
    
    [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
    [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
    [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
    
    
    if porc_nodos_inf >=0.5:
        est_4=gnsum_zip(input_base)
    else:
        if suma_6==0:
            return 0
        est_2aristas=(M-1) * t/suma_9
        est_4=t-1/(M-1)*  est_2aristas*suma_11*suma_12/suma_6
        
    return est_4


def estimator_RDS_I_V1_seed_zip(input_base,input_base_2,Muestra_general):
    [tipo_muestreo]=input_base[1]

    [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
    [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
    [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
    
    est_2aristas=(M-1) * t/(suma_9)
    est_C_ab,est_C_ba=estimacion_Cab(input_base_2,Muestra_general)
    
    if est_C_ab ==0 and est_C_ba==0:
        return 0    
    
    if suma_4==0 or suma_6==0:
            return 0
    est_D_A= suma_4 / suma_10
    est_D_B= suma_6 / suma_12
    
    N_A=t* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)
    return N_A

def estimador_rds1_v2_seed_zip(input_base,input_base_2,Muestra_general):  # estima D_B
    [tipo_muestreo]=input_base[1]

    [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
    [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
    [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
    
    est_C_ab,est_C_ba=estimacion_Cab(input_base_2,Muestra_general)
    if est_C_ab ==0 and est_C_ba==0:
        return 0
    
    est_2aristas=(M-1) * t/(suma_9)
    est_D_B= suma_6 / suma_12
    if suma_6==0:   
        return 0
    beta=est_C_ba*est_D_B
    est_D_A= (est_2aristas * beta) /( t*(beta +est_D_B *est_C_ab )- est_2aristas*est_C_ab )    
    est_rds_1v2=t* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)

    return est_rds_1v2


def estimador_rds1_v2_adapt_seed_zip(input_base,input_base_2,Muestra_general):
    porc_nodos_inf,tipo_muestreo=input_base[:2]
    
    [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
    [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
    [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
        
    if porc_nodos_inf <= 0.5:
        # caso muestreo 6,7,8 y porc bajo infec
        est_rds_1v2_adapt=estimador_rds1_v2_zip(input_base,input_base_2,Muestra_general)

    else:
        # caso muestreo 6,7,8 y porc alto infec
        est_C_ab,est_C_ba=estimacion_Cab(input_base_2,Muestra_general)

        if est_C_ab ==0 and est_C_ba==0:
            return 0    

        est_2aristas=(M-1) * t/(suma_9)
        if suma_4==0:   
            return 0        
        est_D_A= suma_4 / suma_10

        beta=est_C_ab*est_D_A
        est_D_B= (est_2aristas * beta) /( t*(beta +est_D_A *est_C_ba )- est_2aristas*est_C_ba )
        est_rds_1v2_adapt=t* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)            
            
    return est_rds_1v2_adapt


def estimador_rds2_seed_zip(input_base):
    tipo_muestreo=input_base[1]

    [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
    [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
    [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
    
    est_2aristas=(M-1) * t/(suma_9)
    RDS_2=1/(M-1)*est_2aristas*suma_10
    return RDS_2

def estimador_rds2_adapt_seed_zip(input_base):
    porc_nodos_inf,tipo_muestreo=input_base[:2]
    
    [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
    [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
    [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
    
    if porc_nodos_inf >= 0.5: # caso muestreo 6,7,8 y porc alto infec
        RDS_2_adapt=estimador_rds2_zip(input_base)
    else: # caso muestreo 6,7,8 y porc bajo infec
        est_2aristas=(M-1) * t/(suma_9)
        RDS_2_adapt=t- 1/(M-1)*est_2aristas*suma_12
    return RDS_2_adapt






def alg_estimador_manager_zip(n_metodo_est,input_base,input_base_2,Muestra_general):
    est_N_B=-10    
    # definicion del diccionario con los algs est
    L1=["GNSUM_WS","GNSUM_ADAPT_WS","RDS_I_v1_WS","RDS_I_v2_WS","RDS_I_v2_ADAPT_WS","RDS_II_WS","RDS_IIv2_ADAPT_WS"]
    L_name_est=["PIMPLE","EMV","GNSUM","GNSUM_ADAPT","RDS_I_v1","RDS_I_v2","RDS_I_v2_ADAPT","RDS_II","RDS_II_ADAPT"] +L1
  

    metodo_est=L_name_est[n_metodo_est - 1]
    if metodo_est=="PIMPLE":
        est_N_u=N_u_PIMPLE_zip(input_base)
    elif metodo_est=="EMV":       
        est_N_u=N_u_EMV_zip(input_base)
    elif metodo_est=="GNSUM":
        est_N_u=gnsum_zip(input_base)
    elif metodo_est=="GNSUM_ADAPT":
         est_N_u=gnsum_adapt_zip(input_base)
    elif metodo_est=="RDS_I_v1":
        est_N_u=estimator_RDS_I_V1_zip(input_base,input_base_2,Muestra_general)
    elif metodo_est=="RDS_I_v2":
        est_N_u=estimador_rds1_v2_zip(input_base,input_base_2,Muestra_general)
    elif metodo_est=="RDS_I_v2_ADAPT":
        est_N_u=estimador_rds1_v2_adapt_zip(input_base,input_base_2,Muestra_general)
    elif metodo_est=="RDS_II":      
        est_N_u=estimador_rds2_zip(input_base)
    elif metodo_est=="RDS_II_ADAPT":     
        est_N_u=estimador_rds2_adapt_zip(input_base)
    
    elif metodo_est=="GNSUM_WS":
        est_N_u=gnsum_seed_zip(input_base)
    elif metodo_est=="GNSUM_ADAPT_WS":
         est_N_u=gnsum_adapt_seed_zip(input_base)
    elif metodo_est=="RDS_I_v1_WS":
        est_N_u=estimator_RDS_I_V1_seed_zip(input_base,input_base_2,Muestra_general)
    elif metodo_est=="RDS_I_v2_WS":
        est_N_u=estimador_rds1_v2_seed_zip(input_base,input_base_2,Muestra_general)
    elif metodo_est=="RDS_I_v2_ADAPT_WS":
        est_N_u=estimador_rds1_v2_adapt_seed_zip(input_base,input_base_2,Muestra_general)
    elif metodo_est=="RDS_II_WS":      
        est_N_u=estimador_rds2_seed_zip(input_base)
    elif metodo_est=="RDS_IIv2_ADAPT_WS":     
        est_N_u=estimador_rds2_adapt_seed_zip(input_base)    
    return est_N_u,metodo_est

def calculos_basicos(tipo_muestreo,Muestra_general,input_base_2):
    [t,set_nodos,dict_nei_G,dic_grados_G,porc_nodos_inf,set_infectados]=input_base_2
    M=len(Muestra_general)
    
    if tipo_muestreo == "1" or tipo_muestreo == "2" or tipo_muestreo == "3" or tipo_muestreo == "4" or tipo_muestreo == "5":

        set_noinfectados=set_nodos-set_infectados
        suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8=[0]*8
        for x_i in Muestra_general:
            neig_i=dict_nei_G[x_i]
            d_i=dic_grados_G[x_i]
            yih= len(neig_i & set_infectados)  
            yiJ= d_i-yih
            
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
        input_base=[porc_nodos_inf,tipo_muestreo,t,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]
    
    #return suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,M
        return input_base
    
    else:

        set_noinfectados=set_nodos-set_infectados
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

        x_0=Muestra_general[0]
        neig_0=dict_nei_G[x_0]
        yih= len(neig_0 & set_infectados)  
        d_i=dic_grados_G[x_0]

        # infectados
        suma_01=yih/d_i
        suma_02=yih
        suma_03=d_i

        func_directriz=len({x_0} & set_infectados)
        suma_04=func_directriz
        suma_05=func_directriz*d_i

        # NO infectados
        g_directriz=-func_directriz+1
        suma_06=g_directriz
        suma_07=g_directriz*d_i
        yiJ= len(neig_0 & set_noinfectados)
        suma_08=yiJ
        input_base=[porc_nodos_inf,tipo_muestreo,t,M]+[suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5,suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]
        
        return input_base
    
def N_u_PIMPLE_zip(input_base):
    tipo_muestreo=input_base[1]
    if tipo_muestreo == "1" or tipo_muestreo == "2" or tipo_muestreo == "3" or tipo_muestreo == "4" or tipo_muestreo == "5":
        [porc_nodos_inf,tipo_muestreo,t,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]=input_base
        est_N_u_PIMPLE=t/M *suma_1
    else:
        [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
        [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
        [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
        est_N_u_PIMPLE=t/M *(suma_01 +suma_1)
    return est_N_u_PIMPLE



def N_u_EMV_zip(input_base):
    tipo_muestreo=input_base[1]
    if tipo_muestreo == "1" or tipo_muestreo == "2" or tipo_muestreo == "3" or tipo_muestreo == "4" or tipo_muestreo == "5":
        [porc_nodos_inf,tipo_muestreo,t,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]=input_base
        est_2=t*suma_2/suma_3
    else:
        [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
        [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
        [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
        est_2=t*(suma_02+suma_2)/(suma_03+suma_3)
        
    return est_2


def gnsum_zip(input_base):
    tipo_muestreo=input_base[1]
    
    if tipo_muestreo=="5":
        [porc_nodos_inf,tipo_muestreo,t,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]=input_base
        if suma_4==0:
            return 0
        est_3=t/M*suma_2/suma_5*suma_4
    elif tipo_muestreo == "6" or tipo_muestreo == "7" or tipo_muestreo == "8":
        [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
        [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
        [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
        if suma_4==0:
            return 0
        num=1/M*(t * suma_02 + (M-1) * t/(suma_9) *suma_1)
        den=(suma_05*t + (M-1) * t/(suma_9) * suma_4 ) / (suma_04 * t +  (M-1) * t/suma_9 * suma_10)
        est_3=num/den
        
    return est_3





def gnsum_adapt_zip(input_base):
    [porc_nodos_inf,tipo_muestreo]=input_base[:2]
    
    if tipo_muestreo=="5":
        [porc_nodos_inf,tipo_muestreo,t,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]=input_base
        if porc_nodos_inf >=0.5:
            est_4=gnsum_zip(input_base)
        else:
            if suma_6==0:
                return 0
            est_B=t/M*suma_8/suma_7*suma_6
            est_4=t-est_B
    elif tipo_muestreo == "6" or tipo_muestreo == "7" or tipo_muestreo == "8":
        [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
        [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
        [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
        if porc_nodos_inf >=0.5:
            est_4=gnsum_zip(input_base)
        else:
            if suma_6==0:
                return 0
            num=1/M*(suma_08*t  + (M-1) * t/(suma_9) * suma_11)
            den=(suma_07*t + (M-1) * t/(suma_9) * suma_6 ) / (suma_06 * t +  (M-1) * t/suma_9 * suma_12)
            est_B=num/den
            est_4=t-est_B
    return est_4


def estimator_RDS_I_V1_zip(input_base,input_base_2,Muestra_general):
    [tipo_muestreo]=input_base[1]
    if tipo_muestreo=="5":
        [porc_nodos_inf,tipo_muestreo,t,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]=input_base

        est_C_ab,est_C_ba=estimacion_Cab(input_base_2,Muestra_general)
        if est_C_ab ==0 and est_C_ba==0:
            return 0
        if suma_4==0 or suma_6==0:
            return 0    
        est_D_A=suma_5/suma_4
        est_D_B=suma_7/suma_6

        N_A=t* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)
    elif tipo_muestreo == "6" or tipo_muestreo == "7" or tipo_muestreo == "8":
        [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
        [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
        [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
        
        est_2aristas=(M-1) * t/(suma_9)
        est_C_ab,est_C_ba=estimacion_Cab(input_base_2,Muestra_general)
        if est_C_ab ==0 and est_C_ba==0:
            return 0
        if suma_4==0 or suma_6==0:
            return 0            
        est_D_A= (suma_05*t+ est_2aristas*suma_4)/(suma_04*t+ est_2aristas *suma_10 )
        est_D_B= (suma_07*t+ est_2aristas*suma_6)/(suma_06*t+ est_2aristas *suma_12 )

        N_A=t* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)
    return N_A

def estimador_rds1_v2_zip(input_base,input_base_2,Muestra_general):  # estima D_B
    [tipo_muestreo]=input_base[1]
    if tipo_muestreo=="5":
        [porc_nodos_inf,tipo_muestreo,t,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]=input_base

        est_C_ab,est_C_ba=estimacion_Cab(input_base_2,Muestra_general)
        if est_C_ab ==0 and est_C_ba==0:
            return 0
        
        if suma_6 ==0:
            return 0        
        est_D_B=suma_7/suma_6
        est_2aristas=t*suma_3/M
        
        beta=est_C_ba*est_D_B
        est_D_A= (est_2aristas * beta) /( t*(beta +est_D_B *est_C_ab )- est_2aristas*est_C_ab )    
        est_rds_1v2=t* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)
    elif tipo_muestreo == "6" or tipo_muestreo == "7" or tipo_muestreo == "8":
        
        [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
        [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
        [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
        est_2aristas=(M-1) * t/(suma_9)
        est_C_ab,est_C_ba=estimacion_Cab(input_base_2,Muestra_general)
        if est_C_ab ==0 and est_C_ba==0:
            return 0        
        if suma_6 ==0:
            return 0         
        est_D_B= (suma_07*t+ est_2aristas*suma_6)/(suma_06*t+ est_2aristas *suma_12 )
        
        beta=est_C_ba*est_D_B
        est_D_A= (est_2aristas * beta) /( t*(beta +est_D_B *est_C_ab )- est_2aristas*est_C_ab )    
        est_rds_1v2=t* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)

    return est_rds_1v2


def estimador_rds1_v2_adapt_zip(input_base,input_base_2,Muestra_general):
    porc_nodos_inf,tipo_muestreo=input_base[:2]
    
    if tipo_muestreo=="5":
        [porc_nodos_inf,tipo_muestreo,t,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]=input_base
        if porc_nodos_inf <= 0.5:
            # caso muestreo 5 y porc bajo infec
            est_rds_1v2_adapt=estimador_rds1_v2_zip(input_base,input_base_2,Muestra_general)

        else:
            # caso muestreo 5 y porc alto infec
            est_C_ab,est_C_ba=estimacion_Cab(input_base_2,Muestra_general)
            if est_C_ab ==0 and est_C_ba==0:
                return 0
            
            if suma_4 ==0:
                return 0
            est_D_A=suma_5/suma_4
            est_2aristas=t*suma_3/M

            beta=est_C_ab*est_D_A
            est_D_B= (est_2aristas * beta) /( t*(beta +est_D_A *est_C_ba )- est_2aristas*est_C_ba )
            est_rds_1v2_adapt=t* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)
                
    elif tipo_muestreo == "6" or tipo_muestreo == "7" or tipo_muestreo == "8": # caso muestreo 6,7,8
                
        [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
        [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
        [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
        
        if porc_nodos_inf <= 0.5:
            # caso muestreo 6,7,8 y porc bajo infec
            est_rds_1v2_adapt=estimador_rds1_v2_zip(input_base,input_base_2,Muestra_general)
            
        else:
            # caso muestreo 6,7,8 y porc alto infec
            est_C_ab,est_C_ba=estimacion_Cab(input_base_2,Muestra_general)
            if est_C_ab ==0 and est_C_ba==0:
                return 0            
            est_2aristas=(M-1) * t/(suma_9)
            if suma_4 ==0:
                return 0            
            est_D_A= (suma_05*t+ est_2aristas*suma_4)/(suma_04*t+ est_2aristas *suma_10 )
            
            beta=est_C_ab*est_D_A
            est_D_B= (est_2aristas * beta) /( t*(beta +est_D_A *est_C_ba )- est_2aristas*est_C_ba )
            est_rds_1v2_adapt=t* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)            
            
    return est_rds_1v2_adapt


def estimador_rds2_zip(input_base):
    tipo_muestreo=input_base[1]
    if tipo_muestreo=="5":
        [porc_nodos_inf,tipo_muestreo,t,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]=input_base
        RDS_2=t/M*suma_4
    elif tipo_muestreo == "6" or tipo_muestreo == "7" or tipo_muestreo == "8":
                
        [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
        [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
        [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
        est_2aristas=(M-1) * t/(suma_9)
        RDS_2=1/M*(suma_04*t +est_2aristas*suma_10)
    return RDS_2

def estimador_rds2_adapt_zip(input_base):
    porc_nodos_inf,tipo_muestreo=input_base[:2]
    if tipo_muestreo=="5":
        [porc_nodos_inf,tipo_muestreo,t,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]=input_base
        if porc_nodos_inf >= 0.5: # caso muestreo 5 y porc alto infec
            RDS_2_adapt=t/M*suma_4
        else: # caso muestreo 5 y porc bajo infec
            est_NJ=t/M*suma_6
            RDS_2_adapt=t-est_NJ
    elif tipo_muestreo == "6" or tipo_muestreo == "7" or tipo_muestreo == "8":
                
        [porc_nodos_inf,tipo_muestreo,t,M]=input_base[0:4]
        [suma_01,suma_1,suma_02,suma_2,suma_03,suma_3,suma_04,suma_4,suma_05,suma_5]=input_base[4:14]
        [suma_06,suma_6,suma_07,suma_7,suma_08,suma_8,suma_9,suma_10,suma_11,suma_12]=input_base[14:24]
        if porc_nodos_inf >= 0.5: # caso muestreo 6,7,8 y porc alto infec
            RDS_2_adapt=estimador_rds2_zip(input_base)
        else: # caso muestreo 6,7,8 y porc bajo infec
            est_2aristas=(M-1) * t/(suma_9)
            RDS_2_adapt=t- 1/M*(suma_06*t +est_2aristas*suma_12)
    return RDS_2_adapt











def alg_estimador_manager(n_metodo_est,G,cant_nodos, dict_nei_G, list_nei_G, dic_grados_G, set_infectados,Muestra_general, tipo_muestreo=1,est_2E=1, cant_aprox_hop_need=1 ):
    est_N_B=-10    
    # definicion del diccionario con los algs est
    L_name_est=["PIMPLE_prom","EMV","gnsum","gnsum_ajustado","RDS I v1","RDS I v2","RDS I v2_ajustado","RDS II","RDS IIv2_ajustado"] 

    metodo_est=L_name_est[n_metodo_est - 1]
    if metodo_est=="PIMPLE_prom":
        est_N_u=N_u_PIMPLE_prom(Muestra_general, G, set_infectados,list_nei_G,dic_grados_G)  # 12/10 /2022

    elif metodo_est=="EMV":       
        est_N_u=N_u_EMV(G, set_infectados, Muestra_general,list_nei_G,dic_grados_G)  # 12 / 10 /2022

    elif metodo_est=="gnsum":
        input_base=[cant_nodos, dict_nei_G,dic_grados_G, set_infectados]
        input_muestra=[Muestra_general,tipo_muestreo,est_2E,cant_aprox_hop_need]
        input_comb="2"
        est_N_u=gnsum(input_base,input_muestra,input_comb)
        
    elif metodo_est=="gnsum_ajustado":
        if 2*len(set_infectados) >= cant_nodos:
            input_base=[cant_nodos, dict_nei_G,dic_grados_G, set_infectados]
            input_muestra=[Muestra_general,tipo_muestreo,est_2E,cant_aprox_hop_need]
            input_comb="2"
            est_N_u=gnsum(input_base,input_muestra,input_comb)
        else:
            set_nodes=set(G.nodes())
            set_noinfectados=set_nodes- set_infectados
            input_base=[cant_nodos, dict_nei_G,dic_grados_G, set_noinfectados]
            input_muestra=[Muestra_general,tipo_muestreo,est_2E,cant_aprox_hop_need]
            input_comb="2"
            est_N_B=gnsum(input_base,input_muestra,input_comb)            
            est_N_u=cant_nodos-est_N_B            

        
    elif metodo_est=="RDS I v1":
        input_base=[G,cant_nodos,dict_nei_G,dic_grados_G,set_infectados]
        input_muestra=[Muestra_general,tipo_muestreo ]
        est_N_u=estimator_RDS_I_V1(input_base,input_muestra )        
    elif metodo_est=="RDS I v2":
        input_base=[G,cant_nodos,dict_nei_G,dic_grados_G,set_infectados]
        input_muestra= [Muestra_general,tipo_muestreo ]
        est_N_u=estimador_rds1_v2(input_base, input_muestra)
        
        
    elif metodo_est=="RDS I v2_ajustado":
        if 2*len(set_infectados) >= cant_nodos:
            input_base=[G,cant_nodos,dict_nei_G,dic_grados_G,set_infectados]
            input_muestra= [Muestra_general,tipo_muestreo ]
            est_N_u=estimador_rds1_v2(input_base, input_muestra)
        else:
            
            set_nodes=set(G.nodes())
            set_noinfectados=set_nodes- set_infectados
            input_base=[G,cant_nodos,dict_nei_G,dic_grados_G,set_noinfectados]
            input_muestra= [Muestra_general,tipo_muestreo ]
            est_N_B=estimador_rds1_v2(input_base, input_muestra)        
            est_N_u=cant_nodos-est_N_B
        
    elif metodo_est=="RDS II":      
        input_base=[cant_nodos,set_infectados,dic_grados_G]
        input_muestra=[tipo_muestreo, Muestra_general]
        est_N_u=estimador_rds2(input_base,input_muestra)
        print("RDS II",est_N_u)
    elif metodo_est=="RDS IIv2_ajustado":     
        if 2*len(set_infectados) >= cant_nodos:        
            input_base=[cant_nodos,set_infectados,dic_grados_G]
            input_muestra=[tipo_muestreo, Muestra_general]
            est_N_u=estimador_rds2(input_base,input_muestra)        
        else:
            
            set_nodes=set(G.nodes())
            set_noinfectados=set_nodes- set_infectados
            input_base=[cant_nodos,set_noinfectados,dic_grados_G]
            input_muestra=[tipo_muestreo, Muestra_general]
            est_N_B=estimador_rds2(input_base, input_muestra)        
            est_N_u=cant_nodos-est_N_B            
        print("RDS II ajustado",est_N_u)
    if est_N_B == 0:
        return 0   
    
    return est_N_u,metodo_est
    
def N_u_PIMPLE_parcial(x_i,G, set_infectados,list_nei_G,dic_grados_G):  # 12/10 /2022
    # i=nodo del muestreo (any)
    # G: grafo G
    # lista_infectados: nodos infectados del grafo G
    N=len(G.nodes)
    neig_i=list_nei_G[x_i]
    yiu= len(neig_i & set_infectados)
    #Calcular EMV d_i:
    d_i= dic_grados_G[x_i]
    
    #Calcular el esti N_i^u:
    N_iu=N * yiu / d_i
    return N_iu

def N_u_PIMPLE_prom(Muestra_general, G, set_infectados,list_nei_G,dic_grados_G):  # 12/10 /2022
    list_N_iu=[]
    for x_i in Muestra_general:
        N_iu=N_u_PIMPLE_parcial(x_i,G, set_infectados,list_nei_G,dic_grados_G)
        list_N_iu.append(N_iu)

    N_u_PIMPLE=np.mean(list_N_iu)   
    
    return N_u_PIMPLE



def N_u_EMV(G, set_infectados, Muestra_general,list_nei_G,dic_grados_G):  # 12 / 10 /2022
    list_d_i=[]
    list_yiu=[]
    N=len(G.nodes)

    for x_i in Muestra_general:
        neig_i=list_nei_G[x_i]
        yiu= len(neig_i & set_infectados) 
        list_yiu.append(yiu)    

        #Calcular EMV d_i:
        d_i= dic_grados_G[x_i]    
        list_d_i.append(d_i)
    #Calcular el esti N^u:
    N_u=N * sum(list_yiu)/sum(list_d_i)

    return N_u


def Hansen_Hurwitz_est_num(input_base, input_muestra):
    cant_nodos, dict_nei_G,dic_grados_G, set_infectados=input_base[0],input_base[1],input_base[2],input_base[3]
    Muestra_general,tipo_muestreo =input_muestra[0],input_muestra[1]
    
    tipo_sample=tipo_muestreo
    # SOn todos muestreos con repeticion, plt len(Muestra_general) = al tamaño del sample
    if tipo_sample=="5":
        
        #prob_sel_i=Prob_sel_sample_tip_5(cant_nodos)  # cte
        suma=0
        for x_i in Muestra_general:
            y_ih=len(dict_nei_G[x_i] & set_infectados)
            suma=suma+y_ih
        est_yih=suma*cant_nodos/len(Muestra_general)       
    elif tipo_sample=="6" or tipo_sample=="7" or tipo_sample=="8":
        Muestra_general=Muestra_general[1:]
        suma=0
        for x_i in Muestra_general:
            suma=suma + 1/dic_grados_G[x_i]
        est_2aristas=cant_nodos*len(Muestra_general)/suma

        #cte,prob_sel_i=Prob_sel_sample_tip_678(cant_nodos , dic_grados_G , est_2aristas, Muestra_general)  # cte, dict
        
        for x_i in Muestra_general:
            y_ih=len(dict_nei_G[x_i] & set_infectados)
            suma=suma+y_ih/dic_grados_G[x_i]
        est_yih=suma*est_2aristas/len(Muestra_general)       
    else:
        print("Tipo de muestreo no corresponde")
        est_yih=-1
        return est_yih
    
    
    return est_yih

def Hansen_Hurwitz_est_den(input_base, input_muestra):
    cant_nodos, dict_nei_G,dic_grados_G, set_infectados=input_base[0],input_base[1],input_base[2],input_base[3]
    Muestra_general,tipo_muestreo =input_muestra[0],input_muestra[1]

    tipo_sample=tipo_muestreo
    # SOn todos muestreos con repeticion, plt len(Muestra_general) = al tamaño del sample
    if tipo_sample=="5":
        suma_1=0
        suma_2=0    
        for x_i in Muestra_general:
            func_directriz=len({x_i} & set_infectados)
            suma_1=suma_1+func_directriz*dic_grados_G[x_i]
            suma_2=suma_2 + func_directriz
    elif tipo_sample=="6" or tipo_sample=="7" or tipo_sample=="8":
        Muestra_general=Muestra_general[1:]
        suma_1=0
        suma_2=0    
        for x_i in Muestra_general:
            func_directriz=len({x_i} & set_infectados)
            suma_1=suma_1+func_directriz
            suma_2=suma_2 + func_directriz/dic_grados_G[x_i]
    else:
        print("Tipo de muestreo no corresponde")
        delta_A=-1
        return delta_A
    if suma_2 !=0:
        delta_A=suma_1/ suma_2   
    else:
        delta_A=-1
    return delta_A

def gnsum(input_base,input_muestra,input_comb):   # muestreo tipo 5,6,7,8
    num=input_comb
    [cant_nodos, dict_nei_G,dic_grados_G, set_infectados]=input_base[:]
    [Muestra_general,tipo_muestreo]= input_muestra[0], input_muestra[1]
    if tipo_muestreo =="2" or tipo_muestreo =="6":
        est_2E,cant_aprox_hop_need=input_muestra[2],input_muestra[3]
    
    est_yih=Hansen_Hurwitz_est_num(input_base, input_muestra) # muestreo tipo 5,6,7,8
    
    est_V_HF_mean=Hansen_Hurwitz_est_den(input_base, input_muestra) # muestreo tipo 5,6,7,8
    if est_V_HF_mean!=-1:
        est_Nu=est_yih/est_V_HF_mean
    else:
        est_Nu=0
    return est_Nu


#  EStimador RDS I   

def estimacion_Cab(input_base,Muestra_general):
    [dict_nei_G,dic_grados_G,set_infectados]=input_base    
    # Estimacion de C_ab, C_ba
    suma_aa,suma_ab,suma_ba,suma_bb=0,0,0,0
    for x_i in Muestra_general:
        func_directriz=len( {x_i} & set_infectados )
        if func_directriz==1: # se obtiene raa y rab
            nei_xi=dict_nei_G[x_i]
            r_aa= len(nei_xi & set_infectados)  
            r_ab= dic_grados_G[x_i]-r_aa
            suma_aa= suma_aa + r_aa
            suma_ab= suma_ab + r_ab        
        else:
            nei_xi=dict_nei_G[x_i]
            r_ba= len(nei_xi & set_infectados)  
            r_bb= dic_grados_G[x_i]-r_ba
            suma_ba= suma_ba + r_ba
            suma_bb= suma_bb + r_bb        
    if suma_aa + suma_ab ==0 or suma_bb + suma_ba==0:
        return -1,-1
    est_C_ab=suma_ab/(suma_aa + suma_ab)
    est_C_ba=suma_ba/(suma_bb + suma_ba)
    
    return est_C_ab,est_C_ba

def estimator_RDS_I_V1(input_base,input_muestra):
    [G,cant_nodos,dict_nei_G,dic_grados_G,set_infectados]=input_base
    [Muestra_general,tipo_muestreo ]=input_muestra

    input_base=[dict_nei_G,dic_grados_G,set_infectados]
    est_C_ab,est_C_ba=estimacion_Cab(input_base,Muestra_general)

    input_base=[cant_nodos, dict_nei_G,dic_grados_G, set_infectados]
    input_muestra= [Muestra_general,tipo_muestreo ]
    est_D_A=Hansen_Hurwitz_est_den(input_base, input_muestra)

    set_nodos=set(G.nodes())
    set_noinf=set_nodos - set_infectados
    input_base=[cant_nodos, dict_nei_G,dic_grados_G, set_noinf]
    input_muestra= [Muestra_general,tipo_muestreo ]
    est_D_B=Hansen_Hurwitz_est_den(input_base, input_muestra)
    
    if est_C_ab ==-1 or est_C_ba==-1:
        return 0
    if est_D_B!=-1 and est_D_A!=-1:
        N_A=cant_nodos* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)
        N_B=cant_nodos* (est_D_A*est_C_ab)/(est_D_A*est_C_ab + est_D_B*est_C_ba)
    else:
        N_A=0
    return N_A


##  -----------------------   Estimador RDS I version 2

def estimador_2arista(input_base,input_muestra):    # input con nodo semilla, estimacion sin nodo semilla
    [cant_nodos,dic_grados_G]=input_base
    [tipo_muestreo,Muestra_general]=input_muestra
    # estimar estimar est_D_A # estimar total_aristas  
    # total_aristas=len(G.edges)
    # caso muestreo tipo 5:
    if tipo_muestreo== "5":
        suma=0
        for x_i in Muestra_general:
            d_i=dic_grados_G[x_i]
            suma=suma+d_i
        avg_degree=suma/len(Muestra_general)
        est_2aristas=avg_degree*cant_nodos

    if tipo_muestreo== "6" or tipo_muestreo== "7" or tipo_muestreo== "8":    
        Muestra_general=Muestra_general[1:]
        suma=0
        for x_i in Muestra_general:
            d_i=dic_grados_G[x_i]
            suma=suma+1/d_i
        avg_degree=len(Muestra_general)/suma
        est_2aristas=avg_degree*cant_nodos

    return est_2aristas

        

def estimador_rds1_v2(input_base, input_muestra):
    [G,cant_nodos,dict_nei_G,dic_grados_G,set_infectados]=input_base
    [Muestra_general,tipo_muestreo ]=input_muestra
    
    input_base=[dict_nei_G,dic_grados_G,set_infectados]
    est_C_ab,est_C_ba=estimacion_Cab(input_base,Muestra_general)

    set_nodos=set(G.nodes())
    set_noinf=set_nodos - set_infectados
    input_base=[cant_nodos, dict_nei_G,dic_grados_G, set_noinf]
    input_muestra= [Muestra_general,tipo_muestreo ]
    est_D_B=Hansen_Hurwitz_est_den(input_base, input_muestra)

    if est_C_ab ==-1 or est_C_ba==-1:
        return 0
    
    if est_D_B==-1:
        return 0
    
    input_base=[cant_nodos,dic_grados_G]
    input_muestra=[tipo_muestreo,Muestra_general]
    est_2aristas=estimador_2arista(input_base,input_muestra) # input con nodo semilla, estimacion sin nodo semilla
    beta=est_C_ba*est_D_B
    est_D_A= (est_2aristas * beta) /( cant_nodos*(beta +est_D_B *est_C_ab )- est_2aristas*est_C_ab )

    N_A=cant_nodos* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)
    # N_B=cant_nodos* (est_D_A*est_C_ab)/(est_D_A*est_C_ab + est_D_B*est_C_ba)
    return N_A

##  -----------------------   Estimador RDS II version 
def estimador_rds2(input_base,input_muestra):
    cant_nodos,set_infectados,dic_grados_G=input_base[0],input_base[1],input_base[2]
    tipo_muestreo,Muestra_general=input_muestra[0],input_muestra[1]
    
    suma_1=0
    suma_2=0
    if tipo_muestreo=="5":
        for x_i in Muestra_general:
            func_directriz=len( {x_i} & set_infectados )
            suma_1=func_directriz+suma_1
        est_Nu= cant_nodos* suma_1/len(Muestra_general)
    else:        
        Muestra_general=Muestra_general[1:]
        for x_i in Muestra_general:
            func_directriz=len( {x_i} & set_infectados )
            delta_i=dic_grados_G[x_i]
            suma_1=func_directriz/delta_i+suma_1
            suma_2=1/delta_i+suma_2

        est_Nu=cant_nodos*suma_1/suma_2
    return est_Nu

def gnsum_completo(input_base, input_muestra):
    cant_nodos, dict_nei_G,dic_grados_G, set_infectados=input_base[0],input_base[1],input_base[2],input_base[3]
    Muestra_general,tipo_muestreo =input_muestra[0],input_muestra[1]
    
    suma_1,suma_2,suma_3,suma_4=0,0,0,0
    M=len(Muestra_general)
    x_0=Muestra_general[0]
    y_0h=cant_nodos*len(dict_nei_G[x_0] & set_infectados)
    
    f_directriz=len( {x_0} & set_infectados)
    num_h0=f_directriz*dic_grados_G[x_0]*cant_nodos
    den_h0=f_directriz*cant_nodos
    for x_i in Muestra_general[1:]:
        suma_1=suma_1+len(dict_nei_G[x_i] & set_infectados)/dic_grados_G[x_i]
        suma_2=suma_2 + 1/dic_grados_G[x_i]   # estimacion 2 |E|
        
        f_directriz=len( {x_i} & set_infectados)
        suma_3=suma_3+ f_directriz
        suma_4=suma_4+ f_directriz/dic_grados_G[x_i]
        
        
    est_2aristas=cant_nodos*(M-1)/suma_2
    est_yuh= (y_0h+est_2aristas*suma_1)/M
    
    v_hu_mean=(num_h0+est_2aristas*suma_3)/(den_h0+est_2aristas*suma_4)
    est_N_u=est_yuh/v_hu_mean

    return est_N_u

def estimator_RDS_I_V1_completo(input_base,input_muestra):
    [cant_nodos,dict_nei_G,dic_grados_G,set_infectados, set_noinfectados]=input_base    
    [Muestra_general ]=input_muestra

    input_base=[dict_nei_G,dic_grados_G,set_infectados]
    est_C_ab,est_C_ba=estimacion_Cab(input_base,Muestra_general)
    
    suma_2,suma_3,suma_4,suma_5,suma_6=0,0,0,0,0  #
    M=len(Muestra_general)
    
    x_0=Muestra_general[0]
    
    f_directriz=len( {x_0} & set_infectados)
    num_h0_infec=f_directriz*dic_grados_G[x_0]*cant_nodos
    den_h0_infec=f_directriz*cant_nodos
    
    
    g_directriz=len( {x_0} & set_noinfectados)
    num_h0_noinf=g_directriz*dic_grados_G[x_0]*cant_nodos
    den_h0_noinf=g_directriz*cant_nodos    
    for x_i in Muestra_general[1:]:
        suma_2=suma_2 + 1/dic_grados_G[x_i]   # estimacion 2 |E|
        
        f_directriz=len( {x_i} & set_infectados)
        suma_3=suma_3+ f_directriz
        suma_4=suma_4+ f_directriz/dic_grados_G[x_i]  
        
        g_directriz=len( {x_i} & set_noinfectados)    
        suma_5=suma_5+ g_directriz
        suma_6=suma_6+ g_directriz/dic_grados_G[x_i]  
    
    est_2aristas=cant_nodos*(M-1)/suma_2
    est_D_A=(num_h0_infec+est_2aristas*suma_3)/(den_h0_infec+est_2aristas*suma_4)
    est_D_B=(num_h0_noinf+est_2aristas*suma_5)/(den_h0_noinf+est_2aristas*suma_6)

    if est_C_ab ==-1 or est_C_ba==-1:
        return 0
    if est_D_B!=-1 and est_D_A!=-1:
        N_A=cant_nodos* (est_D_B*est_C_ba)/(est_D_A*est_C_ab + est_D_B*est_C_ba)
        
    else:
        N_A=0
    return N_A


