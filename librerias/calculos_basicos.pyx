import random

def calculos_basicos(tipo_muestreo,Muestra_general,input_base_2):
    [N,set_nodos,dict_nei_G,dic_grados_G,porc_nodos_inf,set_infectados]=input_base_2
    M=len(Muestra_general)
    set_noinfectados=set_nodos-set_infectados

    if tipo_muestreo == "5":
        suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8=[0]*8
        num_c_ab,den_c_ab=0,0
        num_c_ba,den_c_ba=0,0
        
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
            
            
            l_u=random.sample(neig_i, 1)
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
        input_base=[porc_nodos_inf,tipo_muestreo,N,M,suma_1,suma_2,suma_3,suma_4,suma_5,suma_6,suma_7,suma_8]
        print("calculos basicos sin cython")
        return input_base,est_cab,est_cba
    



def calculos_basicos_cython(str tipo_muestreo, list Muestra_general, list input_base_2):
    cdef int N, M
    cdef set set_nodos, set_infectados, set_noinfectados
    cdef dict dict_nei_G, dic_grados_G
    cdef float suma_1, suma_2, suma_3, suma_4, suma_5, suma_6, suma_7, suma_8
    cdef int num_c_ab, den_c_ab, num_c_ba, den_c_ba
    cdef int x_i, d_i, yih, yiJ, func_directriz, g_directriz
    cdef list l_u
    cdef int u, u_func_directriz, u_g_directriz
    cdef float est_cab, est_cba
    cdef list input_base

    [N, set_nodos, dict_nei_G, dic_grados_G, porc_nodos_inf, set_infectados] = input_base_2
    M = len(Muestra_general)
    set_noinfectados = set_nodos - set_infectados

    if tipo_muestreo == "5":
        suma_1 = suma_2 = suma_3 = suma_4 = suma_5 = suma_6 = suma_7 = suma_8 = 0
        num_c_ab = den_c_ab = num_c_ba = den_c_ba = 0

        for x_i in Muestra_general:
            neig_i = dict_nei_G[x_i]
            d_i = dic_grados_G[x_i]
            yih = len(neig_i & set_infectados)
            yiJ = d_i - yih

            func_directriz = len({x_i} & set_infectados)
            g_directriz = -func_directriz + 1

            suma_1 += yih / d_i
            suma_2 += yih
            suma_3 += d_i
            suma_4 += func_directriz
            suma_5 += func_directriz * d_i
            suma_6 += g_directriz
            suma_7 += g_directriz * d_i
            suma_8 += yiJ

            l_u = random.sample(neig_i, 1)
            u = l_u[0]
            u_func_directriz = len({u} & set_infectados)
            u_g_directriz = -u_func_directriz + 1
            num_c_ab += d_i * func_directriz * u_g_directriz
            den_c_ab += d_i * func_directriz
            num_c_ba += d_i * g_directriz * u_func_directriz
            den_c_ba += d_i * g_directriz

        if den_c_ab == 0 or den_c_ba == 0:
            est_cab = est_cba = 0
        else:
            est_cab = num_c_ab / den_c_ab
            est_cba = num_c_ba / den_c_ba

        #input_base = [porc_nodos_inf, tipo_muestreo, N, M, suma_1, suma_2, suma_3, suma_4, suma_5, suma_6, suma_7, suma_8]
        print("calculos basicos cython v_1")
        #return input_base, est_cab, est_cba
        
        return porc_nodos_inf, tipo_muestreo, N, M, suma_1, suma_2, suma_3, suma_4, suma_5, suma_6, suma_7, suma_8,est_cab, est_cba


def calculos_basicos_cython_v2(str tipo_muestreo, list Muestra_general, int N,set set_nodos,dict dict_nei_G,dict dic_grados_G,float porc_nodos_inf,set set_infectados):
    
    cdef set set_noinfectados
    cdef float suma_1, suma_2, suma_3, suma_4, suma_5, suma_6, suma_7, suma_8
    cdef int num_c_ab, den_c_ab, num_c_ba, den_c_ba
    cdef int x_i, d_i, yih, yiJ, func_directriz, g_directriz
    cdef list l_u
    cdef int u, u_func_directriz, u_g_directriz
    cdef float est_cab, est_cba
    

    M = len(Muestra_general)
    set_noinfectados = set_nodos - set_infectados

    if tipo_muestreo == "5":
        suma_1 = suma_2 = suma_3 = suma_4 = suma_5 = suma_6 = suma_7 = suma_8 = 0
        num_c_ab = den_c_ab = num_c_ba = den_c_ba = 0

        for x_i in Muestra_general:
            neig_i = dict_nei_G[x_i]
            d_i = dic_grados_G[x_i]
            yih = len(neig_i & set_infectados)
            yiJ = d_i - yih

            func_directriz = len({x_i} & set_infectados)
            g_directriz = -func_directriz + 1

            suma_1 += yih / d_i
            suma_2 += yih
            suma_3 += d_i
            suma_4 += func_directriz
            suma_5 += func_directriz * d_i
            suma_6 += g_directriz
            suma_7 += g_directriz * d_i
            suma_8 += yiJ

            l_u = random.sample(neig_i, 1)
            u = l_u[0]
            u_func_directriz = len({u} & set_infectados)
            u_g_directriz = -u_func_directriz + 1
            num_c_ab += d_i * func_directriz * u_g_directriz
            den_c_ab += d_i * func_directriz
            num_c_ba += d_i * g_directriz * u_func_directriz
            den_c_ba += d_i * g_directriz

        if den_c_ab == 0 or den_c_ba == 0:
            est_cab = est_cba = 0
        else:
            est_cab = num_c_ab / den_c_ab
            est_cba = num_c_ba / den_c_ba
        print("calculos basicos cython v_2")
        #input_base = [porc_nodos_inf, tipo_muestreo, N, M, suma_1, suma_2, suma_3, suma_4, suma_5, suma_6, suma_7, suma_8]
        #return input_base, est_cab, est_cba
        return porc_nodos_inf, tipo_muestreo, N, M, suma_1, suma_2, suma_3, suma_4, suma_5, suma_6, suma_7, suma_8,est_cab, est_cba
    