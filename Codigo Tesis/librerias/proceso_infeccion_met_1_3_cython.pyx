# Archivo: proceso_infeccion_met_1_3.pyx

import numpy as np
cimport numpy as np
from scipy.stats import bernoulli

# Agrega estas líneas para acceder a los encabezados de NumPy
np.import_array()

def proceso_infeccion_met_1_3(list last_list_aceptados, list list_nei_G, dict dic_grados_G, double prob_acept):
    # resetear lista de pendientes:
    cdef list last_list_aceptados_new = []
    cdef int d_xi, ii
    cdef list nei_xi
    cdef np.ndarray eva_list_pend
    
    for x_i in last_list_aceptados:
        # calcular vecinos de ultimos nodos aceptados
        nei_xi = list_nei_G[x_i]
        d_xi = dic_grados_G[x_i]
        eva_list_pend = bernoulli.rvs(prob_acept, size=d_xi)
        
        for ii in range(d_xi):
            if eva_list_pend[ii] == 1:
                last_list_aceptados_new.append(nei_xi[ii])
    
    return last_list_aceptados_new
