import random
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def proceso_infeccion_met_1_3_random(list last_list_aceptados, list list_nei_G, double prob_acept):
    cdef list last_list_aceptados_new = []
    cdef int x_i, u
    cdef list nei_xi

    for x_i in last_list_aceptados:
        nei_xi = list_nei_G[x_i]
        for u in nei_xi:
            if random.random() < prob_acept:
                last_list_aceptados_new.append(u)
    
    return last_list_aceptados_new
