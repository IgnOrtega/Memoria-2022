def compute_c_ab_sample_7_8_cython(list sample_aristas, set set_infectados):
    cdef int v_f_directriz, v_g_directriz
    cdef int u_f_directriz, u_g_directriz
    cdef double est_cab,est_cba

    cdef int suma_num_c_ab = 0, suma_den_c_ab = 0
    cdef int suma_num_c_ba = 0, suma_den_c_ba = 0

    for v, u in sample_aristas:
        v_f_directriz = int(v in set_infectados)
        v_g_directriz = 1 - v_f_directriz

        u_f_directriz = int(u in set_infectados)
        u_g_directriz = 1 - u_f_directriz

        suma_num_c_ab += v_f_directriz * u_g_directriz
        suma_den_c_ab += v_f_directriz
        suma_num_c_ba += v_g_directriz * u_f_directriz
        suma_den_c_ba += v_g_directriz

    if suma_den_c_ab == 0 or suma_den_c_ba == 0:
        est_cab = 0
        est_cba = 0
    else:
        est_cab = suma_num_c_ab / suma_den_c_ab
        est_cba = suma_num_c_ba / suma_den_c_ba

    return est_cab, est_cba
