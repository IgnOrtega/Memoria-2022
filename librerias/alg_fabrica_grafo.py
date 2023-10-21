import networkx as nx
import random
import pickle 
from datetime import datetime

def made_lot_graph(densidad,tipo_grafo,n,parametro,inicio,termino):
        save_direccion=nombre_ubi_save(densidad, tipo_grafo,parametro)
        if tipo_grafo=="escala":
            tipo_grafo="algoritmo_escala"
            alpha_1="m"
            alpha=alpha_1+str(int(parametro))        
        else:
            tipo_grafo="algoritmo_aleatorio"
            alpha_1="p"
            alpha=alpha_1+str(int(1/parametro))    

        config=[n]
        config.append(parametro)

        
        for num_graph in range(inicio,termino):
                    date_1=datetime.now()                         
                    #print(tipo_grafo, config)
                    G= made_graph(tipo_grafo, config)

                    nuevo_nombre ="densidad_"+densidad+"graph_n"+str(n)+alpha+"n°"+str(num_graph)

                    var_save= G # variable a guardar
                    # creamos fichero binario 
                    fichero_bin= open( save_direccion+nuevo_nombre, "wb")
                    # metemos la var_wsave en el fichero binario
                    pickle.dump(var_save, fichero_bin)
                    fichero_bin.close()    
                    #print(save_direccion+nuevo_nombre)
                    date_2=datetime.now()     
                    print(num_graph, date_2,date_1,date_2-date_1)                    

def nombre_ubi_save(densidad, tipo_grafo,parametro):
    save_direccion="Datos_grafos/"
    if densidad=="pesada":
        carpeta_1="densidad_pesada/"
    elif densidad=="media":
        carpeta_1="densidad_media/"
        
    elif densidad=="ligera":    
        carpeta_1="densidad_ligera/"
    else:    
        carpeta_1="densidad_super_ligera/"        

    carpeta_2=tipo_grafo+"/"
    if tipo_grafo== "escala":
        carpeta_3="m "+str(int(parametro))+"/"
    else:
        carpeta_3="p "+str(int(1/parametro))+"/"

    save_direccion=save_direccion+carpeta_1+carpeta_2+carpeta_3
    return save_direccion




def _random_subset(seq, m):
    """Return m unique elements from seq.

    This differs from random.sample which can return repeated
    elements if seq holds repeated elements.

    Note: rng is a random.Random or numpy.random.RandomState instance.
    """
    targets = set()
    while len(targets) < m:
        #x = rng.choice(seq)
        x=random.choice(seq)
        targets.add(x)
    return targets

def barabasi_albert_graph(n, m, seed=None):

# Por cada iteracion agregra m aristas
# Si source=m al inicio y llega hasta n, entonces son (n-m) iteracion
# La cantidad de aristas es: (cant de iteraciones) x (cant de aristas agregadas)
# Conclusion: G tiene (n-m)m aristas 


    if m < 1 or  m >=n:
        raise nx.NetworkXError("Barabási–Albert network must have m >= 1"
                               " and m < n, m = %d, n = %d" % (m, n))
    if seed is not None:
        random.seed(seed)
 
    # Add m initial nodes (m0 in barabasi-speak)
    G=nx.empty_graph(m)

    # Target nodes for new edges
    targets=list(range(m)) # [0,1,2,...,m-1]
    # List of existing nodes, with nodes repeated once for each adjacent edge
    repeated_nodes=[]
    # Start adding the other n-m nodes. The first node is m.
    source=m   # el nodo m el nodo que se le agregaran aristas en el iteracion
    while source<n:
        # Add edges to m nodes from the source.
        G.add_edges_from(zip([source]*m,targets)) # [m,...,m][0,...,source=m-1] iteracion 1
        
        ## Guardar en una lista los nodos que se le agregaron aristas en esta iteracion
        # guardar en una lista los nodos que se le agregaron aristas hasta el nodo source (repeated_nodes=repeated_nodes+targets, esta sumando 2 listas)
        repeated_nodes.extend(targets)
        # And the new node "source" has m edges to add to the list.
        repeated_nodes.extend([source]*m)
        
        # Elegir m nodos que se enlazarán al nuevo nodo source en la siguiente iteración ;
        # los nodos se eligen con prob proporcional a la cant de aristas que tienen hasta el momento;
        # Pick uniformly from repeated_nodes (preferential attachment)
        targets = _random_subset(repeated_nodes, m)   # es un subset [0,1,2,...,source(al inicio iter )]
        source += 1 # cambiar el nodo fuente
    return G



def made_graph(Fuente_grafo, config):
    # Revisar de donde cargar el grafo
    #Opcion: grafo de bbdd
    if Fuente_grafo== "bbdd":    
        if config[1]== "internet":
            name=config[0]
            # cargar fichero binario
            fh = open(name, 'rb' )
            G = nx.read_edgelist(fh)     
            G = nx.convert_node_labels_to_integers(G)  
            G = nx.to_undirected(G)
        else:
            name=config[0]
            file = open(name, 'rb')
            # dump information to that filea
            G = pickle.load(file)
            #close the file
            file.close()    

    #Opcion: Creacion grafo de algoritmos
    elif Fuente_grafo=="algoritmo_escala":
        n=config[0]
        m=config[1]    
        G= barabasi_albert_graph(n, m, seed=None)

    elif Fuente_grafo=="algoritmo_aleatorio":
        n=config[0]
        p=config[1]        
        G=nx.erdos_renyi_graph(n, p, seed=None, directed=False)    
    return G

    