import networkx as nx
import random
from librerias  import alg_fabrica_grafo as alg_grafo
import pickle 

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

    if m < 1 or  m >=n:
        raise nx.NetworkXError("Barabási–Albert network must have m >= 1"
                               " and m < n, m = %d, n = %d" % (m, n))
    if seed is not None:
        random.seed(seed)
 
    # Add m initial nodes (m0 in barabasi-speak)
    G=nx.empty_graph(m)

    # Target nodes for new edges
    targets=list(range(m))
    # List of existing nodes, with nodes repeated once for each adjacent edge
    repeated_nodes=[]
    # Start adding the other n-m nodes. The first node is m.
    source=m   # nodo que se agrega en la iteracion
    while source<n:
        # Add edges to m nodes from the source.
        G.add_edges_from(zip([source]*m,targets)) # [m,...,m][0,...,source=m-1] iteracion 1
        
        
        # guardar en una lista los nodos que se le agregaron aristas hasta el nodo source
        repeated_nodes.extend(targets)
        # And the new node "source" has m edges to add to the list.
        repeated_nodes.extend([source]*m)
        
        # Elegir nodos que se enlazarán al nuevo nodo source en la siguiente iteración
        # Pick uniformly from repeated_nodes (preferential attachment)
        targets = alg_grafo._random_subset(repeated_nodes, m)
        source += 1
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
        G=alg_grafo.barabasi_albert_graph(n, m, seed=None)

    elif Fuente_grafo=="algoritmo_aleatorio":
        n=config[0]
        p=config[1]        
        G=nx.erdos_renyi_graph(n, p, seed=None, directed=False)    
    return G
    