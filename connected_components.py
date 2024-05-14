

import pandas as pd
import networkx as nx
from sklearn.neighbors import radius_neighbors_graph
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import numpy as np



def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb

IF1_cell_mapping = {"other": rgb_to_hex((190, 190, 190)),
                    "CD15+Tumor": rgb_to_hex((73, 176, 248)),
                    "CD15-Tumor": rgb_to_hex((138, 79, 45)),
                    "Tcell": rgb_to_hex((235, 74, 148)),
                    "Bcell": rgb_to_hex((204, 49, 31)),
                    "BnTcell": rgb_to_hex((236, 95, 42)),
                    "Neutrophil": rgb_to_hex((0, 40, 245)),
                    "Macrophage": rgb_to_hex((97, 209, 62)),
                    "DC": rgb_to_hex((49, 113, 30))}

datapath = 'if_data_short/0316_IF1.csv'
data = pd.read_csv(datapath)
print('Data loaded')

def create_proportion_plot(data, radius, min_connection, chosen_cell_types):
    chosen_data = data[data['cell type'].str.contains(chosen_cell_types)]
    coordinates = chosen_data[['nucleus.x', 'nucleus.y']].values
    graph = radius_neighbors_graph(coordinates, radius=radius, mode='distance')
    converted = nx.convert_matrix.from_scipy_sparse_array(graph)
    connected = list(nx.connected_components(converted))
    interesting_connections = [connection for connection in connected if len(connection) > min_connection]
    chosen_data = chosen_data.reset_index(drop=True)
    # Calculate percent proportions of cell types in each connected component
    proportions = []
    for connection in interesting_connections:
        sub_data = chosen_data.iloc[list(connection)]
        counts = sub_data['cell type'].value_counts()
        total = counts.sum()
        percent = (counts / total * 100).to_dict()
        proportions.append(percent)
    print(proportions)
    # proportions gives the percent proportions for each connected component in interesting_components
    # Plotting
    df = pd.DataFrame(proportions)
    df.plot(kind="bar", stacked=True, color=IF1_cell_mapping, width=0.7)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.xlabel(f'Connected Components')
    plt.ylabel('Percentage')
    plt.title(f'Proportions of Cell Types in Each Component: Each component with minimal length {min_connection}, neighbor radius {radius}')
    plt.show()



print(create_proportion_plot(data, 30, 50, "Tcell|T cell|Treg|Bcell|DC|NK"))
#I picked cell types most likely for TLS