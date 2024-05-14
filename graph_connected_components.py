

import pandas as pd
import networkx as nx
from sklearn.neighbors import radius_neighbors_graph
import numpy as np
import plotly.graph_objects as go

def rgb_to_hex(rgb):
    '''
    expected input is like: (190, 255, 255) etc
    expected output is a hex color: '#beffff' or whatever
    '''
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

datapath = 'if_data_short/0197_IF1.csv'
data = pd.read_csv(datapath)
print('Data loaded')

def create_phenotype_graph(data, radius, min_connection, chosen_cell_types):
    chosen_data = data[data['cell type'].str.contains(chosen_cell_types)]
    coordinates = chosen_data[['nucleus.x', 'nucleus.y']].values
    graph = radius_neighbors_graph(coordinates, radius=radius, mode='distance')
    converted = nx.convert_matrix.from_scipy_sparse_array(graph)
    connected = list(nx.connected_components(converted))
    interesting_connections = [connection for connection in connected if len(connection) > min_connection]
    chosen_data = chosen_data.reset_index(drop=True)
    chosen_data["cluster_index"] = None
    for cluster_i, connection in enumerate(interesting_connections):
        for number in connection:
            chosen_data.loc[number, 'cluster_index'] = cluster_i
    print(chosen_data)


    fig = go.Figure()

    grouped = chosen_data.groupby('cluster_index')
    for name, group in sorted(grouped):
        fig.add_trace(
            go.Scattergl(
                mode='markers',
                x=group['nucleus.x'],
                y=group['nucleus.y'],
                marker=dict(
                    color=group['color'],
                    size=10  # Increased marker size
                ),
                name=name,
                hovertemplate=(
                        "Cell Type: %{customdata[0]}<br>" +
                        "Cell ID: %{customdata[1]}<br>" +
                        "Cluster index: %{customdata[2]}<br>" +
                        "Phenotype:%{customdata[3]}<br>" +
                        "X:%{customdata[4]}<br>" +
                        "Y:%{customdata[5]}"
                ),
                customdata=np.stack((group['cell type'], group['cell.ID'],
                                     group['cluster_index'],group['phenotype'],
                                    group['nucleus.x'], group['nucleus.y']), axis=-1),
                showlegend=True
            )
        )

        # Add lines between connected cells


    fig.update_layout(
        title=f'Interactive scatter plot of cell data neighbor radius {radius}, showing at least {min_connection}-cell connected components',
        width=1600,
        height=900,
        template='plotly_dark'
    )

    fig.show()
    fig.write_html("plot_clustered.html")
    print('plot loaded')

print(create_phenotype_graph(data, 30, 20, "Tcell|T cell|Treg|Bcell|DC|NK"))
#picked cell types most likely for TLS
