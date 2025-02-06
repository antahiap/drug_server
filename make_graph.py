
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from neo4j import GraphDatabase
from fa2 import ForceAtlas2


class MakeGraph:

    def __init__(self):
        return None
    
    def _scale_to_range(self, values, min_range=10, max_range=100):
        """
        Scales a list of values to a defined range [min_range, max_range].

        :param values: List of values to scale.
        :param min_range: Minimum value of the new scale (default is 10).
        :param max_range: Maximum value of the new scale (default is 100).
        :return: A list of scaled values.
        """
        if not values:
            return []

        min_value = min(values)
        max_value = max(values)

        if min_value == max_value:
            return [min_range] * len(values)
        scaled_values = [
            min_range + (value - min_value) / (max_value - min_value) * (max_range - min_range)
            for value in values
        ]

        return scaled_values

    def _node_heatmap(self, color_list, color_select, measure):

        selected_sizes = [measure[i] for i, c in enumerate(color_list) if c == color_select]
        print(selected_sizes, color_list, color_select)
        norm = plt.Normalize(min(selected_sizes), max(selected_sizes))
        cmap = cm.jet #Reds
        heat_colors = [cmap(norm(size)) for size in selected_sizes]

        updated_color = []
        j=0
        for i, ci in enumerate(color_list):
            if ci == color_select:
                updated_color.append(heat_colors[j])
                j+=1
            else:
                updated_color.append(ci)
        return updated_color

    def _calculate_pagerank_radius(self, G, alpha=0.85):
        """Calculates PageRank for nodes and assigns it as radius with linear scaling."""
        # Calculate the PageRank scores for each node
        pagerank_scores = nx.pagerank(G, alpha=alpha)

        return self._scale_to_range(pagerank_scores.values(), 0.1, 100)

    def _graph_to_json(self, G, positions):    

        """Visualizes the graph with ForceAtlas2 layout and SimRank-based node sizes"""
        node_sizes = self._calculate_pagerank_radius(G)
        edge_widths = self._scale_to_range(
                [G[u][v]['weight'] for u, v in G.edges()],  # Edge thickness based on weight
                1, 1000)

        color_dict = {
            'source': 'blue',
            'pathway': '#bab0ab',
            'drug': '#77b7b2', #'#769AB3',
            'cellular_component': 'green',
            'molecular_function': '#9c755f',
            'effect/phenotype': '#58a14e',
            'anatomy': '#4f79a7',
            'exposure': '#af7aa1',
            'biological_process': '#edc949',
            'gene/protein': '#ff9ea7',
            'disease': '#f28e2c'#'#FF8C00',
            }

        node_labels0 = nx.get_node_attributes(G, 'type')  
        node_colors = [color_dict.get(node_labels0.get(node, ''), 'gray') for node in G.nodes()]
        node_colors = self._node_heatmap(node_colors, 'blue', node_sizes)
        node_labels ={ node: (G.nodes[node].get('name') if G.nodes[node].get('type') == 'source' else '') for node in G.nodes}
        #node_labels ={ node: G.nodes[node].get('name') for node in G.nodes}


        # Json graph
        min_pos = np.min(list(positions.values()), axis=0)
        max_pos = np.max(list(positions.values()), axis=0)

        scale_factor = max(max_pos - min_pos)
        scaled_positions = {node: (pos[0] * 800 / scale_factor, pos[1] * 800 / scale_factor) for node, pos in positions.items()}    

        graph_data = {
        "nodes": [
            {
                "id": str(node),  # Ensure the ID is a string
                "label": list(node_labels.values())[i], # attrs.get("label", node),
                "type": list(node_labels0.values())[i],
                "size": node_sizes[i], #attrs.get("size", 300),
                "color": node_colors[i], # attrs.get("color", "gray"),
                "x": float(scaled_positions[node][0]),
                "y": float(scaled_positions[node][1])
            }
            for i, node in enumerate(G.nodes())  # Correct way to iterate over nodes with attributes
        ],
        "edges": [
            {"from": str(u), "to": str(v), "weight": d.get("weight", 1)}
            for u, v, d in G.edges(data=True)  # Correctly iterating over edges
        ]
        }

        return graph_data

    def _add_source_nodes_G(self, G, db):

        node_ids = list(G.nodes())
        # query = f"""
        #     WITH {node_ids} AS id_list
        #     MATCH (n)-[:BELONGS_TO]-(s:Source)
        #     WHERE n.id IN id_list
        #     RETURN distinct n.id as node_id, s.name as source_name, s.id as source_id
        # """
        # source_nodes = self._query_to_dataframe(db, query)
        # print(source_nodes)
        source_nodes = db.query_source_nodes(node_ids)
        source_nodes_info = source_nodes[['source_name', 'source_id']].drop_duplicates()

        for i, si in source_nodes_info.iterrows():
            G.add_node(si['source_id'], name=si['source_name'], type='source')

        for i, si in source_nodes.iterrows():
            G.add_edge(si['node_id'], si['source_id'], type='BELONS_TO', weight=1)

        return G

    def _query_to_dataframe(self, database, driver, query):
        with driver.session(database=database) as session:
            result = session.run(query)
            records = [record.data() for record in result]
            return pd.DataFrame(records)
        
    def _fa2(self, G):

        # Initialize ForceAtlas2
        forceatlas2 = ForceAtlas2(
            outboundAttractionDistribution=True,  # Dissuade hubs
            linLogMode=False,  
            adjustSizes=False,  
            edgeWeightInfluence=1.0,
            jitterTolerance=1.0,  # Tolerance
            barnesHutOptimize=True,
            barnesHutTheta=1.2,
            scalingRatio=2.0,
            strongGravityMode=False,
            gravity=1.0,
            verbose=True
        )

        # Compute the ForceAtlas2 layout
        pos = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=2000)
        return pos

    def source_graph(self, db, drug, disease):

        data = db.query_attention_pair(disease, drug)#['paths']
        rows = []

        for item in data['paths']:
            row = {
                'nodeIds': [x['nodeId'] for x in item['nodes']],
                'nodeTypes': [x['nodeType'] for x in item['nodes']],
                'edgeInfos': [x['edgeInfo'] for x in item['edges']],
                'edgeScore': [x['score'] for x in item['edges']],
                'avg_score': item['avg_score']
            }
            rows.append(row)


        # Create a DataFrame
        paths_df = pd.DataFrame(rows)

        nodesIds = [item for sublist in paths_df['nodeIds'] for item in sublist]
        nodeTypes = [item for sublist in paths_df['nodeTypes'] for item in sublist] 

        nodes_set = set(zip(nodesIds, nodeTypes))
        Gp = nx.Graph()
        for ni in nodes_set:
            Gp.add_node(ni[0], type=ni[1])

        # make edge list
        edge_lists = []
        for j, pj in enumerate(paths_df['nodeIds']):
            for i in range(0, len(pj)-1):
                edge_lists.append([pj[i], pj[i+1], paths_df['edgeScore'][j][i]])
        for ei in edge_lists:
            Gp.add_edge(ei[0], ei[1], weight=float(ei[2]))

        Gp = self._add_source_nodes_G(Gp, db)
        #positions = nx.kamada_kawai_layout(Gp)
        positions = self._fa2(Gp)

        return self._graph_to_json(Gp, positions)


