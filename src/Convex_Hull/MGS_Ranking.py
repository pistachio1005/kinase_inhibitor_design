import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import animation
import glob
import random
import math
import itertools
import copy

import numpy as np

from resources.lookup_dicts import SingleAAtoTripleAA, GraphColors


def get_doublet_scores(logfile: str, score_cutoff: float):

    doublet_scores = []

    # pair_set_kstar/match*-kstar/kstar-[*_*]/submit.out
    kname = logfile.split('/')[2].split('-')[1].split('_')
    resnums = []
    for i in kname:
        if '[' in i:
            new = i.replace('[', '')
            resnums.append(new)
        elif ']' in i:
            new = i.replace(']', '')
            resnums.append(new)

    on_header = False

    f = open(logfile, 'r')
    for line in f:
        if "error" in line or "Error" in line or "Java" in line or "java" in line:
            print("ERROR! Found logfile issue with %s" % logfile)
            exit()
        if "seconds" in line:
            break
        if on_header:
            fields = line.split(',')

            # get the kstar score
            score = fields[2]
            if score == 'none' or '-' in score:
                continue
            score = float(score)
            if score <= score_cutoff:
                continue

            # get the mutations in order
            mut_field = fields[1].split(' ')
            mut_index = []
            counter = 0
            for item in mut_field:
                if '=' in item:
                    counter += 1
                elif item == resnums[0]:
                    counter += 1
                    mut_index.append(counter)
                else:
                    counter += 1
            counter = 0
            for item in mut_field:
                if '=' in item:
                    counter += 1
                elif item == resnums[1]:
                    counter += 1
                    mut_index.append(counter)
                else:
                    counter += 1

            # convert the mutations to single letter
            mut1 = mut_field[mut_index[0]].split('=')[1]
            mut2 = mut_field[mut_index[1]].split('=')[1]

            mut1_single = SingleAAtoTripleAA[mut1]
            mut2_single = SingleAAtoTripleAA[mut2]

            new_entry = (resnums[0], resnums[1], mut1_single, mut2_single, score)

            doublet_scores.append(new_entry)

        if "Assignments" in line:
            on_header = True

    ordered_scores = sorted(doublet_scores, key=lambda tup: tup[4], reverse=True)

    return ordered_scores


def kstar_to_g(log_kstar):

    raw_kstar = pow(10, log_kstar)

    dG = -1 * 8.3145 * 298 * math.log(raw_kstar)

    kcal = dG * (1 / 4184)

    return kcal


def energy_window_doublets(log_scores, window_cutoff: float):

    window_scores = []

    best_node_dG = kstar_to_g(log_scores[0][4])

    for node in log_scores:

        curr_dG = kstar_to_g(node[4])

        ddG = abs(best_node_dG - curr_dG)

        if ddG < window_cutoff:

            window_scores.append(node)

    ordered_scores = sorted(window_scores, key=lambda tup: tup[4], reverse=False)

    return ordered_scores


def add_all_nodes(scores, graph):
    print("--- Adding nodes ---")
    used_colors = [data["color"] for v, data in graph.nodes(data=True)]
    node_color = random.choice(GraphColors)
    while node_color in used_colors:
        node_color = random.choice(GraphColors)
    for node in scores:
        ID = ("%s%s%s%s" % (node[0], node[2], node[1], node[3]))
        graph.add_node(ID, identifier=ID, doublet=[node[0], node[1]], resmuts={node[0]: node[2], node[1]: node[3]}, score=node[4], color=node_color)


def need_doublet_edge(n1_node, n2_node):

    n1_doublet = n1_node["doublet"]
    n2_doublet = n2_node["doublet"]

    shared_res = [x for x in n1_doublet if x in n2_doublet]

    # doublets with no shared index are always compatible
    if len(shared_res) == 0:
        return True

    # if do share index, check AA type
    if len(shared_res) == 1:
        shared_res = shared_res[0]
        n1_share = n1_node["resmuts"][shared_res]
        n2_share = n2_node["resmuts"][shared_res]
        if n1_share == n2_share:
            return True

    return False


def add_all_edges(graph):

    print("--- Adding edges ---")
    all_nodes = list(graph.nodes)

    for n1 in range(0, len(all_nodes)-1):
        for n2 in range(n1+1, len(all_nodes)):
            n1_node = graph.nodes[all_nodes[n1]]
            n2_node = graph.nodes[all_nodes[n2]]
            if need_doublet_edge(n1_node, n2_node):
                graph.add_edge(n1_node["identifier"], n2_node["identifier"])


def visualize_MGS(graph, visual_type: str):
    layers = {}
    all_doublets = []

    for node in graph.nodes():
        doub = graph.nodes[node]["doublet"]
        if doub not in all_doublets:
            all_doublets.insert(0, doub)

    counter = 1
    for d in all_doublets:
        related_nodes = []
        for node in graph.nodes():
            doub = graph.nodes[node]["doublet"]
            if d == doub:
                related_nodes.append(node)
        layers[counter] = related_nodes
        counter += 1

    column_colors = [data["color"] for v, data in graph.nodes(data=True)]

    if visual_type == '2D':
        pos = nx.multipartite_layout(graph, subset_key=layers)
        plt.figure(figsize=(12, 12))
        nx.draw(graph, pos, with_labels=True, node_color=column_colors, font_size=6)
        plt.axis("equal")
        plt.show()
        return

    elif visual_type == '3D':
        pos_2d = nx.multipartite_layout(graph, subset_key=layers)
        pos_3d = {}
        for id, coord in pos_2d.items():
            descr = [float(coord[0]), 0, float(graph.nodes[id]['score'])]
            array = np.array(descr)
            pos_3d[id] = array

        edges = np.array([(pos_3d[u], pos_3d[v]) for u, v in graph.edges()])

        fig = plt.figure()

        ax = fig.add_subplot(111, projection="3d")
        for id, coord in pos_3d.items():
            ax.scatter(coord[0], coord[1], coord[2], alpha=0.8, color=graph.nodes[id]['color'], s=20)

        for vizedge in edges:
            ax.plot(*vizedge.T, color="gray")

        for node in graph.nodes(data=True):
            muts = list(node[1]['resmuts'].values())
            node_label = ''.join(muts)
            ax.text(*pos_3d[node[0]], node_label, size=10)

        for members in layers.values():
            top_member_ID = members[len(members)-1]
            top_node = graph.nodes[top_member_ID]
            column_name = '_'.join(top_node['doublet'])
            coords = pos_3d[top_member_ID]
            title_coords = [coords[0], coords[1], coords[2]+3]
            ax.text(*title_coords, column_name)

        ax.grid(False)
        for dim in (ax.xaxis, ax.yaxis, ax.zaxis):
            dim.set_ticks([])
        ax.set_xlabel('Pair Sets')
        ax.set_zlabel('K* Score')

        plt.tight_layout()
        plt.show()
        return

    elif visual_type == '3D_animated':
        def init():
            pos_2d = nx.multipartite_layout(graph, subset_key=layers)
            pos_3d = {}
            for id, coord in pos_2d.items():
                descr = [float(coord[0]), 0, float(graph.nodes[id]['score'])]
                array = np.array(descr)
                pos_3d[id] = array

            edges = np.array([(pos_3d[u], pos_3d[v]) for u, v in graph.edges()])

            for id, coord in pos_3d.items():
                ax.scatter(coord[0], coord[1], coord[2], alpha=0.8, color=graph.nodes[id]['color'], s=20)

            for vizedge in edges:
                ax.plot(*vizedge.T, color="gray")

            for node in graph.nodes(data=True):
                muts = list(node[1]['resmuts'].values())
                node_label = ''.join(muts)
                ax.text(*pos_3d[node[0]], node_label, size=5)

            for members in layers.values():
                top_member_ID = members[len(members)-1]
                top_node = graph.nodes[top_member_ID]
                column_name = '_'.join(top_node['doublet'])
                coords = pos_3d[top_member_ID]
                title_coords = [coords[0], coords[1], coords[2]+3]
                ax.text(*title_coords, column_name)

            ax.grid(False)
            for dim in (ax.xaxis, ax.yaxis, ax.zaxis):
                dim.set_ticks([])
            ax.set_xlabel('Pair Sets')
            ax.set_zlabel('K* Score')

            plt.tight_layout()
            return

        def _frame_update(index):
            ax.view_init(index * 0.2, index * 0.5)
            return

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")

        ani = animation.FuncAnimation(
            fig,
            _frame_update,
            init_func=init,
            interval=100,
            cache_frame_data=False,
            frames=100,
        )
        plt.show()


def is_valid_clique(graph, path: list, tracker: dict):

    for node in path:
        doublet_list = graph.nodes()[node]['doublet']
        doublet_key = ("%s_%s" % (doublet_list[0], doublet_list[1]))
        if tracker[doublet_key]:
            return False
        else:
            tracker[doublet_key] = True

    for status in tracker.values():
        if not status:
            return False

    return True

def clique_to_sequence(graph, path: list):

    chain_length = 0

    for node in path:
        doublet_list = graph.nodes()[node]['doublet']
        if int(doublet_list[0]) > chain_length:
            chain_length = int(doublet_list[0])
        if int(doublet_list[1]) > chain_length:
            chain_length = int(doublet_list[1])

    sequence = [''] * chain_length
    total_dG = 0

    for node in path:
        mutations = graph.nodes()[node]['resmuts']
        for index, identity in mutations.items():
            new_index = int(index) - 1
            sequence[new_index] = identity
        score = graph.nodes()[node]['score']
        dG = kstar_to_g(score)
        total_dG += dG

    return sequence, total_dG


# faster if you have a small # nodes
def find_sequences_from_all_cliques(full_graph):

    print("Finding all fully connected subgraphs...")

    columns = []

    for node in full_graph.nodes():
        curr_doublet = full_graph.nodes[node]['doublet']
        if curr_doublet not in columns:
            columns.insert(0, curr_doublet)

    number_partitions = len(columns)
    number_nodes = len(list(full_graph.nodes))
    number_edges = len(list(full_graph.edges))

    print("This graph has %s partitions, %s nodes, and %s edges" % (number_partitions, number_nodes, number_edges))

    visit_column = {}

    for curr_column in columns:
        new_key = ("%s_%s" % (curr_column[0], curr_column[1]))
        visit_column[new_key] = False

    print("Searching for all cliques...")

    all_subgraphs = [c for c in nx.algorithms.clique.find_cliques(full_graph) if len(c) == len(columns)]

    print('Found %s cliques, verifying each sequence...' % len(all_subgraphs))

    total_clique = 0

    all_sequences = {}

    for g in all_subgraphs:
        curr_visits = copy.deepcopy(visit_column)
        if is_valid_clique(full_graph, g, curr_visits):
            total_clique += 1
            new_seq, energy = clique_to_sequence(full_graph, g)
            all_sequences[tuple(new_seq)] = energy

    print("Found %s valid sequences:" % total_clique)
    sorted_sequences = dict((sorted(all_sequences.items(), key=lambda item: item[1])))
    for seq, score in sorted_sequences.items():
        print(seq, score)


class Partition:

    def __init__(self, posns, member_nodes):
        self.posns = posns
        self.member_nodes = member_nodes

    def __str__(self):
        return "Pos: %s, Member Nodes: %s" % (self.posns, self.member_nodes)


def node_compatibility(node1, node2):

    shared_residues = set(node1[1]['doublet']) & set(node2[1]['doublet'])

    if not shared_residues:
        return True

    shared_residue = list(shared_residues)[0]

    mut1 = node1[1]['resmuts'][shared_residue]
    mut2 = node2[1]['resmuts'][shared_residue]

    return mut1 == mut2


def recursive_DFS(pinfo: list, nodes):

    # base case
    if len(nodes) == len(pinfo):
        return [nodes]

    sequences = []

    current_partition = pinfo[len(nodes)]

    for node in current_partition.member_nodes:
        is_compatible = True
        for path_node in nodes:
            if not node_compatibility(node, path_node):
                is_compatible = False
                break
        if is_compatible:
            new_path = (*nodes, node)
            sequences += recursive_DFS(pinfo, new_path)

    return sequences


def path_to_sequence(path, chain_length):

    sequence = [''] * chain_length
    total_dG = 0

    for node in path:
        mutations = node[1]['resmuts']
        for index, identity in mutations.items():
            new_index = int(index) - 1
            sequence[new_index] = identity
        score = node[1]['score']
        dG = kstar_to_g(score)
        total_dG += dG

    return sequence, total_dG


def find_sequences(graph, chain_length, outfile):

    partitions = set()
    for node in graph.nodes(data=True):
        partitions.add(tuple(node[1]['doublet']))

    partition_info = []
    for p in partitions:
        new_partition = Partition(p, [x for x in graph.nodes(data=True) if x[1]['doublet'] == list(p)])
        partition_info.append(new_partition)

    # future work: order partitions for better pruning + create mappings between partitions for quick lookup

    print("Searching for full sequences...")
    full_paths = recursive_DFS(partition_info, [])

    print("Found %s full sequences- now sorting by energy" % len(full_paths))
    full_sequences = []
    for p in full_paths:
        new_seq = path_to_sequence(p, chain_length)
        full_sequences.append(new_seq)

    print("Full sequences:")
    full_sequences.sort(key=lambda x: x[1])
    for f in full_sequences:
        print(f)

    with open(outfile, 'w') as file:
        for n in range(1, chain_length+1):
            file.write('Res' + str(n) + ',')
        file.write("deltaG\n")
        for seq in full_sequences:
            for amino in seq[0]:
                file.write(str(amino) + ',')
            file.write(str(seq[1]) + '\n')
    print("Saved sequences to %s" % outfile)


# create an empty graph
doublet_graph = nx.Graph()

# visualize doublets for each match (example code shown for 1 match)
for logfile in glob.glob("pair_set_kstar/match1-kstar/kstar-*/submit.out"):
    print("Getting K* results from %s" % logfile)

    # get the Pair Set K* Scores from the log files (submit.out)
    # optional: ignore nodes below a user-specified K* score
    doublet_scores = get_doublet_scores(logfile, 0)

    # reduce nodes via energy window (kcal/mol)
    # if you don't want a window, set this cutoff to math.inf
    window_scores = energy_window_doublets(doublet_scores, math.inf)

    # construct nodes into an MG
    add_all_nodes(window_scores, doublet_graph)

# recommended: for many designs, a recursive tree search of ordered partitions is needed
# input the networkx graph and the length of design chain
find_sequences(doublet_graph, 18, 'math1-fullseqs.csv')

# alternative: for small chains, brute force using networkx methods is ok
# add_all_edges(doublet_graph)
# find_sequences_from_all_cliques(doublet_graph)

# visualize the search space. Animated visuals must be run on command line (intellij does not support)
# options: 2D, 3D, 3D_animated
# visualize_MGS(doublet_graph, '2D')
