'''

disruption_targeting_analysis_functions.py

Functions for running native targeting and proxy targeting approaches for cohesion and disruption analyses

'''

### Importing prerequisites

import networkx as nx
import itertools
import random

### Alternative path centrality function for targeting metric

def alternative_path_centrality(net):
    
    '''
    
    Implementation of Alternative Path Centrality (APC) metric from: Shavitt, Y., & Singer, Y. (2007). Beyond centrality—Classifying topological significance using backup efficiency and alternative 
    paths. New Journal of Physics, 9(8), 266–266. https://doi.org/10.1088/1367-2630/9/8/266
    
    This implementation has currently only been validated for unweighted and undirected networks
    
    Args:
        -net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
    
    Returns:
        -dict: A dictionary with nodes as keys and APC values as values
    
    '''
    
    # Setting up required sub-functions for calculating reciprocal sums
    
    def reciprocal_sum(l):
        return sum(1 / i for i in l if i != 0)
    
    def calculate_reciprocal_sums(net1, lengths_dict, descendant_dict, exclude_node):
        holding_lengths_list_excl_node_endpoints = [lengths_dict[i][j] for i in net1.nodes() for j in descendant_dict[i] 
                                                    if i != exclude_node and j != exclude_node]
        recip_length_sum = reciprocal_sum(holding_lengths_list_excl_node_endpoints)
        
        return recip_length_sum
    
    # Calculating initial dictionaries based on input network
    
    shortest_path_lengths = dict(nx.all_pairs_shortest_path_length(net)) # Dict with nodes in net as keys and a dict as values, with nodes as keys and the path length between the first node and each node as values
    
    net_descendants = {node : nx.descendants(net, node) for node in net.nodes()} # Dict with nodes in net as keys and a list of all reachable nodes as values
    
    # Calculating node-level APC values
    
    centrality = {}
    
    for node in net.nodes():
        
        holding_subgraph = net.copy()
        holding_subgraph.remove_node(node)
        
        subgraph_lengths = dict(nx.all_pairs_shortest_path_length(holding_subgraph))
        
        subgraph_descendants = {node1 : nx.descendants(holding_subgraph, node1) for node1 in holding_subgraph.nodes()}
        
        orig_recip_sum = calculate_reciprocal_sums(net, shortest_path_lengths, net_descendants, node)
        
        subgraph_recip_sum = calculate_reciprocal_sums(holding_subgraph, subgraph_lengths, subgraph_descendants, node)
        
        centrality[node] = orig_recip_sum - subgraph_recip_sum

    return centrality

### Performance metric functions

# Borgatti fragmentation

def borgatti_fragment(net):
    
    '''
    
    Implementation of fragmentation metric from: Borgatti, S. P. (2006). Identifying sets of key players in a social network. Computational and Mathematical Organization Theory, 12(1), 21–34. https://doi.org/10.1007/s10588-006-7084-x

    This implementation has currently only been validated for unweighted and undirected networks
    
    Args:
        -net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
    
    Returns:
        -float: Borgatti fragmentation metric for the input network. Ranges from 0 to 1
    
    '''
    
    node_count = len(net.nodes())
    holding_dist_list = []
    for i in net.nodes():
        for j in net.nodes():
            if nx.has_path(net, i, j) == False:
                pass
            elif i == j:
                pass
            else:
                holding_dist_val = nx.shortest_path_length(net, i, j)
                holding_dist_val_recip = 1 / holding_dist_val
                holding_dist_list.append(holding_dist_val_recip)
    holding_dist_sum = sum(holding_dist_list)
    frag_val = 1 - (holding_dist_sum / (node_count * (node_count - 1)))
    return frag_val

# Average local efficiency

def avg_local_efficiency(net):
    
    '''
    
    Based on the borgatti_fragment function. Average local efficiency is equal to the average of 1 minus the Borgatti fragmentation for each egonet in a given network

    This implementation has currently only been validated for unweighted and undirected networks
    
    Args:
        -net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
    
    Returns:
        -float: Average local efficiency metric for the input network. Ranges from 0 to 1
    
    '''
    
    holding_egonets_list = [nx.ego_graph(net, i) for i in list(net.nodes())]
    holding_local_efficiency_list = []
    for egonet in holding_egonets_list:
        if len(egonet.nodes()) == 1:
            holding_efficiency = 0
        else:
            holding_efficiency = 1 - borgatti_fragment(egonet)
        holding_local_efficiency_list.append(holding_efficiency)
    holding_result = sum(holding_local_efficiency_list) / len(holding_local_efficiency_list)
    return holding_result

### LCC size

def lcc_size(net):
    
    '''
    
    Share of a network's nodes that are part of its largest connected component 

    This implementation has currently only been validated for unweighted and undirected networks
    
    Args:
        -net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
    
    Returns:
        -float: Share of nodes in the input network that are part of the largest connected component. Ranges from 0 to 1
    
    '''
    
    holding_components_list = sorted(nx.connected_components(net), key = len, reverse = True)
    holding_result = len(holding_components_list[0]) / len(net.nodes())
    return holding_result

# 3LCC size

def three_lcc_size(net):
    
    '''
    
    Share of a network's nodes that are part of its three largest connected components 

    This implementation has currently only been validated for unweighted and undirected networks
    
    Args:
        -net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
    
    Returns:
        -float: Share of nodes in the input network that are part of the three largest connected components. Ranges from 0 to 1
    
    '''
    
    holding_components_list = sorted(nx.connected_components(net), key = len, reverse = True)
    index_num = min(3, len(holding_components_list))
    holding_result = 0
    for i in range(0, index_num):
        holding_result = holding_result + len(holding_components_list[i])
    holding_result = holding_result / len(net.nodes())
    return holding_result

### Cohesion analysis functions

# Native cohesion analysis function

def native_network_cohesion(net, path_length, set_size):
    
    '''
    
    For a given input network, path length, and number of nodes, this function finds the set of nodes that maximizes the share of other nodes in the network that are reachable within the path length 
    and returns that share  

    This implementation has currently only been validated for unweighted and undirected networks
    
    Args:
        -net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
        -path_length (int): Maximum number of hops to consider nodes neighbors
        -set_size (int): Number of seed nodes selected
    
    Returns:
        -float: Maximum share of nodes in the input network that are reachable from a set of nodes of size set_size with maximum path length path_length. Ranges from 0 to 1
    
    '''
    
    # Setting up analysis

    holding_net = net.copy()
    nodelist = list(holding_net.nodes())
    percent_reached = 0
    
    # Creating all possible combinations of nodes in net of size set_size
    
    node_combos = list(itertools.combinations(nodelist, set_size))
    
    # Looping through node combinations 
    
    for combo in node_combos:

        # Looping through each node in the node combination
        
        combo_nodes = [combo[i] for i in range(0, len(combo))]
        combo_nodelist = []
        for node in combo:
            
            # Creating induced subgraph on the seed node that includes all neighbors at most path_length away
            
            holding_egonet = nx.ego_graph(holding_net, node, radius = path_length)
            
            # Adding neighbors to combo_nodelist
            
            for alter in list(holding_egonet.nodes()):
                combo_nodelist.append(alter)
                
        # Removing duplicated and seed nodes from combo_nodelist
        
        combo_nodelist_no_combo_nodes_deduped = list(set([i for i in combo_nodelist if i not in combo_nodes]))
        
        # Calcuating the number of reached nodes by the total number of non-seed nodes in the network; replaces current value with new value if new value is higher
        
        holding_percent_reached = len(combo_nodelist_no_combo_nodes_deduped) / (len(nodelist) - set_size)
        if holding_percent_reached > percent_reached:
            percent_reached = holding_percent_reached
        else:
            pass
    
    return percent_reached

# Proxy cohesion analysis function

def proxy_network_cohesion(target_net, cohesion_net, target_func, path_length, set_size):
    
    '''
    
    For a given input targeting network, cohesion network, targeting function, path length, and number of nodes, this function calculates the targeting function on the targeting network, takes the 
    highest-valued nodes, and tests combinations of targeted nodes to find the set of nodes that maximizes the share of other nodes in the cohesion_network that are reachable within the path length and 
    returns that share. Note that this implementation only considers nodes that are shared between the target_net and cohesion_net
    
    This implementation has currently only been validated for unweighted and undirected networks
    
    Args:
        -target_net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
        -cohesion_net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
        -target_func (func): One of the following node-level centrality metric functions - nx.degree_centrality, nx.betweenness_centrality, nx.closeness_centrality, 
        nx.eigenvector_centrality, alternative_path_centrality
        -path_length (int): Maximum number of hops to consider nodes neighbors
        -set_size (int): Number of seed nodes selected
    
    Returns:
        -float: Maximum share of nodes in the input network that are reachable from a set of nodes of size set_size with maximum path length path_length. Ranges from 0 to 1
    
    '''
    
    # Creating common nodeset between target_net and cohesion_net
    
    target_net_nodes = list(target_net.nodes)
    cohesion_net_nodes = list(cohesion_net.nodes)
    common_nodes = list(set([node for node in target_net_nodes if node in cohesion_net_nodes]))
    
    # Initial errors based on values passed in function call
    
    targeting_metric_functions = [nx.degree_centrality, nx.betweenness_centrality, nx.closeness_centrality, 
                                  nx.eigenvector_centrality, alternative_path_centrality]
    
    if target_func not in targeting_metric_functions:
        raise ValueError("Targeting metric must be one of: degree, betweenness, closeness, eigenvector, alternative_path_centrality")
    
    if len(common_nodes) < set_size:
        raise ValueError(f"Number of nodes to be removed must not be more than the total number of common nodes ({len(common_nodes)})")

    # Getting list of targeted nodes by running the target_func on the target_net, only considering shared nodes
        # To account for cases where multiple nodes may have the same value on the targeting metric, this list is constructed by choosing all nodes with the highest value, then all nodes with the 2nd 
        # highest value and so on until the list of targeted nodes is longer than the set_size input. In cases where all nodes have different values for the targeting metric, the top set_size nodes are chosen
    
    target_net_target_metric = target_func(target_net)
    common_nodes_target_metric = {k : target_net_target_metric[k] for k in common_nodes}
    top_x_target_metric_vals = sorted(list(set(common_nodes_target_metric.values())), reverse = True)[:set_size]
    
    targeted_nodes = []
    
    for i in range(0, set_size):
        if len(targeted_nodes) < set_size:
            nodes_w_target_metric_val = [key for key in common_nodes_target_metric 
                                         if common_nodes_target_metric[key] == top_x_target_metric_vals[i]]
            for i in nodes_w_target_metric_val:
                targeted_nodes.append(i)
    
    # Initializing max cohesion value holder

    percent_reached = 0

    # Creating all possible combinations of targeted nodes of size set_size
    
    node_combos = list(itertools.combinations(targeted_nodes, set_size))
    
    # Looping through node combinations
    
    for combo in node_combos:
        
        # Looping through each node in the node combination
        
        combo_nodes = [combo[i] for i in range(0, len(combo))]
        combo_nodelist = []
        for node in combo:
        
            # Creating induced subgraph on the seed node that includes all neighbors at most path_length away

            holding_egonet = nx.ego_graph(cohesion_net, node, radius = path_length)
            
            # Adding neighbors to combo_nodelist
            
            for alter in list(holding_egonet.nodes()):
                combo_nodelist.append(alter)
        
        # Removing duplicated and seed nodes from combo_nodelist
        
        combo_nodelist_no_combo_nodes_deduped = list(set([i for i in combo_nodelist if i not in combo_nodes]))
        
        # Calculating the number of reached nodes in cohesion_net by the total number of non-seed nodes in cohesion_net; replaces current value with new value if new value is higher
        
        holding_percent_reached = len(combo_nodelist_no_combo_nodes_deduped) / (len(cohesion_net_nodes) - set_size)
        if holding_percent_reached > percent_reached:
            percent_reached = holding_percent_reached
        else:
            pass
    
    return percent_reached

### Disruption analysis functions

# Native disruption analysis function

def native_network_disruption(net, target_func, n_remove, seed_list):
    
    '''
    
    Function for running a native targeting disruption analysis. This function iteratively calculates node centralities based on target_func, selects a node to remove based on the highest target_func value,
    and calculates four performance metrics for the new subgraph (Borgatti fragmentation, average local efficiency, LCC size, and 3LCC size). In cases where multiple nodes have the highest target_func value,
    one is chosen at random. This repeats until n_remove nodes have been removed. The function returns a dictionary that contains lists of values for each performance metric for each seed
    
    Note that there is a separate function for when using alternative path centrality (APC) as the targeting metric - see below
    
    This implementation has currently only been validated for unweighted and undirected networks
    
    Args:
        -net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
        -target_func (func): One of the following node-level centrality metric functions - nx.degree_centrality, nx.betweenness_centrality, nx.closeness_centrality, 
        nx.eigenvector_centrality
        -n_remove (int): Maximum number of nodes to be removed from the network as part of the disruption analysis
        -seed_list (list): List of integers to be used as random seeds
    
    Returns:
        -dict: A dictionary with performance metric abbreviations as keys ("BF", "ALE", "LS", "3LS") and lists-of-lists as values, where each sub-list contains the performance metric values for each seeded run
    
    '''
    
    # Initial errors based on values passed in function call
    
    if target_func not in [nx.degree_centrality, nx.betweenness_centrality, nx.closeness_centrality, nx.eigenvector_centrality]:
        raise ValueError("Targeting metric must be one of: degree, betweenness, closeness, eigenvector")
        
    if n_remove > len(net.nodes()) - 2:
        raise ValueError("Number of nodes to be removed must be 2 lower than the total number of nodes in the original network")
    
    if len(seed_list) == 0:
        raise ValueError("Must provide at least one seed")
    
    # Creating aggregate results lists
    
    holding_agg_results_bf = []
    holding_agg_results_ale = []
    holding_agg_results_ls = []
    holding_agg_results_3ls = []
    
    # Looping through the different seeds
    
    for s in seed_list:
        random.seed(s)
        
        # Setting up disruption analysis
        
        holding_net = net.copy()
        removal_counter = 0
        
        # Initializing run's results lists with original network's values for each performance metric 
        
        holding_results_bf = [borgatti_fragment(holding_net)]
        holding_results_ale = [avg_local_efficiency(holding_net)]
        holding_results_ls = [lcc_size(holding_net)]
        holding_results_3ls = [three_lcc_size(holding_net)]
        
        # Network disruption analysis
        
        while removal_counter < n_remove:
            removal_counter += 1
            
            # Calculating node-level targeting metric
            
            if target_func == nx.eigenvector_centrality:
                holding_target_dict = target_func(holding_net, max_iter = 100000)
            else:
                holding_target_dict = target_func(holding_net)
            
            # Getting all nodes in the current version of the network that have the highest value for the targeting metric
            
            holding_max_val = max(holding_target_dict.values())
            max_val_nodes = [key for key in holding_target_dict if holding_target_dict[key] == holding_max_val]
            
            # Where multiple nodes have the maximum value, choose one at random
            
            removal_node = random.choice(max_val_nodes)
            
            # Removing chosen removal_node and calculating the four performance metrics on the new subgraph
            
            holding_net.remove_node(removal_node)
            holding_results_bf.append(borgatti_fragment(holding_net))
            holding_results_ale.append(avg_local_efficiency(holding_net))
            holding_results_ls.append(lcc_size(holding_net))
            holding_results_3ls.append(three_lcc_size(holding_net))
        
        # Appending results lists to aggregate results lists
        
        holding_agg_results_bf.append(holding_results_bf)
        holding_agg_results_ale.append(holding_results_ale)
        holding_agg_results_ls.append(holding_results_ls)
        holding_agg_results_3ls.append(holding_results_3ls)
        
    # Returning dictionary of aggregate results lists
        
    holding_dict = {'BF' : holding_agg_results_bf, 'ALE' : holding_agg_results_ale, 'LS' : holding_agg_results_ls, '3LS' : holding_agg_results_3ls}
        
    return holding_dict

# Native disruption analysis using alternative path centrality (APC) function

def apc_native_network_disruption(net, n_remove, seed_list):
    
    '''
    
    Function for running a native targeting disruption analysis. This function iteratively calculates alternative path centrality values for nodes, selects a node to remove based on these values,
    and calculates four performance metrics for the new subgraph (Borgatti fragmentation, average local efficiency, LCC size, and 3LCC size). In cases where multiple nodes have the highest APC value,
    one is chosen at random. This repeats until n_remove nodes have been removed. The function returns a dictionary that contains lists of values for each performance metric for each seed
    
    Note that this function only uses alternative path centrality (APC) as a targeting metric. This was done to allow for optimizations when updating APC values after removing a node, as APC is a computationally
    expensive function to run
    
    This implementation has currently only been validated for unweighted and undirected networks
    
    Args:
        -net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
        -n_remove (int): Maximum number of nodes to be removed from the network as part of the disruption analysis
        -seed_list (list): List of integers to be used as random seeds
    
    Returns:
        -dict: A dictionary with performance metric abbreviations as keys ("BF", "ALE", "LS", "3LS") and lists-of-lists as values, where each sub-list contains the performance metric values for each seeded run
    
    '''
        
    # Initial errors based on values passed in function call
        
    if n_remove > len(net.nodes()) - 2:
        raise ValueError("Number of nodes to be removed must be 2 lower than the total number of nodes in the original network")
    
    if len(seed_list) == 0:
        raise ValueError("Must provide at least one seed")
    
    # Creating aggregate results lists
    
    holding_agg_results_bf = []
    holding_agg_results_ale = []
    holding_agg_results_ls = []
    holding_agg_results_3ls = []
    
    # Looping through different seeds
    
    for s in seed_list:
        random.seed(s)
        
        # Setting up disruption analysis
        
        holding_net = net.copy()
        apc_dict = alternative_path_centrality(holding_net)
        removal_counter = 0
        
        # Initializing run's results lists with original network's values for each performance metric 
        
        holding_results_bf = [borgatti_fragment(holding_net)]
        holding_results_ale = [avg_local_efficiency(holding_net)]
        holding_results_ls = [lcc_size(holding_net)]
        holding_results_3ls = [three_lcc_size(holding_net)]
        
        # Network disruption analysis
        
        while removal_counter < n_remove:
            removal_counter += 1
            
            # Getting all nodes in the current version of the network that have the highest value for APC based on apc_dict
            
            holding_max_val = max(apc_dict.values())
            max_val_nodes = [key for key in apc_dict if apc_dict[key] == holding_max_val]
            
            # Where multiple nodes have the maximum value, choose one at random
            
            removal_node = random.choice(max_val_nodes)
            
            # Getting a list of all other nodes in the same component as removal_node
            
            removal_node_component = nx.node_connected_component(holding_net, removal_node)
            removal_node_component = [i for i in removal_node_component if i != removal_node]
            
            # Removing chosen removal_node and calculating the four performance metrics on the new subgraph
            
            holding_net.remove_node(removal_node)
            holding_results_bf.append(borgatti_fragment(holding_net))
            holding_results_ale.append(avg_local_efficiency(holding_net))
            holding_results_ls.append(lcc_size(holding_net))
            holding_results_3ls.append(three_lcc_size(holding_net))
            
            # Updating apc_dict values for nodes in the same component as the removed node
                # APC values for nodes in other components are not impacted by changes to the network in other components
                # This reduces the cost of running the alternative_path_centrality function by creating a smaller input network; this has little effect for networks that are comprised of a single component,
                # but materially improves performance as the number of components increases
            
            apc_dict.pop(removal_node)
            holding_subgraph = nx.subgraph(holding_net, removal_node_component)
            holding_subgraph_apc = alternative_path_centrality(holding_subgraph)
            for i in removal_node_component:
                apc_dict[i] = holding_subgraph_apc[i]
        
        # Appending results lists to aggregate results lists
        
        holding_agg_results_bf.append(holding_results_bf)
        holding_agg_results_ale.append(holding_results_ale)
        holding_agg_results_ls.append(holding_results_ls)
        holding_agg_results_3ls.append(holding_results_3ls)
        
    # Returning dictionary of aggregate results lists
        
    holding_dict = {'BF' : holding_agg_results_bf, 'ALE' : holding_agg_results_ale, 'LS' : holding_agg_results_ls, '3LS' : holding_agg_results_3ls}
        
    return holding_dict

# Proxy disruption analysis function

def proxy_network_disruption(target_net, disrupt_net, target_func, n_remove, seed_list):
    
    '''
    
    Function for running a proxy targeting disruption analysis. This function iteratively calculates node centralities based on target_func for target_net, selects a node to remove from both target_net and
    disrupt_net based on the highest target_func value, and calculates four performance metrics for the new subgraph of disrupt_net (Borgatti fragmentation, average local efficiency, LCC size, and 3LCC size). 
    In cases where multiple nodes have the highest target_func value, one is chosen at random. This repeats until n_remove nodes have been removed. The function returns a dictionary that contains lists of 
    values for each performance metric for each seed
    
    Note that there is a separate function for when using alternative path centrality (APC) as the targeting metric - see below
    
    This implementation has currently only been validated for unweighted and undirected networks
    
    Args:
        -target_net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
        -disrupt_net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
        -target_func (func): One of the following node-level centrality metric functions - nx.degree_centrality, nx.betweenness_centrality, nx.closeness_centrality, 
        nx.eigenvector_centrality
        -n_remove (int): Maximum number of nodes to be removed from the network as part of the disruption analysis
        -seed_list (list): List of integers to be used as random seeds
    
    Returns:
        -dict: A dictionary with performance metric abbreviations as keys ("BF", "ALE", "LS", "3LS") and lists-of-lists as values, where each sub-list contains the performance metric values for each seeded run
        
    '''
    
    # Creating common nodeset between target_net and disrupt_net
    
    target_net_nodes = list(target_net.nodes())
    disrupt_net_nodes = list(disrupt_net.nodes())
    common_nodes = list(set([node for node in target_net_nodes if node in disrupt_net_nodes]))
    
    # Initial errors based on values passed in function call

    if target_func not in [nx.degree_centrality, nx.betweenness_centrality, nx.closeness_centrality, nx.eigenvector_centrality]:
        raise ValueError("Targeting metric must be one of: degree, betweenness, closeness, eigenvector")
        
    if n_remove > len(disrupt_net.nodes()) - 2:
        raise ValueError("Number of nodes to be removed must be 2 lower than the total number of nodes in the original network")
    
    if n_remove > len(common_nodes):
        raise ValueError(f"Number of nodes to be removed must not be more than the total number of common nodes ({len(common_nodes)})")
    
    if len(seed_list) == 0:
        raise ValueError("Must provide at least one seed")

    # Creating aggregate results lists
    
    holding_agg_results_bf = []
    holding_agg_results_ale = []
    holding_agg_results_ls = []
    holding_agg_results_3ls = []
    
    # Looping through different seeds
    
    for s in seed_list:
        random.seed(s)
        
        # Setting up disruption analysis
        
        holding_target_net = target_net.copy()
        holding_disrupt_net = disrupt_net.copy()
        removal_counter = 0
        
        # Initializing run's results lists with original disrupt_net's values for each performance metric
        
        holding_results_bf = [borgatti_fragment(holding_disrupt_net)]
        holding_results_ale = [avg_local_efficiency(holding_disrupt_net)]
        holding_results_ls = [lcc_size(holding_disrupt_net)]
        holding_results_3ls = [three_lcc_size(holding_disrupt_net)]
        
        # Network disruption analysis
        
        while removal_counter < n_remove:
            removal_counter += 1
            holding_common_nodes = list(set([node for node in list(holding_target_net.nodes())
                                             if node in list(holding_disrupt_net.nodes())]))
            
            # Calculating node-level targeting metric based on the current version of target_net for common nodes
            
            if target_func == nx.eigenvector_centrality:
                holding_target_dict = target_func(holding_target_net, max_iter = 100000)
            else:
                holding_target_dict = target_func(holding_target_net)
            holding_common_node_target_val_dict = {k : holding_target_dict[k] for k in holding_common_nodes}
            
            # Getting all nodes shared between the current versions of target_net and disrupt_net that have the highest value for the targeting metric
            
            holding_max_val = max(holding_common_node_target_val_dict.values())
            max_val_nodes = [key for key in holding_common_node_target_val_dict
                             if holding_common_node_target_val_dict[key] == holding_max_val]
            
            # Where multiple nodes have the maximum value, choose one at random
            
            removal_node = random.choice(max_val_nodes)
            
            # Removing chosen removal_node from the current versions of target_net and disrupt_net and calculating the four performance metrics on the new disrupt_net subgraph
            
            holding_target_net.remove_node(removal_node)
            holding_disrupt_net.remove_node(removal_node)
            holding_results_bf.append(borgatti_fragment(holding_disrupt_net))
            holding_results_ale.append(avg_local_efficiency(holding_disrupt_net))
            holding_results_ls.append(lcc_size(holding_disrupt_net))
            holding_results_3ls.append(three_lcc_size(holding_disrupt_net))
       
        # Appending results lists to aggregate results lists
        
        holding_agg_results_bf.append(holding_results_bf)
        holding_agg_results_ale.append(holding_results_ale)
        holding_agg_results_ls.append(holding_results_ls)
        holding_agg_results_3ls.append(holding_results_3ls)
        
    # Returning dictionary of aggregate results lists
        
    holding_dict = {'BF' : holding_agg_results_bf, 'ALE' : holding_agg_results_ale, 'LS' : holding_agg_results_ls, '3LS' : holding_agg_results_3ls}
        
    return holding_dict

# Proxy disruption analysis using alternative path centrality (APC) function

def apc_proxy_network_disruption(target_net, disrupt_net, n_remove, seed_list):
    
    '''
    
    Function for running a proxy targeting disruption analysis. This function iteratively calculates alternative path centrality values for target_net, selects a node to remove from both target_net and
    disrupt_net based on the highest APC value, and calculates four performance metrics for the new subgraph of disrupt_net (Borgatti fragmentation, average local efficiency, LCC size, and 3LCC size). 
    In cases where multiple nodes have the highest APC value, one is chosen at random. This repeats until n_remove nodes have been removed. The function returns a dictionary that contains lists of 
    values for each performance metric for each seed
    
    Note that this function only uses alternative path centrality (APC) as a targeting metric. This was done to allow for optimizations when updating APC values after removing a node, as APC is a computationally
    expensive function to run
    
    This implementation has currently only been validated for unweighted and undirected networks
    
    Args:
        -target_net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
        -disrupt_net (networkx.Graph): A NetworkX graph object. The current implementation has only been validated for unweighted and undirected networks
        -n_remove (int): Maximum number of nodes to be removed from the network as part of the disruption analysis
        -seed_list (list): List of integers to be used as random seeds
    
    Returns:
        -dict: A dictionary with performance metric abbreviations as keys ("BF", "ALE", "LS", "3LS") and lists-of-lists as values, where each sub-list contains the performance metric values for each seeded run

    
    '''
    
    # Creating common nodeset between target_net and disrupt_net
    
    target_net_nodes = list(target_net.nodes())
    disrupt_net_nodes = list(disrupt_net.nodes())
    common_nodes = list(set([node for node in target_net_nodes if node in disrupt_net_nodes]))
    
    # Initial errors based on values passed in function call
        
    if n_remove > len(disrupt_net.nodes()) - 2:
        raise ValueError("Number of nodes to be removed must be 2 lower than the total number of nodes in the original network")
    
    if n_remove > len(common_nodes):
        raise ValueError(f"Number of nodes to be removed must not be more than the total number of common nodes ({len(common_nodes)})")
    
    if len(seed_list) == 0:
        raise ValueError("Must provide at least one seed")

    # Creating aggregate results lists
    
    holding_agg_results_bf = []
    holding_agg_results_ale = []
    holding_agg_results_ls = []
    holding_agg_results_3ls = []
    
    # Looping through different seeds
    
    for s in seed_list:
        random.seed(s)
        
        # Setting up disruption analysis
        
        holding_target_net = target_net.copy()
        holding_disrupt_net = disrupt_net.copy()
        target_apc_dict = alternative_path_centrality(holding_target_net)
        removal_counter = 0
        
        # Initializing run's results lists with original disrupt_net's values for each performance metric
        
        holding_results_bf = [borgatti_fragment(holding_disrupt_net)]
        holding_results_ale = [avg_local_efficiency(holding_disrupt_net)]
        holding_results_ls = [lcc_size(holding_disrupt_net)]
        holding_results_3ls = [three_lcc_size(holding_disrupt_net)]
        
        # Network disruption analysis
        
        while removal_counter < n_remove:
            removal_counter += 1
            holding_common_nodes = list(set([node for node in list(holding_target_net.nodes())
                                             if node in list(holding_disrupt_net.nodes())]))
            
            # Getting APC values based on target_net for common nodes
            
            holding_common_node_apc_val_dict = {k : target_apc_dict[k] for k in holding_common_nodes}
            
            # Getting all nodes shared between the current versions of target_net and disrupt_net that have the highest APC value
            
            holding_max_val = max(holding_common_node_apc_val_dict.values())
            max_val_nodes = [key for key in holding_common_node_apc_val_dict
                             if holding_common_node_apc_val_dict[key] == holding_max_val]
            
            # Where multiple nodes have the maximum value, choose one at random
            
            removal_node = random.choice(max_val_nodes)
            
            # Getting a list of all other nodes in the same component as removal_node
            
            removal_node_component_target_net = nx.node_connected_component(holding_target_net, removal_node)
            removal_node_component_target_net = [i for i in removal_node_component_target_net if i != removal_node]
            
            # Removing chosen removal_node from the current versions of target_net and disrupt_net and calculating the four performance metrics on the new disrupt_net subgraph
            
            holding_target_net.remove_node(removal_node)
            holding_disrupt_net.remove_node(removal_node)
            holding_results_bf.append(borgatti_fragment(holding_disrupt_net))
            holding_results_ale.append(avg_local_efficiency(holding_disrupt_net))
            holding_results_ls.append(lcc_size(holding_disrupt_net))
            holding_results_3ls.append(three_lcc_size(holding_disrupt_net))
            
            # Updating target_apc_dict values for nodes in the same component as the removed node
                # APC values for nodes in other components are not impacted by changes to the network in other components
                # This reduces the cost of running the alternative_path_centrality function by creating a smaller input network; this has little effect for networks that are comprised of a single component,
                # but materially improves performance as the number of components increases
            
            target_apc_dict.pop(removal_node)
            holding_subgraph = nx.subgraph(holding_target_net, removal_node_component_target_net)
            holding_subgraph_apc = alternative_path_centrality(holding_subgraph)
            for i in removal_node_component_target_net:
                target_apc_dict[i] = holding_subgraph_apc[i]
       
        # Appending results lists to aggregate results lists
        
        holding_agg_results_bf.append(holding_results_bf)
        holding_agg_results_ale.append(holding_results_ale)
        holding_agg_results_ls.append(holding_results_ls)
        holding_agg_results_3ls.append(holding_results_3ls)
        
    # Returning dictionary of aggregate results lists
        
    holding_dict = {'BF' : holding_agg_results_bf, 'ALE' : holding_agg_results_ale, 'LS' : holding_agg_results_ls, '3LS' : holding_agg_results_3ls}
        
    return holding_dict