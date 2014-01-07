# Copyright (c) 2007 Pedro Matiello <pmatiello@gmail.com>
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.


"""
python-graph

A library for working with graphs in Python.
"""


# Module metadata
__authors__ = "Pedro Matiello"
__license__ = "MIT"


# Imports
import searching
import accessibility
import sorting
import minmax


class graph:
	"""
	Graph class.

	Basic operations are defined in this file.
	Algorithms should refer to external files.
	"""


	def __init__(self):
		"""
		Initialize a graph.
		"""
		self.nodes = {}		# Arrow/Edge lists	(like an adjacency list)
		self.weights = {}	# Arrow/Edge weight list


	def __str__(self):
		"""
		Return a string representing the graph when requested by str() (or print).

		@rtype:  string
		@return: String representing the graph.
		"""
		return "<graph object " + str(self.get_nodes()) + " " + str(self.weights) + ">"


	def __len__(self):
		"""
		Return the size of the graph when requested by len().

		@rtype:  number
		@return: Size of the graph.
		"""
		return len(self.nodes)


	def get_nodes(self):
		"""
		Return node list.

		@rtype:  list
		@return: Node list.
		"""
		return self.nodes.keys()


	def get_node(self, node):
		"""
		Return requested node.

		@type  node: *
		@param node: Node identifier

		@rtype:  list
		@return: List of nodes directly accessible from given node.
		"""
		return self.nodes[node]


	def add_nodes(self, nodelist):
		"""
		Create num nodes.

		@type  nodelist: list
		@param nodelist: List of nodes to be added to the graph.
		"""
		for each in nodelist:
			self.nodes[each] = []


	def add_edge(self, u, v, wt=1):
		"""
		Add an edge (u,v) to the graph connecting nodes u and v.

		@attention: this function should not be used in directed graphs: use add_arrow() instead.

		@type  u: *
		@param u: One node.

		@type  v: *
		@param v: Other node.
		
		@type  wt: number
		@param wt: Edge weight.
		
		"""
		if (v not in self.nodes[u]):
			self.nodes[u].append(v)
			self.nodes[v].append(u)
			self.weights[(u, v)] = wt
			self.weights[(v, u)] = wt


	def add_arrow(self, u, v, wt=1):
		"""
		Add an arrow (u,v) to the directed graph connecting node u to node v.

		@type  u: *
		@param u: One node.

		@type  v: *
		@param v: Other node.

		@type  wt: number
		@param wt: Arrow weight.
		"""
		if (v not in self.nodes[u]):
			self.nodes[u].append(v)
			self.weights[(u, v)] = wt


	def del_edge(self, u, v):
		"""
		Remove an edge (u, v) from the graph.

		@attention: this function should not be used in directed graphs: use del_arrow() instead.

		@type  u: *
		@param u: One node.

		@type  v: *
		@param v: Other node.
		"""
		if (v in self.nodes[u]):
			self.nodes[u].remove(v)
			self.nodes[v].remove(u)
			del(self.weights[(u,v)])
			del(self.weights[(v,u)])


	def del_arrow(self, u, v):
		"""
		Remove an arrow (u, v) from the directed graph.

		@type  u: *
		@param u: One node.

		@type  v: *
		@param v: Other node.
		"""
		if (v in self.nodes[u]):
			self.nodes[u].remove(v)
			del(self.weights[(u,v)])


	def get_arrow_weight(self, u, v):
		"""
		Get the weight of an arrow.

		@type  u: *
		@param u: One node.

		@type  v: *
		@param v: Other node.
		
		@rtype:  number
		@return: Arrow weight
		"""
		return(self.weights[(u, v)])


	def get_edge_weight(self, u, v):
		"""
		Get the weight of an arrow.

		@type  u: *
		@param u: One node.

		@type  v: *
		@param v: Other node.
		
		@rtype:  number
		@return: Edge weight
		"""
		return(self.weights[(u, v)])


	def depth_first_search(self):
		"""
		Depht-first search.

		@rtype:  tuple
		@return:  tupple containing a dictionary and two lists:
			1. Generated spanning tree
			2. Graph's preordering
			3. Graph's postordering
		"""
		return searching.depth_first_search(self)


	def breadth_first_search(self):
		"""
		Breadth-first search.

		@rtype:  dictionary
		@return: Generated spanning_tree
		"""
		return searching.breadth_first_search(self)


	def accessibility(self):
		"""
		Accessibility matrix (transitive closure).

		@rtype:  dictionary
		@return: Accessibility information for each node.
		"""
		return accessibility.accessibility(self)


	def mutual_accessibility(self):
		"""
		Mutual-accessibility matrix (strongly connected components).

		@rtype:  list
		@return: Mutual-accessibility information for each node.
		"""
		return accessibility.mutual_accessibility(self)


	def topological_sorting(self):
		"""
		Topological sorting.

		@attention: Topological sorting is meaningful only for directed acyclic graphs.

		@rtype:  list
		@return: Topological sorting for the graph.
		"""
		return sorting.topological_sorting(self)


	def connected_components(self):
		"""
		Connected components.

		@attention: Indentification of connected components is meaningful only for non-directed graphs.

		@rtype:  dictionary
		@return: Pairing that associates each node to its connected component.
		"""
		return accessibility.connected_components(self)

	def minimal_spanning_tree(self):
		"""
		Minimal spanning tree.

		@attention: Minimal spanning tree meaningful only for weighted graphs.

		@type  graph: graph
		@param graph: Graph.

		@rtype:  list
		@return: Generated spanning tree.
		"""
		return minmax.minimal_spanning_tree(self)

    ###############################################################
    ##### from here on, AvB added code ############################
    ##### from here on, AvB added code ############################
    ###############################################################

    def get_node_weights(self,node):
        """
        Returns list of weights of edges connected to a node

		@type  node: *
		@param node: One Node

		@rtype:  list
		@return: List of nodes directly accessible from given node.
        """
        return [ self.weights[(u,node)] for u in self.get_node(node) ]

    # end of function get_node_weights


    def get_node_weighted_connectivity(self,node):
        """
        Returns list of weights of edges connected to a node

		@type  node: *
		@param node: Node Identifier

		@rtype:  float
		@return: relative weight of outgoing edges
        """
        return float(sum(self.get_node_weights(node))) /\
                 ( float(sum(self.weights.values()))/2.0 )

    # end of function get_node_weighted_connectivity


    def weakest_connected_node(self):
        """
        Return (relatively) weakest connected node

        @attention: only usefull if weights != 1 are assigned

		@rtype:  number
		@return: Node number
        """
        tmp = []
        for node in self.get_nodes():
            tmp.append( ( self.get_node_weighted_connectivity(node), node ) )
        tmp.sort()
        return tmp[0][1]

    # end of function weakest_connected_node


	def add_node(self, node):
		"""
		Create a single node.

		@type  node: *
		@param node: One Node identifier
		"""
        self.nodes[node] = []

    # end of function add_node


	def del_node(self, node):
        """
        Remove a node from the graph.

		@attention: this function removes all concerning edges, arrows and weigths.

		@type  node: *
		@param node: One Node Identifier
        """
		if node in self.get_nodes():
			# delete all edges, weights and outward
			# pointing arrows from this node
			for u in self.nodes[node]:
				self.nodes[u].remove(node)
				try:
					del(self.weights[(node,u)])
					del(self.weights[(u,node)])
				except:
					pass
			# scan all nodes to see if there is an
			# inwards pointing arrow
			for u in self.get_nodes():
				if node in self.nodes[u]:
					self.nodes[u].remove(node)
					try:    self.del_arrow(u,node)
					except: pass
            # and finally delete the node n it self
            del( self.nodes[node] )

    # end of function del_node


    def has_edge(self, u, v):
		"""
		Does edge (u,v) connecting nodes u and v exist?

		@type  u: *
		@param u: One node Identifier.

		@type  v: *
		@param v: Other node Identifier.
		"""
        if (v in self.nodes[u]) and (u in self.nodes[v]):
            if self.weights.has_key((u, v)) and self.weights.has_key((v, u)):
                return True
            else:
                return False
        else:
            return False

    # end of function has_edge


	def set_edge_weight(self, u, v, w):
		"""
		Set the weight of an edge.

		@type  u: *
		@param u: One node Identifier.

		@type  v: *
		@param v: Other node Identifier.

		@type  w: number (integer or float)
		@param w: Weight of this edge.
		"""
		self.weights[(u, v)] = w
		self.weights[(v, u)] = w

    # end of function set_edge_weight


    def get_subgraph_count(self):
        """
        Get number of non-connected subgraphs.
        """
        return len( Set( self.connected_components() ) )

    # end of function get_subgraph_count


    def node_count(self):
        """ """
        return len(self)

    # end of function node_count


    def edge_count(self):
        """
        Get number of edges in the graph.

		@attention: Function is only meaningfull for graphs without arrows.

		@rtype:  number
		@return: Number of edges in the graph
        """
        return len(self.weights)/2

    # end of function edge_count


    def max_edge_count(self):
        """
        Get maximal number of edges possible.

		@rtype:  number
		@return: Maximum number of possible edges
        """
        return sum(range(0,len(self)))

    # end of function max_edge_count


    def connectivitysaturation(self):
        """
        Get relative connectivity saturation.

		@rtype:  float
		@return: ratio of maximum number of possible edges
        """
        return float( self.edge_count() ) / float( self.max_edge_count() )

    # end of function connectivitysaturation


    def degree(self,u):
        """
        Return the degree of node u.

		@type  u: *
		@param u: One node Identifier.

		@rtype:  number
		@return: Number of edges connected to this node.
        """
        return len(self.get_node(u))

    # end of function get_max_edge_count


    def node_connectivity(self):
        """
        Count for each node the number of edges and summarize how
        many nodes have a specific number of edges.

		@rtype:  tuple
		@return: A tuple containing a count number of edges per node
        """
        connectivity = [0]*len(self)
        for node in self.get_nodes():
            connectivity[ len(self.get_node(node)) ] += 1
        return tuple(connectivity)

    # end of function node_connectivity
