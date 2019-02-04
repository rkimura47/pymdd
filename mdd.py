from itertools import chain # used in various places
from collections import deque # used in prune_recursive
from json import dump, load # used in dumpJSON and loadJSON
from random import sample # used in compile_top_down

class MDDArc(object):
    """MDDArc represents a single arc in the MDD.

    MDDArc represents a single arc in the MDD.  An MDDArc is uniquely
    identified by its head/tail nodes, label, and weight.

    Args:
        label (object): label of arc (e.g., assigned value)
        weight (float): weight of arc (e.g., coefficient)
        tail (MDDNode): tail/source node
        head (MDDNode): head/destination node
    """

    def __init__(self, label, weight, tail, head):
        """Construct a new 'MDDArc' object."""
        self.label = label
        self.weight = weight
        self.tail = tail
        self.head = head

    # Allows MDDArcs to be used as dictionary keys.
    def __hash__(self):
        """Return the hash value."""
        return hash((self.label, self.tail, self.head))

    # Rich comparison methods: here the latter four are automatically
    # derived from the first two
    def __eq__(self, other):
        """Return self == other."""
        return self.tail == other.tail and self.label == other.label and self.head == other.head and self.weight == other.weight
    def __lt__(self, other):
        """Return self < other."""
        if self.tail != other.tail:
            return self.tail < other.tail
        elif self.label != other.label:
            return self.label < other.label
        elif self.head != other.head:
            return self.head < other.head
        else:
            return self.weight < other.weight
    def __ne__(self, other):
        """Return self != other."""
        return not self.__eq__(other)
    def __gt__(self, other):
        """Return self > other."""
        return not(self.__eq__(other) or self.__lt__(other))
    def __le__(self, other):
        """Return self <= other."""
        return self.__eq__(other) or self.__lt__(other)
    def __ge__(self, other):
        """Return self >= other."""
        return self.__eq__(other) or not self.__lt__(other)

    def __str__(self):
        return 'A(' + str(self.label) + ',' + str(self.weight) + ':' + str(self.tail) + ',' + str(self.head) + ')'

    def __repr__(self):
        s = 'MDDArc('
        s += repr(self.label) + ', '
        s += repr(self.weight) + ', '
        s += repr(self.tail) + ', '
        s += repr(self.head) + ')'
        return s

class MDDNode(object):
    """MDDNode represents a single node in the MDD.

    MDDNode represents a single node in the MDD.  An MDDNode is uniquely
    identified by its layer and state.  The (node) state must be a hashable
    object.

    Args:
        layer (int): layer the node is in
        state (object): state associated with node
    """

    def __init__(self, layer, state):
        """Construct a new 'MDDNode' object."""
        self.layer = layer
        self.state = state

    # Allows MDDNodes to be used as dictionary keys.
    def __hash__(self):
        """Return the hash value."""
        return hash((self.layer, self.state))

    # Rich comparison methods: here the latter four are automatically
    # derived from the first two
    def __eq__(self, other):
        """Return self == other."""
        return self.layer == other.layer and self.state == other.state
    def __lt__(self, other):
        """Return self < other."""
        if self.layer != other.layer:
            return self.layer < other.layer
        else:
            return self.state < other.state
    def __ne__(self, other):
        """Return self != other."""
        return not self.__eq__(other)
    def __gt__(self, other):
        """Return self > other."""
        return not(self.__eq__(other) or self.__lt__(other))
    def __le__(self, other):
        """Return self <= other."""
        return self.__eq__(other) or self.__lt__(other)
    def __ge__(self, other):
        """Return self >= other."""
        return self.__eq__(other) or not self.__lt__(other)

    def __str__(self):
        return 'N_' + str(self.layer) + '(' + str(self.state) + ')'

    def __repr__(self):
        return 'MDDNode(' + repr(self.layer) + ', ' + repr(self.state) + ')'

class MDDNodeInfo(object):
    """MDDNodeInfo represents information associated with an MDDNode.

    MDDNodeInfo represents information associated with an MDDNode.

    Args:
        incoming (set): set of incoming arcs (default: set())
        outgoing (set): set of outgoing arcs (default: set())
    """

    def __init__(self, incoming=None, outgoing=None):
        """Construct a new 'MDDNode' object."""
        # NOTE: the sets incoming and outgoing are NOT logically linked!!!
        # This means that it is the PROGRAMMER'S responsibility to ensure that
        # each arc in the MDD is represented twice, in both head.incoming and
        # tail.outgoing!
        if incoming is None:
            self.incoming = set()
        if outgoing is None:
            self.outgoing = set()
        # _tmp is an internal attribute that is used as temporary storage for
        # various MDD calculations (e.g., shortest/longest path)
        self._tmp = None

    def __str__(self):
        return '<in=' + str(self.incoming) + ', out=' + str(self.outgoing) + '>'

    def __repr__(self):
        return 'MDDNodeInfo(' + repr(self.incoming) + ', ' + repr(self.outgoing) + ')'

class MDD(object):
    """MDD represents a multivalued decision diagram (MDD).

    MDD represents a multivalued decision diagram, or MDD.

    Args:
        name (str): name of MDD (default: 'mdd')
        nodes (List[Dict[MDDNode, MDDNodeInfo]]): nodes of MDD;
            if None (default), set to empty list
    """

    def __init__(self, name='mdd', nodes=None):
        """Construct a new 'MDD' object."""
        # 'nodes' is a list of dicts (one for each node layer),
        # and each dict stores the nodes in that layer;
        # each node is represented as a (MDDNode, MDDNodeInfo) key-value pair
        self.nodes = nodes
        self.name = name
        if self.nodes is None:
            self.nodes = []

    @property
    def numNodeLayers(self):
        """Number of node layers; equal to number of 'variables' + 1."""
        return len(self.nodes)
    @property
    def numArcLayers(self):
        """Number of arc layers; equal to number of 'variables'."""
        return len(self.nodes)-1

    @property
    def widthList(self):
        """Number of nodes in each layer"""
        return list(len(lyr) for lyr in self.nodes)
    @property
    def maxWidth(self):
        """Maximum number of nodes in a single node layer."""
        return max(len(lyr) for lyr in self.nodes)

    def __str__(self, showLong=False, showIncoming=False):
        """Return str(self).
        
        Return a (human-readable) string representation of the MDD.

        Args:
            showLong (bool): use more vertical space (default: False)
            showIncoming (bool): show incoming arcs (default: False)

        Returns:
            str: string representation of MDD
        """
        s = '== MDD (' + self.name + ', ' + str(self.numArcLayers) + ' layers) ==\n'
        if showLong:
            # Long form
            s += '# Nodes\n'
            for j in range(len(self.nodes)):
                s += 'L' + str(j) + ':\n'
                for v in self.nodes[j]:
                    s += '\t' + str(v) + ': <'
                    s += 'in={' + ', '.join(str(a) for a in self.nodes[j][v].incoming) + '}, '
                    s += 'out={' + ', '.join(str(a) for a in self.nodes[j][v].outgoing) + '}'
                    s += '>\n'
            s += '# (Outgoing) Arcs\n'
            s += '\n'.join(str(a) for a in self.alloutgoingarcs())
            if showIncoming:
                s += '\n# (Incoming) Arcs\n'
                s += '\n'.join(str(a) for a in self.allincomingarcs())
        else:
            # Short form
            s += '# Nodes\n'
            for j in range(len(self.nodes)):
                s += 'L' + str(j) + ': '
                s += ', '.join(str(v) for v in self.allnodes_in_layer(j)) + '\n'
            s += '# (Outgoing) Arcs\n'
            s += ', '.join(str(a) for a in self.alloutgoingarcs())
            if showIncoming:
                s += '\n# (Incoming) Arcs\n'
                s += ', '.join(str(a) for a in self.allincomingarcs())
        return s

    def __repr__(self):
        return 'MDD(' + repr(self.name) + ', ' + repr(self.nodes) + ')'

    def _get_node_info(self, node):
        """Get 'MDDNodeInfo' corresponding to 'node'.

        Get the 'MDDNodeInfo' object corresponding to the 'MDDNode'
        object 'node'. Note this function can *not* be used to populate
        the underlying dictionary; it can only be used to reference
        the object.

        NOTE: In general, you should use allnodeitems_in_layer(...) if you
        want to update node info in a systematic manner. The author
        recommends only using this function if allnodeitems_in_layer(...)
        cannot be used.
        """
        return self.nodes[node.layer][node]

    def _add_arc(self, newarc):
        """Add an arc to the MDD, without sanity checks."""
        self._get_node_info(newarc.tail).outgoing.add(newarc)
        self._get_node_info(newarc.head).incoming.add(newarc)

    def _remove_arc(self, rmvarc):
        """Remove an arc from the MDD, without sanity checks."""
        self._get_node_info(rmvarc.tail).outgoing.remove(rmvarc)
        self._get_node_info(rmvarc.head).incoming.remove(rmvarc)

    def _add_node(self, newnode):
        """Add a node to the MDD, without sanity checks."""
        # NOTE: If an identical node already exists, its incoming and outgoing
        # arcs will be ERASED!!!
        self.nodes[newnode.layer][newnode] = MDDNodeInfo()

    def _remove_node(self, rmvnode):
        """Remove a node from the MDD, without sanity checks."""
        for arc in self._get_node_info(rmvnode).incoming:
            self._get_node_info(arc.tail).outgoing.remove(arc)
        for arc in self._get_node_info(rmvnode).outgoing:
            self._get_node_info(arc.head).incoming.remove(arc)
        del self.nodes[rmvnode.layer][rmvnode]

    def _remove_nodes(self, rmvnodes):
        """Remove a list of nodes from the MDD, without sanity checks."""
        for v in rmvnodes:
            self._remove_node(v)

    # Default inarcfun, outarcfun methods
    @staticmethod
    def _default_inarcfun(mgnode, inarc, lyr):
        return MDDArc(inarc.label, inarc.weight, inarc.tail, mgnode)
    @staticmethod
    def _default_outarcfun(mgnode, outarc, lyr):
        return MDDArc(outarc.label, outarc.weight, mgnode, outarc.head)

    #
    def _merge_nodes_internal(self, mnodes, mlayer, nodefun, inarcfun=None, outarcfun=None):
        """Merge specified nodes into a new node, with MDDNode/MDDArc functions.

        Merge specified nodes into new supernode and modify arcs appropriately.
        The difference between this function and _merge_nodes is that nodefun,
        inarcfun, and outarcfun directly return MDDNodes and MDDArcs (as opposed
        to returning the new node state and new arc weights).

        Args:
            mnodes (List[MDDNode]): nodes to be merged together
            mlayer (int): layer containing merged nodes
                NOTE: all nodes in mnodes must be in layer mlayer
            nodefun (Callable[[List[MDDNode], int], MDDNode]):
                nodefun(vlist, j) returns the node resulting from merging nodes
                in 'vlist' in layer 'j'
            inarcfun (Callable[[MDDNode, MDDArc, int], MDDArc]):
                inarcfun(mgnode, inarc, j) returns the arc (corresponding to
                'inarc') incoming to the new merged node 'mgnode' in layer 'j';
                if inarcfun is None (default), the original 'inarc' data is
                used unchanged
                NOTE: head of returned arc must be 'mgnode'
            outarcfun (Callable[[MDDNode, MDDArc, int], MDDArc]):
                outarcfun(mgnode, outarc, j) returns the arc (corresponding to
                'outarc') outgoing from the new merged node 'mgnode'
                in layer 'j';
                if outarcfun is None (default), the original 'outarc' data is
                used unchanged
                NOTE: tail of returned arc must be 'mgnode'

        Returns:
            MDDNode: new merged supernode

        Raises:
            ValueError: cannot merge < 2 nodes
        """
        # Basic check
        if len(mnodes) < 2:
            raise ValueError('Cannot merge < 2 nodes: %s' % str(mnodes))
        # Use default inarcfun/outarcfun if unspecified
        if inarcfun is None:
            inarcfun = self._default_inarcfun
        if outarcfun is None:
            outarcfun = self._default_outarcfun

        # Enumerate incoming/outgoing arcs
        mIncoming = set(chain.from_iterable(self.nodes[mlayer][v].incoming for v in mnodes))
        mOutgoing = set(chain.from_iterable(self.nodes[mlayer][v].outgoing for v in mnodes))

        # Create new supernode, and new incoming/outgoing arcs
        mNode = nodefun(mnodes, mlayer)
        newIncoming = [inarcfun(mNode, arc, mlayer) for arc in mIncoming]
        newOutgoing = [outarcfun(mNode, arc, mlayer) for arc in mOutgoing]
        # Delete merged nodes
        self._remove_nodes(mnodes)
        # Add supernode and its arcs to MDD
        self._add_node(mNode)
        for arc in newIncoming:
            self._add_arc(arc)
        for arc in newOutgoing:
            self._add_arc(arc)

        # Return new merged supernode
        return mNode

    # Default awfun method
    @staticmethod
    def _default_awfun(w,ns,nt,j):
        return w

    def _merge_nodes(self, mnodes, mlayer, nsfun, awinfun=None, awoutfun=None):
        """Merge specified nodes into a new node, with state/weight functions.

        Merge specified nodes into new supernode and modify arcs appropriately.
        The difference between this function and _merge_nodes_internal is that
        nsfun, awinfun, and awoutfun return the new node state and new arc
        weights for the merged supernode (as opposed to directly returning
        MDDNodes and MDDArcs).

        Args:
            mnodes (List[MDDNode]): nodes to be merged together
            mlayer (int): layer containing merged nodes
            nsfun (Callable[[List[object], int], object]): nsfun(slist,j) returns the
                node state resulting from merging node states in 'slist'
                in layer 'j'
            awinfun (Callable[[float, object, object, int], float]):
                awinfun(w,os,ms,j) returns the adjusted weight of an arc with
                weight 'w', old head node state 'os', and new head node (i.e.,
                merged supernode in layer 'j') state 'ms';
                if awinfun is None (default), the original weight is used
            awoutfun (Callable[[float, object, object, int], float]):
                awoutfun(w,os,ms,j) returns the adjusted weight of an arc with
                weight 'w', old tail node state 'os', and new tail node (i.e.,
                merged supernode in layer'j') state 'ms';
                if awoutfun is None (default), the original weight is used

        Returns:
            MDDNode: new merged supernode
        """
        # Use default awfun if unspecified
        if awinfun is None:
            awinfun = self._default_awfun
        if awoutfun is None:
            awoutfun = self._default_awfun

        def nodefun(vlist, lyr):
            return MDDNode(mlayer, nsfun([v.state for v in vlist], lyr))
        def inarcfun(mgnode, inarc, lyr):
            return MDDArc(inarc.label, awinfun(inarc.weight, inarc.head.state, mgnode.state, lyr), inarc.tail, mgnode)
        def outarcfun(mgnode, outarc, lyr):
            return MDDArc(outarc.label, awoutfun(outarc.weight, outarc.tail.state, mgnode.state, lyr), mgnode, outarc.head)
        return self._merge_nodes_internal(mnodes, mlayer, nodefun, inarcfun, outarcfun)

    def _redirect_incoming_arcs(self, old_head, new_head):
        """Redirect all incoming arcs of one node to another node.

        Redirect all incoming arcs of node 'old_head' to another node
        'new_head'.  The 'old_head' node is then pruned from the MDD (in
        particular, its outgoing arcs are removed).

        Args:
            old_head (MDDNode): node whose incoming arcs are redirected
            new_head (MDDNode): new head node of redirected arcs
        """
        # Redirect incoming arcs
        for arc in self._get_node_info(old_head).incoming:
            arc.head = new_head
            self._get_node_info(new_head).incoming.add(arc)
            # Note: because we are operating on arcs in the MDD, this also
            # updates the corresponding arcs in the outgoing arc sets,
            # i.e., we do not need to manually update them!

        # Prune old_head node
        self.nodes[old_head.layer][old_head].incoming = set()
        self.prune_recursive([old_head])

    def _append_new_layer(self):
        """Append a new layer to the MDD."""
        self.nodes.append(dict())
        
    def _clear(self):
        """Reset the MDD."""
        self.nodes = []

    def _reset_tmp(self):
        """Reset tmp attribute."""
        for j in range(self.numNodeLayers):
            for (u, ui) in self.allnodeitems_in_layer(j):
                ui._tmp = None

    def allnodes(self):
        """Return all MDDNodes in the MDD."""
        return chain.from_iterable(l.keys() for l in self.nodes)

    def allnodeitems_in_layer(self, layer):
        """Return all (MDDNode, MDDNodeInfo) pairs in a particular layer."""
        return self.nodes[layer].items()

    def allnodes_in_layer(self, layer):
        """Return all MDDNodes in a particular layer."""
        return self.nodes[layer].keys()

    def alloutgoingarcs(self):
        """Return all outgoing arcs in the MDD."""
        return chain.from_iterable(ui.outgoing for j in range(self.numArcLayers) for ui in self.nodes[j].values())

    def allincomingarcs(self):
        """Return all incoming arcs in the MDD."""
        return chain.from_iterable(ui.incoming for j in range(self.numArcLayers) for ui in self.nodes[j+1].values())

    def add_arc(self, newarc):
        """Add an arc to the MDD.

        Add an arc to the MDD (with sanity checks).

        Args:
            newarc (MDDArc): arc to be added

        Raises:
            RuntimeError: head/tail node of arc does not exist
            ValueError: head and tail nodes must be one layer apart
        """
        if not newarc.tail in self.allnodes_in_layer(newarc.tail.layer):
            raise RuntimeError('tail node of arc does not exist')
        if not newarc.head in self.allnodes_in_layer(newarc.head.layer):
            raise RuntimeError('head node of arc does not exist')
        if newarc.head.layer != newarc.tail.layer + 1:
            raise ValueError('head and tail must be one layer apart (%d != %d + 1)' % (newarc.head.layer, newarc.tail.layer))
        self._add_arc(newarc)

    def remove_arc(self, rmvarc):
        """Remove an arc from the MDD.

        Remove an arc from the MDD (with sanity checks).

        Args:
            rmvarc (MDDArc): arc to be removed

        Raises:
            RuntimeError: head/tail node of arc does not exist
            KeyError: no such incoming/outgoing arc exists in the MDD
        """
        if not rmvarc.tail in self.allnodes_in_layer(rmvarc.tail.layer):
            raise RuntimeError('tail node of arc does not exist')
        if not rmvarc.head in self.allnodes_in_layer(rmvarc.head.layer):
            raise RuntimeError('head node of arc does not exist')
        if not rmvarc in self._get_node_info(rmvarc.tail).outgoing:
            raise KeyError('cannot remove non-existent outgoing arc')
        if not rmvarc in self._get_node_info(rmvarc.head).incoming:
            raise KeyError('cannot remove non-existent incoming arc')
        self._remove_arc(rmvarc)
    
    def add_node(self, newnode):
        """Add a new node to the MDD.

        Add a new node to the MDD (with sanity checks).
        
        Args:
            newnode (MDDNode): node to be added

        Raises:
            IndexError: the MDD does not contain the specified node layer
            ValueError: a duplicate node already exists in the MDD
        """
        if newnode.layer >= self.numNodeLayers or newnode.layer < 0:
            raise IndexError('node layer %d does not exist' % newnode.layer)
        if newnode in self.allnodes_in_layer(newnode.layer):
            raise ValueError('cannot add proposed node; duplicate node already exists')
        self._add_node(newnode)

    def remove_node(self, rmvnode):
        """Remove a node from the MDD.

        Remove a node from the MDD (with sanity checks).

        Args:
            rmvnode (MDDNode): node to be removed

        Raises:
            IndexError: the MDD does not contain the specified node layer
            KeyError: no such node exists in the MDD
        """
        if rmvnode.layer >= self.numNodeLayers or rmvnode.layer < 0:
            raise IndexError('node layer %d does not exist' % rmvnode.layer)
        if not rmvnode in self.allnodes_in_layer(rmvnode.layer):
            raise KeyError('cannot remove non-existent node')
        self._remove_node(rmvnode)

    def merge_nodes(self, mnodes, nsfun, awinfun=None, awoutfun=None):
        """Merge specified nodes into a new node.

        Merge nodes in mnodes into a new supernode.  The state of the new
        merged node is specified by nsfun, while the weights of incoming
        and/or outgoing arcs are specified by awinfun/awoutfun.

        Args:
            mnodes (List[MDDNode]): nodes to be merged together
            nsfun (Callable[[List[object], int], object]):
                nsfun(slist, j) returns the node state resulting from merging
                node states in 'slist' in layer 'j'
            awinfun (Callable[[float, object, object, int], float]):
                awinfun(w,os,ms,j) returns the adjusted weight of an arc with
                weight 'w', old head node state 'os', and new head node (i.e.,
                merged supernode in layer 'j') state 'ms';
                if awinfun is None (default), the original weight is used
            awoutfun (Callable[[float, object, object, int], float]):
                awoutfun(w,os,ms,j) returns the adjusted weight of an arc with
                weight 'w', old tail node state 'os', and new tail node (i.e.,
                merged supernode in layer 'j') state 'ms';
                if awoutfun is None (default), the original weight is used

        Returns:
            MDDNode: new merged supernode

        Raises:
            ValueError: cannot merge nodes in different layers
        """
        # Check all nodes in mnodes are on same layer
        mlayer = [v.layer for v in mnodes]
        if len(set(mlayer)) > 1:
            raise ValueError('cannot merge nodes in different layers')
        return self._merge_nodes(mnodes, mlayer[0], nsfun, awinfun, awoutfun)

    def prune_all(self):
        """Prune all dead nodes with no root/terminal path.

        Prune all nodes from the MDD that do not have a path to the first or
        last layer.  This procedure employs two passes to prune all dead nodes
        in the MDD.  To prune more selectively (e.g., after removing a small
        number of nodes/arcs), use prune_recusive(...).
        """
        # Go up the MDD, from second to last layer to top
        # Note last node layer is numNodeLayers-1
        for j in range(self.numNodeLayers-2, -1, -1):
            # Delete nodes without any outgoing arcs
            # bool(set()) = False, so if x is a set, 'x is empty' == not x
            prnnodes = [u for (u, ui) in self.allnodeitems_in_layer(j) if not ui.outgoing]
            for u in prnnodes:
                self._remove_node(u)

        # Go down the MDD, from second layer to bottom
        for j in range(1, self.numNodeLayers):
            # Delete nodes without any incoming arcs
            prnnodes = [u for (u, ui) in self.allnodeitems_in_layer(j) if not ui.incoming]
            for u in prnnodes:
                self._remove_node(u)

    def prune_recursive(self, origNodeList):
        """Recursively prune dead nodes with no root/terminal path.

        Recursively prune nodes from the MDD that do not have a path to the
        first or last layer.  Specifically, this procedure first initializes a
        FIFO queue of nodes to be pruned from 'origNodeList'.  Then, it pops a
        node 'u' from the queue, removes it from the MDD, then adds any nodes
        'v' neighboring 'u' to the queue if the removal of 'u' results in 'v'
        having no incoming or outgoing arcs.  The process is repeated until
        the queue is empty.  Note that if one wishes to prune all dead nodes,
        prune_all() may be more efficient.

        Args:
            origNodeList (List[MDDNode]): nodes to be potentially pruned
        """
        # Recursively delete nodes that cannot reach the first or last layer
        def prunable(u):
            if not self.nodes[u.layer][u].incoming and u.layer > 0:
                return True
            elif not self.nodes[u.layer][u].outgoing and u.layer < self.numNodeLayers-1:
                return True
            else:
                return False
        prnnodes = deque(u for u in origNodeList if prunable(u))
        while len(prnnodes) > 0:
            u = prnnodes.popleft()
            for arc in self._get_node_info(u).incoming:
                if len(self._get_node_info(arc.tail).outgoing) <= 1:
                    prnnodes.append(arc.tail)
            for arc in self._get_node_info(u).outgoing:
                if len(self._get_node_info(arc.head).incoming) <= 1:
                    prnnodes.append(arc.head)
            self._remove_node(u)

    def compile_top_down(self, numLayers, domainFunc, trFunc, costFunc, rootState, isFeas, maxWidth=None, mergeFunc=None, adjFunc=None, nodeSelFunc=None):
        """Compile the MDD top-down according to a DP formulation.

        Perform a top-down compilation of the MDD according to a dynamic
        programming (DP) formulation.  For the DP-specifying functions
        domainFunc, trFunc, and costFunc, layers should be numbered
        0, 1, ..., numLayers-1.

        Note that 'mergeFunc', 'adjFunc', and 'nodeSelFunc' are only used
        when compiling restricted/relaxed MDDs, so if 'maxWidth' is None,
        (i.e., compiling an exact MDD), they do not need to be set.

        Args:
            numLayers (int): number of (arc) layers (i.e., variables)
            domainFunc (Callable[[int], List[int]]): domainFunc(j) returns the
                domain of layer 'j'
            trFunc (Callable[[object, int, int], object]): trFunc(s,d,j)
                returns the state transitioned to when domain value 'd' is
                selected at current state 's', current layer 'j'
            costFunc (Callable[[object, int, int, object], float]):
                costFunc(s,d,j,ns) returns the cost of selecting domain value
                'd' at current state 's', current layer 'j', resulting in next
                state 'ns' (i.e., the output of trFunc(s,d,j))
            rootState (object): state of the root node
            isFeas (Callable[[object, int], bool]): isFeas(s,j) returns True if
                node state 's' in layer 'j' is feasible and False otherwise
            maxWidth (Callable[[int], int]): maxWidth(j) returns the maximum
                allowable width for layer 'j' of the MDD; if None (default),
                the maximum width is set to +inf and compile_top_down(...)
                returns an exact MDD
            mergeFunc (Callable[[List[object], int], object]):
                mergeFunc(slist,j) returns the node state resulting from
                merging node states in 'slist' in layer 'j' (default: None)
            adjFunc (Callable[[float, object, object, int], float]):
                adjFunc(w,os,ms,j) returns the adjusted weight of an arc with
                weight 'w', old head node state 'os', and new head node (i.e.,
                merged supernode in layer 'j') state 'ms' (default: None)
            nodeSelFunc (Callable[[List[MDDNode], int], List[MDDNode]]):
                nodeSelFunc(vlist,j) returns a list of nodes selected from
                'vlist' in layer 'j' to be either merged (if mergeFunc and
                adjFunc are defined) or removed (if mergeFunc and adjFunc
                are None); if nodeSelFunc is None (default), it picks either
                two (if merging) or one (if removing) node(s) at random

        Raises:
            RuntimeError: mergeFunc and adjFunc must be defined together
            RuntimeError: no more nodes to remove/merge but currWidth > maxWidth
        """
        # Basic parameter checks and default settings
        if maxWidth is None:
            maxWidth = lambda j: float('inf')
        else:
            if mergeFunc is None != adjFunc is None:
                raise RuntimeError('mergeFunc and adjFunc must be defined together')
            if nodeSelFunc is None:
                nodeSelFunc = lambda vlist,j: sample(vlist, 1 + int(mergeFunc is not None))
        # First, clear the MDD
        self._clear()
        # Create first layer, containing only the root
        self._append_new_layer()
        self._add_node(MDDNode(0, rootState))
        for j in range(numLayers):
            # Merge/Remove until current layer is under maxWidth
            currLayer = [u for u in self.allnodes_in_layer(j)]
            currWidth = len(currLayer)
            while currWidth > maxWidth(j):
                mnodes = nodeSelFunc(currLayer,j)
                if mergeFunc is None:   # Remove
                    self._remove_nodes(mnodes)
                else:                   # Merge
                    self._merge_nodes(mnodes, j, mergeFunc, adjFunc)
                currLayer = [u for u in self.allnodes_in_layer(j)]
                if len(currLayer) >= currWidth:
                    raise RuntimeError('no more nodes to remove/merge but width of layer %d > %d' % (j,maxWidth(j)))
                currWidth = len(currLayer)

            # Create the next layer of nodes
            self._append_new_layer()
            # For each node in the current layer and each possible assignment...
            for u in self.allnodes_in_layer(j):
                for d in domainFunc(j):
                    # Apply transition function
                    vstate = trFunc(u.state, d, j)
                    # Add if node is feasible
                    if isFeas(vstate, j+1):
                        v = MDDNode(j+1, vstate)
                        # Check if equivalent node exists
                        if v not in self.allnodes_in_layer(j+1):
                            self._add_node(v)
                        # Add appropriate arc
                        self._add_arc(MDDArc(d, costFunc(u.state, d, j, vstate), u, v))

    def compile_trivial(self, numLayers, domainFunc, costFunc, nodeStateFunc):
        """Compile a trivial MDD.

        Compile a trivial MDD, i.e., one that represents all possible solutions
        based on the given domains and costs.  This MDD contains one node in
        each layer, and an arc for every possible domain value.  For domainFunc
        costFunc, and nodeStateFunc, the layers should be numbered 0, 1, ...,
        numLayers-1.

        Args:
            numLayers (int): number of (arc) layers (i.e., variables)
            domainFunc (Callable[[int], List[int]]): domainFunc(j) returns the
                domain of layer 'j'
            costFunc (Callable[[int, int], float]): costFunc(d,j) returns the
                cost of selecting domain value 'd' at current layer 'j'
            nodeStateFunc (Callable[[int], object]): nodeStateFunc(j) returns
                the node state of the node in layer 'j'
        """
        # First, clear the MDD
        self._clear()
        # Create root node
        self._append_new_layer()
        self._add_node(MDDNode(0, nodeStateFunc(0)))
        # Create each layer
        for j in range(numLayers):
            self._append_new_layer()
            for u in self.allnodes_in_layer(j):
                v = MDDNode(j+1, nodeStateFunc(j+1))
                self._add_node(v)
                for d in domainFunc(j):
                    self._add_arc(MDDArc(d, costFunc(d, j), u, v))

    def filter_and_refine_constraint(self, trFunc, rootState, isFeas, nodeStateFunc, maxWidth=None):
        """Filter and refine MDD for a constraint.

        Perform incremental refinement of a particular constraint on the MDD.

        Args:
            trFunc (Callable[[object, int, int], object]): trFunc(s,d,j) returns
                the state transitioned to when domain value 'd' is selected at
                current state 's', current layer 'j'
            rootState (object): state assigned to root node
            isFeas (Callable[[object, int], bool]): isFeas(s,j) returns True if
                node state 's' in layer 'j' is feasible and False otherwise
            nodeStateFunc (Callable[[object, int], object]): nodeStateFunc(s,j)
                returns the node state resulting from transitioning to state
                's' in layer 'j'
            maxWidth (Callable[[int], int]): maxWidth(j) returns the maximum
                allowable width for layer 'j' of the MDD; if None (default),
                the maximum width is set to 100 for all layers
        """
        if maxWidth is None:
            maxWidth = lambda j: 100
        # Initialize tmpdict
        tmpdict = dict()
        # Set up root node
        for u in self.allnodes_in_layer(0):
            tmpdict[u] = rootState
        for j in range(self.numArcLayers):
            # Filtering
            for (u, ui) in self.allnodeitems_in_layer(j):
                for a in list(ui.outgoing):
                    if not isFeas(trFunc(tmpdict[u], a.label, j), j):
                        self._remove_arc(a)

            # Prune nodes that are no longer reachable
            for (u, ui) in list(self.allnodeitems_in_layer(j)):
                if not ui.outgoing:
                    self._remove_node(u)
            for (u, ui) in list(self.allnodeitems_in_layer(j+1)):
                if not ui.incoming:
                    self._remove_node(u)

            # Refinement
            numNodes = len(self.allnodes_in_layer(j+1))
            for (u, ui) in self.allnodeitems_in_layer(j):
                for a in list(ui.outgoing):
                    ns = trFunc(tmpdict[u], a.label, j)
                    v = a.head
                    if v not in tmpdict:
                        tmpdict[v] = ns
                    elif tmpdict[v] != ns and numNodes < maxWidth(j+1):
                        # Redirect arc to a new node
                        w = MDDNode(j+1, nodeStateFunc(ns, j+1))
                        self._add_node(w)
                        tmpdict[w] = ns
                        self._add_arc(MDDArc(a.label, a.weight, u, w))
                        self._remove_arc(a)
                        # Copy outgoing arcs from v, to w
                        for aa in self._get_node_info(v).outgoing:
                            self._add_arc(MDDArc(aa.label, aa.weight, w, aa.head))
                        numNodes += 1
                    else:
                        # Update
                        tmpdict[v] = ns

        # Reset tmp attribute
        tmpdict.clear()

    def compile_pathlist(self, pathList, minArcCostFunc=None, nodeStateFunc=None):
        """Compile an MDD from a list of paths.

        Compile an MDD based on a list of paths (and their associated weights).
        The weight of each arc is determined according to the canonical arc
        cost construction. If 'minArcCost(k)' is 'L', then in the constructed
        MDD the weight of every arc in layer 'k' will be at least 'L'.

        Args:
            pathList (List[Tuple[float, List[object]]]): list of tuples of
                path weights and paths
            minArcCostFunc (Callable[[int], float]): minArcCostFunc(j) returns
                the desired minimum arc weight in layer 'j' of the constructed
                MDD (layer 0 is not considered); if None (default), returns
                0.0 for every layer
            nodeStateFunc (Callable[[int,int], object]): nodeStateFunc(j,k)
                returns the node state of the 'k'th node in layer 'j'; if
                None (default), returns the tuple '(j,k)'

        Raises:
            RuntimeError: numLayers is not consistent with pathList
        """
        # Sanity check and count number of layers needed
        numLayers = max(len(p[1]) for p in pathList)
        if numLayers != min(len(p[1]) for p in pathList):
            raise RuntimeError('number of layers in pathList is not consistent')
        if nodeStateFunc is None:
            nodeStateFunc = lambda j,k: (j,k)
        if minArcCostFunc is None:
            minArcCostFunc = lambda j: 0.0

        # First, clear the MDD
        self._clear()
        # Create root node and associate it with full pathList
        self._append_new_layer()
        rootNode = MDDNode(0, nodeStateFunc(0,0))
        self._add_node(rootNode)
        tmpdict = {rootNode: pathList}
        # Create initial tree
        for j in range(numLayers):
            self._append_new_layer()
            nodeIndex = 0
            for (u, ui) in self.allnodeitems_in_layer(j):
                currDomain = frozenset(p[1][j] for p in tmpdict[u])
                for d in currDomain:
                    v = MDDNode(j+1, nodeStateFunc(j+1,nodeIndex))
                    self._add_node(v)
                    tmpdict[v] = [p for p in tmpdict[u] if p[1][j] == d]
                    arcWeight = tmpdict[v][0][0] if j+1 == numLayers else 0.0
                    self._add_arc(MDDArc(d, arcWeight, u, v))
                    nodeIndex += 1
        tmpdict.clear()

        # Define canonical arc costs, bottom up
        for j in range(numLayers-1, 0, -1):
            for (u, ui) in self.allnodeitems_in_layer(j):
                minArcWeight = min(a.weight for a in ui.outgoing)
                for a in ui.outgoing:
                    a.weight = a.weight - (minArcWeight - minArcCostFunc(j))
                for a in ui.incoming:
                    a.weight = a.weight + minArcWeight + minArcCostFunc(j)

    def reduce_bottom_up(self, mergeFunc, adjInFunc=None, adjOutFunc=None, ignoreLastLayer=False):
        """Reduce the MDD bottom-up by merging equivalent nodes.
        
        Merge all equivalent nodes in the MDD, i.e., nodes which have the same
        suffix set.  The state of the new node is determined by mergeFunc, and
        the weight of incoming and outgoing arcs of the new node is adjusted by
        adjInFunc and adjOutFunc respectively.

        Args:
            mergeFunc (Callable[[List[object], int], object]):
                mergeFunc(slist,j) returns the node state resulting from merging
                node states in 'slist' in layer 'j'
            adjInFunc (Callable[[float, object, object, int], float]):
                adjInFunc(w,os,ms,j) returns the adjusted weight of an arc with
                weight 'w', old head node state 'os', and new head node (i.e.,
                merged supernode in layer 'j') state 'ms' (default: None)
            adjOutFunc (Callable[[float, object, object, int], float]):
                adjOutFunc(w,os,ms,j) returns the adjusted weight of an arc with
                weight 'w', old tail node state 'os', and new tail node (i.e.,
                merged supernode in layer 'j') state 'ms' (default: None)
            ignoreLastLayer (bool): whether to avoid merging last layer (i.e.,
                merge all terminal nodes into one) or not (default: False)
        """
        # Merge from bottom up
        for j in range(self.numArcLayers - int(ignoreLastLayer), 0, -1):
            # Cluster nodes by their outNeighbors
            outDict = dict()
            for v in self.allnodes_in_layer(j):
                outNeighbors = frozenset([(a.label, a.head, a.weight) for a in self.nodes[j][v].outgoing])
                if outNeighbors not in outDict.keys():
                    outDict[outNeighbors] = []
                outDict[outNeighbors].append(v)

            # Nodes that have the same outNeighbors can be merged together
            for mnodes in outDict.values():
                if len(mnodes) >= 2:
                    self._merge_nodes(mnodes, j, mergeFunc, adjInFunc, adjOutFunc)

    def _find_optimal_path(self, srcNds, destNds, longest):
        """Find an 'optimal' src-dest path.

        Find an 'optimal' path between a set of source nodes and a set of
        destination nodes using bidirectional search.  Note 'srcNds'/'destNds'
        must only contain nodes from the same layer.  Note that if srcNds is
        on a later layer than destNodes, then the path returned is traversed
        in *reverse*, i.e., from bottom to top.

        Also, the MDD cannot contain any long arcs.

        Args:
            srcNds (List[MDDNode]): list of source nodes
            destNds (List[MDDNode]): list of destination nodes
            longest (bool): True/False if computing longest/shortest path resp

        Returns:
            Tuple[float, List[object]]: optimal weight and optimal path
        """
        # Check that all nodes in srcNds/destNds are on the same layer
        if len(set(v.layer for v in srcNds)) > 1:
            raise ValueError('srcNds cannot contain nodes from different layers')
        if len(set(v.layer for v in destNds)) > 1:
            raise ValueError('destNds cannot contain nodes from different layers')
        # Determine iteration parameters
        if longest:
            (limVal, limCmp) = (float('-inf'), lambda x,y: x < y)
        else:
            (limVal, limCmp) = (float('inf'), lambda x,y: x > y)
        srcLayer = srcNds[0].layer
        destLayer = destNds[0].layer
        if srcLayer == destLayer:
            if len(set(srcNds) & set(destNds)) > 0:
                return (0, [])
            else:
                return (limVal, [])
        if srcLayer < destLayer:
            # Computing suffixes
            iterDir = 1
            (nextArcs, otherEnd, oppEnd) = ('outgoing', 'head', 'tail')
        else:
            # Computing prefixes
            iterDir = -1
            (nextArcs, otherEnd, oppEnd) = ('incoming', 'tail', 'head')
        # Initialize tmpdict
        tmpdict = dict()
        for src in srcNds:
            tmpdict[src] = (0, None)
        toBeProcessed = set(srcNds)
        # Compute the optimal path, layer by layer
        for j in range(srcLayer, destLayer, iterDir):
            nextToProcess = set()
            for u in toBeProcessed:
                ui = self._get_node_info(u)
                for arc in getattr(ui, nextArcs):
                    aot = getattr(arc, otherEnd)
                    if aot not in tmpdict or limCmp(tmpdict[aot][0], tmpdict[u][0] + arc.weight):
                        tmpdict[aot] = (tmpdict[u][0] + arc.weight, arc)
                        nextToProcess.add(aot)
            toBeProcessed = nextToProcess
        # Identify the optimal path, layer by layer in reverse
        optVal = limVal
        optNode = None
        for u in destNds:
            if limCmp(optVal, tmpdict[u][0]):
                optVal = tmpdict[u][0]
                optNode = u
        lpath = []
        for j in range(destLayer, srcLayer, -iterDir):
            optArc = tmpdict[optNode][1]
            lpath.append(optArc.label)
            optNode = getattr(optArc, oppEnd)

        # Clear tmpdict
        tmpdict.clear()

        return (optVal, list(reversed(lpath)))

    def _find_optimal_ix(self, node, suffixes, longest):
        """Find an 'optimal' prefix/suffix from a node.

        Find an 'optimal' prefix/suffix from a node.

        Also, the MDD cannot contain any long arcs.

        Args:
            node (MDDNode): source node
            suffixes (bool): True/False if optimizing over suffixes/prefixes
                of 'src' node resp
            longest (bool): True/False if computing longest/shortest path resp

        Returns:
            Tuple[float, List[object]]: optimal weight and optimal path
        """
        if longest:
            (limVal, limCmp) = (float('-inf'), lambda x,y: x < y)
        else:
            (limVal, limCmp) = (float('inf'), lambda x,y: x > y)
        if suffixes:
            # Computing suffixes
            (lastNodeLayer, iterDir) = (self.numNodeLayers-1, 1)
            (nextArcs, otherEnd, oppEnd) = ('outgoing', 'head', 'tail')
        else:
            # Computing prefixes
            (lastNodeLayer, iterDir) = (0, -1)
            (nextArcs, otherEnd, oppEnd) = ('incoming', 'tail', 'head')
        # Corner case
        if node.layer == lastNodeLayer:
            return (0, [])
        # Initialize tmpdict with source node
        tmpdict = {node: (0, None)}
        toBeProcessed = set([node])
        # Compute the optimal ix, layer by layer
        for j in range(node.layer, lastNodeLayer, iterDir):
            nextToProcess = set()
            for u in toBeProcessed:
                ui = self._get_node_info(u)
                for arc in getattr(ui, nextArcs):
                    aot = getattr(arc, otherEnd)
                    if aot not in tmpdict or limCmp(tmpdict[aot][0], tmpdict[u][0] + arc.weight):
                        tmpdict[aot] = (tmpdict[u][0] + arc.weight, arc)
                        nextToProcess.add(aot)
            toBeProcessed = nextToProcess
        # Identify the optimal ix, layer by layer in 'reverse'
        optVal = limVal
        optNode = None
        for u in self.allnodes_in_layer(lastNodeLayer):
            if limCmp(optVal, tmpdict[u][0]):
                optVal = tmpdict[u][0]
                optNode = u
        lpath = []
        for j in range(lastNodeLayer, node.layer, -iterDir):
            optArc = tmpdict[optNode][1]
            lpath.append(optArc.label)
            optNode = getattr(optArc, oppEnd)
        # Ensure path is always traversed top to bottom
        if suffixes:
            lpath = list(reversed(lpath))

        # Clear tmpdict
        tmpdict.clear()
        return (optVal, lpath)

    def _find_opt_rtpath(self, longest):
        """Find an 'optimal' root-terminal path in the MDD.

        Args:
            longest (bool): True/False if computing longest/shortest path resp

        Returns:
            Tuple[float, List[object]]: optimal weight and optimal path
        """
        if longest:
            (limVal, limCmp) = (float('-inf'), lambda x,y: x < y)
        else:
            (limVal, limCmp) = (float('inf'), lambda x,y: x > y)
        # Initialize tmpdict
        tmpdict = dict()
        for j in range(self.numNodeLayers):
            for u in self.allnodes_in_layer(j):
                if j == 0:
                    tmpdict[u] = (0, None)
                else:
                    tmpdict[u] = (limVal, None)
        # Compute the optimal path, layer by layer
        for j in range(self.numArcLayers):
            for (u, ui) in self.allnodeitems_in_layer(j):
                for a in ui.outgoing:
                    if limCmp(tmpdict[a.head][0], tmpdict[u][0] + a.weight):
                        tmpdict[a.head] = (tmpdict[u][0] + a.weight, a)
        # Identify the optimal path, layer by layer
        optVal = limVal
        optNode = None
        for u in self.allnodes_in_layer(self.numArcLayers):
            if limCmp(optVal, tmpdict[u][0]):
                optVal = tmpdict[u][0]
                optNode = u
        lpath = []
        for j in range(self.numArcLayers, 0, -1):
            optArc = tmpdict[optNode][1]
            lpath.append(optArc.label)
            optNode = optArc.tail

        tmpdict.clear()
        return (optVal, list(reversed(lpath)))

    def find_longest_path(self):
        """Find a longest root-terminal path in the MDD.

        Returns:
            Tuple[float, List[object]]: maximum weight and longest path
        """
        return self._find_opt_rtpath(True)

    def find_shortest_path(self):
        """Find a shortest root-terminal path in the MDD.

        Returns:
            Tuple[float, List[object]]: minimum weight and shortest path
        """
        return self._find_opt_rtpath(False)

    def enumerate_all_paths(self):
        """Enumerate all root-terminal paths in the MDD.

        Returns:
            List[Tuple[float, List[object]]]: list of path weights/paths
        """
        # Initialize tmpdict
        tmpdict = dict()
        for j in range(self.numNodeLayers):
            for u in self.allnodes_in_layer(j):
                tmpdict[u] = []
        # Set up root node
        for (u, ui) in self.allnodeitems_in_layer(0):
            for a in ui.outgoing:
                tmpdict[a.head].append((a.weight, [a.label]))
        # Compute paths, layer by layer.
        for j in range(1, self.numArcLayers):
            for (u, ui) in self.allnodeitems_in_layer(j):
                for a in ui.outgoing:
                    for x in tmpdict[u]:
                        tmpdict[a.head].append((x[0] + a.weight, x[1] + [a.label]))
        # Enumerate paths
        paths = []
        for u in self.allnodes_in_layer(self.numArcLayers):
            paths.extend(tmpdict[u])

        tmpdict.clear()
        return paths

    def _enumerate_fromnode(self, node, suffixes):
        """Enumerate all prefixes/suffixes from a node.

        Args:
            node (MDDNode): node to examine
            suffixes (bool): if True, enumerate all suffixes of 'node';
                otherwise, enumerate all prefixes of 'node'

        Returns:
            List[Tuple[float, List[object]]]: list of -ix weights/-ixes
        """
        if suffixes:
            # Computing suffixes
            lastNodeLayer = self.numNodeLayers-1
            iterRange = range(node.layer+1, self.numNodeLayers)
            nextArcs = 'outgoing'
            otherEnd = 'head'
        else:
            # Computing prefixes
            lastNodeLayer = 0
            iterRange = range(node.layer-1, -1, -1)
            nextArcs = 'incoming'
            otherEnd = 'tail'
        # Corner case
        if node.layer == lastNodeLayer:
            return []

        # Initialize tmpdict
        tmpdict = dict()
        for j in iterRange:
            for u in self.allnodes_in_layer(j):
                tmpdict[u] = []
        # Set up first arc of path
        for a in getattr(self._get_node_info(node), nextArcs):
            tmpdict[getattr(a, otherEnd)].append((a.weight, [a.label]))
        # Compute paths, layer by layer.
        # (NOTE: Technically the last iteration isn't needed, but it makes the code cleaner.)
        for j in iterRange:
            for (u, ui) in self.allnodeitems_in_layer(j):
                for a in getattr(ui, nextArcs):
                    for x in tmpdict[u]:
                        newWeight = x[0] + a.weight
                        newPath = x[1] + [a.label] if suffixes else [a.label] + x[1]
                        tmpdict[getattr(a, otherEnd)].append((newWeight, newPath))
        # Enumerate paths
        ixes = []
        for u in self.allnodes_in_layer(lastNodeLayer):
            ixes.extend(tmpdict[u])
        # Clear tmpdict
        tmpdict.clear()

        return ixes

    def enumerate_all_prefixes(self, node):
        """Enumerate all prefixes of a node.

        Args:
            node (MDDNode): node to examine

        Returns:
            List[Tuple[float, List[object]]]: list of prefix weights/prefixes
        """
        return self._enumerate_fromnode(node, False)

    def enumerate_all_suffixes(self, node):
        """Enumerate all suffixes of a node.

        Args:
            node (MDDNode): node to examine

        Returns:
            List[Tuple[float, List[object]]]: list of suffix weights/suffixes
        """
        return self._enumerate_fromnode(node, True)

    # Default functions/args for GraphViz output
    @staticmethod
    def _default_ndf(state, lyr):
        return 'label="{}"'.format(state)

    @staticmethod
    def _default_adf(label, weight, lyr):
        if label == 0:
            return 'style=dotted,label="{}"'.format(weight)
        else:
            return 'label="{}"'.format(weight)

    _default_asa =  {'key': lambda a: a.label}
    _default_nsa = {'key': lambda v: v.state, 'reverse': True}

    def output_to_dot(self, nodeDotFunc=None, arcDotFunc=None, arcSortArgs=None, nodeSortArgs=None, reverseDir=False):
        """Write the graphical structure of the MDD to a file.

        Write the graphical structure of the MDD to a file (<MDDName>.gv) in
        the DOT language.  The MDD can then be visualized with GraphViz.

        Args:
            nodeDotFunc (Callable[[object, int], str]): nodeDotFunc(s,j)
                returns a string with the DOT options to use given node state
                's' in layer 'j'; if None (default), a sensible default is used
            arcDotFunc (Callable[[object, float, int], str]): arcDotFunc(l,w,j)
                returns a string with the DOT options to use given arc label
                'l', arc weight 'w', and tail node layer 'j'; if None (default),
                a sensible default is used
            arcSortArgs (dict): arguments specifying how to sort a list of arcs
                via list.sort() (i.e., 'key' and, optionally, 'reverse');
                GraphViz then attempts to order the arcs accordingly in the
                output graph; if arcSortArgs is None (default), no such order
                is enforced
            nodeSortArgs (dict): arguments specifying how to sort a list of
                nodes via list.sort() (i.e., 'key' and, optionally, 'reverse');
                GraphViz then attempts to order the nodes accordingly in the
                output graph; if nodeSortArgs is None (default), no such order
                is enforced
            reverseDir (bool): if True, show the MDD with arcs oriented in the
                opposite direction, so the terminal node appears at the top and
                the root node at the bottom (default: False)
        """

        # Use default output functions if unspecified
        if nodeDotFunc is None:
            nodeDotFunc = self._default_ndf
        if arcDotFunc is None:
            arcDotFunc = self._default_adf
        if reverseDir:
            iterRange = range(self.numArcLayers, 0, -1)
            (nextArcAttr, srcAttr, destAttr) = ('incoming', 'head', 'tail')
        else:
            iterRange = range(self.numArcLayers)
            (nextArcAttr, srcAttr, destAttr) = ('outgoing', 'tail', 'head')

        outf = open('{}.gv'.format(self.name), 'w')
        outf.write('digraph {} {{\n'.format(self.name))
        if reverseDir:
            outf.write('edge [dir=back];\n')
        if arcSortArgs is not None:
            outf.write('ordering=out;\n')
        for v in self.allnodes():
            outf.write('{}[{}];\n'.format(hash(v), nodeDotFunc(v.state, v.layer)))
        for j in iterRange:
            for (u, ui) in self.allnodeitems_in_layer(j):
                arcsinlayer = [a for a in getattr(ui, nextArcAttr)]
                if arcSortArgs is not None:
                    arcsinlayer.sort(**arcSortArgs)
                for arc in arcsinlayer:
                    outf.write('{} -> {}[{}];\n'.format(hash(getattr(arc, srcAttr)), hash(getattr(arc, destAttr)), arcDotFunc(arc.label, arc.weight, arc.tail.layer)))
        if nodeSortArgs is not None:
            for j in range(self.numNodeLayers):
                nodesinlayer = [v for v in self.allnodes_in_layer(j)]
                if len(nodesinlayer) > 1:
                    nodesinlayer.sort(**nodeSortArgs)
                    for i in range(len(nodesinlayer) - 1):
                        outf.write('{} -> {}[style=invis];\n'.format(hash(nodesinlayer[i]), hash(nodesinlayer[i+1])))
                    outf.write('{rank=same')
                    for v in nodesinlayer:
                        outf.write(';{}'.format(hash(v)))
                    outf.write('}\n')
        outf.write('}')
        outf.close()

    def dumpJSON(self, stateDumpFunc=repr, labelDumpFunc=repr):
        """Dump the MDD into a JSON file.

        Dump the contents of the MDD into a JSON file for later retrieval.

        Args:
            stateDumpFunc (Callable[[object], str]): stateDumpFunc(s) returns
                a string representation of the node state 's' (default: repr)
            labelDumpFunc (Callable[[object], str]): labelDumpFunc(l) returns
                a string representation of the arc label 'l' (default: repr)
        """
        dataList = []
        dataList.append({'Type': 'name', 'name': self.name})
        for v in self.allnodes():
            dataList.append({'Type': 'node', 'layer': v.layer, 'state': stateDumpFunc(v.state), 'id': hash(v)})
        for a in self.alloutgoingarcs():
            dataList.append({'Type': 'arc', 'label': labelDumpFunc(a.label), 'weight': float(a.weight), 'tail': hash(a.tail), 'head': hash(a.head)})
        outf = open(self.name + '.json', 'w')
        dump(dataList, outf)
        outf.close()

    def loadJSON(self, fname, stateLoadFunc=eval, labelLoadFunc=eval):
        """Load an MDD from a JSON file.

        Load the contents of an MDD from a JSON file.
        NOTE: Since node states and arc labels can be arbitrary python
        objects, loadJSON uses eval() by default to construct these attributes.
        THIS ALLOWS FOR ARBITRARY CODE EXECUTION! USE RESPONSIBLY!!!

        Args:
            fname (str): name of input file (e.g., mdd.json)
            stateLoadFunc (Callable[[str], object]): stateLoadFunc(s) returns
                the node state corresponding to string 's' (default: eval)
            labelLoadFunc (Callable[[str], object]): labelLoadFunc(s) returns
                the arc label corresponding to string 's' (default: eval)

        Raises:
            ValueError: unknown item type (e.g., incorrect input file format)
        """
        self._clear()
        mddf = open(fname, 'r')
        dataList = load(mddf)
        mddf.close()
        nodeDict = dict()
        for item in dataList:
            if item['Type'] == 'name':
                self.name = item['name']
            elif item['Type'] == 'node':
                while int(item['layer']) >= self.numNodeLayers:
                    self._append_new_layer()
                newnode = MDDNode(int(item['layer']), stateLoadFunc(item['state']))
                self.add_node(newnode)
                nodeDict[item['id']] = newnode
            elif item['Type'] == 'arc':
                newarc = MDDArc(labelLoadFunc(item['label']), float(item['weight']), nodeDict[item['tail']], nodeDict[item['head']])
                self.add_arc(newarc)
            else:
                raise ValueError('Unknown item type: check input file format')
