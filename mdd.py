from itertools import chain # used in various places

class MDDArc(object):
    """MDDArc represents a single arc in the MDD.

    MDDArc represents a single arc in the MDD.  An MDDArc is uniquely
    identified by its head/tail nodes, label, and weight.
    """

    def __init__(self, label, weight, tail, head):
        """Construct a new 'MDDArc' object.

        Args:
            label (object): label of arc (e.g., assigned value)
            weight (float): weight of arc (e.g., coefficient)
            tail (MDDNode): tail/source node
            head (MDDNode): head/destination node
        """
        self.label = label
        self.weight = weight
        self.tail = tail
        self.head = head

    # Allows MDDArcs to be used as dictionary keys.
    def __hash__(self):
        """Return the hash value.

        Return the hash value of the 'MDDArc' object.

        Returns:
            int: hash value
        """
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
        return 'A(' + str(self.tail) + ',' + str(self.head) + ':' + str(self.label) + ')'

    def __repr__(self):
        s = 'MDDArc('
        s += repr(self.label) + ', '
        s += repr(self.weight) + ', '
        s += repr(self.tail) + ', '
        s += repr(self.head) + ')'
        return s

class MDDNode(object):
    """MDDNode represents a single node in the MDD.

    MDDNode represents a single in the MDD.  An MDDNode is uniquely identified
    by its layer and state.  The node state must be a hashable object.
    """

    def __init__(self, layer, state):
        """Construct a new 'MDDNode' object.

        Args:
            layer (int): layer the node is in
            state (object): state associated with node
        """
        self.layer = layer
        self.state = state

    # Allows MDDNodes to be used as dictionary keys.
    def __hash__(self):
        """Return the hash value.

        Return the hash value of the 'MDDArc' object.

        Returns:
            int: hash value
        """
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
    """

    def __init__(self, incoming=None, outgoing=None):
        """Construct a new 'MDDNode' object.

        Args:
            incoming (set): set of incoming arcs
            outgoing (set): set of outgoing arcs
        """
        # NOTE: the sets incoming and outgoing are NOT logically linked!!!
        # This means that it is the PROGRAMMER'S responsibility to ensure that
        # each arc in the MDD is represented twice, in both head.incoming and
        # tail.outgoing!
        if incoming is None:
            self.incoming = set()
        if outgoing is None:
            self.outgoing = set()
        # _tmp is an internal attribute that is used as a placeholder during
        # various MDD calculations (e.g., shortest/longest path)
        self._tmp = None

    def __str__(self):
        return '<in=' + str(self.incoming) + ', out=' + str(self.outgoing) + '>'

    def __repr__(self):
        return 'MDDNodeInfo(' + repr(self.incoming) + ', ' + repr(self.outgoing) + ')'

class MDD(object):
    """MDD represents a multivalued decision diagram (MDD).

    MDD represents a multivalued decision diagram, or MDD.
    """

    def __init__(self, name='mdd'):
        """Construct a new 'MDD' object.

        Args:
            name (str): name of MDD (default: 'mdd')

        """
        self.nodes = []
        # Note numLayers = number of ARC layers
        # Number of NODE layers is numLayers+1,
        # i.e., 0, 1, ..., numLayers!!!
        self.numLayers = -1
        self.name = name

    def __str__(self, showLong=False, showIncoming=False):
        """Return str(self).
        
        Return a (human-readable) string representation of the MDD.

        Args:
            showLong: use more vertical space (default: False)
            showIncoming: show incoming arcs (default: False)

        Returns:
            str: string representation of MDD
        """
        s = '== MDD (' + self.name + ', ' + str(self.numLayers) + ' layers) ==\n'
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
        return repr(self.nodes) + '\nnumLayers = ' + repr(self.numLayers)

    # Add an arc to the MDD, without sanity checks.
    def _add_arc(self, newarc):
        self.nodes[newarc.tail.layer][newarc.tail].outgoing.add(newarc)
        self.nodes[newarc.head.layer][newarc.head].incoming.add(newarc)

    # Remove an arc from the MDD, without sanity checks.
    def _remove_arc(self, rmvarc):
        self.nodes[rmvarc.tail.layer][rmvarc.tail].outgoing.remove(rmvarc)
        self.nodes[rmvarc.head.layer][rmvarc.head].incoming.remove(rmvarc)

    # Add a node to the MDD, without sanity checks.
    # NOTE: If an identical node already exists, its incoming and outgoing
    # arcs will be ERASED!!!
    def _add_node(self, newnode):
        self.nodes[newnode.layer][newnode] = MDDNodeInfo()

    # Remove a node from the MDD, without sanity checks.
    def _remove_node(self, rmvnode):
        rmvlayer = rmvnode.layer
        for arc in self.nodes[rmvlayer][rmvnode].incoming:
            self.nodes[rmvlayer-1][arc.tail].outgoing.remove(arc)
        for arc in self.nodes[rmvlayer][rmvnode].outgoing:
            self.nodes[rmvlayer+1][arc.head].incoming.remove(arc)
        del self.nodes[rmvlayer][rmvnode]

    # Default inarcfun, outarcfun methods
    @staticmethod
    def _default_inarcfun(mgnode, inarc):
        return MDDArc(inarc.label, inarc.weight, inarc.tail, mgnode)
    @staticmethod
    def _default_outarcfun(mgnode, outarc):
        return MDDArc(outarc.label, outarc.weight, mgnode, outarc.head)

    # Merge specified nodes into a new node, without sanity checks.
    # Args:
    #     mnodes (list(MDDNode)): nodes to be merged together
    #     mlayer (int): layer containing merged nodes
    #         NOTE: all nodes in mnodes must be in layer mlayer
    #     nodefun (list(MDDNode) -> MDDNode): nodefun(vlist) returns the node
    #         resulting from merging nodes in vlist
    #     inarcfun ((MDDNode, MDDArc) -> MDDArc): inarcfun(mgnode, inarc)
    #         returns the arc (corresponding to inarc) incoming to the new
    #         merged node mgnode; if inarcfun is None (default), the original
    #         inarc data is used unchanged
    #         NOTE: head of returned arc must be mgnode
    #     outarcfun ((MDDNode, MDDArc) -> MDDArc): outarcfun(mgnode, outarc)
    #         returns the arc (corresponding to outarc) outgoing from the new
    #         merged node mgnode; if outarcfun is None (default), the original
    #         outarc data is used unchanged
    #         NOTE: tail of returned arc must be mgnode
    #
    def _merge_nodes(self, mnodes, mlayer, nodefun, inarcfun=None, outarcfun=None):
        # Use default inarcfun/outarcfun if unspecified
        if inarcfun is None:
            inarcfun = self._default_inarcfun
        if outarcfun is None:
            outarcfun = self._default_outarcfun

        # Enumerate incoming/outgoing arcs
        mIncoming = set(chain.from_iterable(self.nodes[mlayer][v].incoming for v in mnodes))
        mOutgoing = set(chain.from_iterable(self.nodes[mlayer][v].outgoing for v in mnodes))

        # Create the new supernode and its arcs
        mNode = nodefun(mnodes)
        self._add_node(mNode)
        for arc in mIncoming:
            self._add_arc(inarcfun(mNode, arc))
        for arc in mOutgoing:
            self._add_arc(outarcfun(mNode, arc))

        # Delete merged nodes (with different state)
        for v in mnodes:
            if v.state != mNode.state:
                self._remove_node(v)


    def merge_nodes(self, mnodes, mlayer, mergefun, adjfun):
        """Merge specified nodes into a new node.

        Merge all nodes in mnodes (in layer mlayer) into a single new node
        with the appropriate arcs.  The state of the new node is determined
        by mergefun and the weight of arcs coming into the new node is
        adjusted by adjfun.

        Args:
            mnodes (list(MDDNode)): nodes to be merged together
            mlayer (int): layer containing merged nodes
            mergefun (list(object) -> object): mergefun(slist) returns the
                node state resulting from merging node states in 'slist'
            adjfun ((float, object, object) -> float): adjfun(w,os,ms) returns
                the adjusted weight of an arc with weight w, tail node state os,
                and head node (i.e., merged supernode) state ms

        Raises:
            RuntimeError: all nodes in mnodes must be in layer mlayer
        """
        if any(v.layer != mlayer for v in mnodes):
            raise RuntimeError('all nodes in mnodes must be in layer mlayer')
            return
        def nodefun(vlist):
            return MDDNode(mlayer, mergefun([v.state for v in vlist]))
        def inarcfun(mgnode, inarc):
            return MDDArc(inarc.label, adjfun(inarc.weight, inarc.tail.state, mgnode.state), inarc.tail, mgnode)
        self._merge_nodes(mnodes, mlayer, nodefun, inarcfun)

    # Append a new layer to the MDD.
    def _append_new_layer(self):
        self.nodes.append(dict())
        self.numLayers = self.numLayers + 1
        
    # Reset the MDD.
    def _clear(self):
        self.nodes = []
        self.numLayers = -1

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
        return chain.from_iterable(ui.outgoing for j in range(self.numLayers) for ui in self.nodes[j].values())

    def allincomingarcs(self):
        """Return all incoming arcs in the MDD."""
        return chain.from_iterable(ui.incoming for j in range(self.numLayers) for ui in self.nodes[j+1].values())
    
    def add_node(self, newnode):
        """Add a new node to the MDD.

        Add a new node to the MDD.
        
        Args:
            newnode (MDDNode): node to be added

        Raises:
            IndexError: the MDD does not contain the specified node layer
            RuntimeError: a duplicate node already exists in the MDD
        """
        if newnode.layer > self.numLayers or newnode.layer < 0:
            raise IndexError('node layer %d does not exist' % newnode.layer)
            return
        if newnode in self.allnodes_in_layer(newnode.layer):
            raise RuntimeError('cannot add proposed node; duplicate node already exists')
        self._add_node(newnode)

    def prune_dead_nodes(self):
        """Prune nodes from the MDD that cannot reach the last layer.

        Prune nodes from the MDD that do not have a path to the last layer.
        """
        # Go up the MDD, from second to last layer to top
        for j in range(self.numLayers-1, 0, -1):
            # Delete nodes without any outgoing arcs
            prnnodes = [u for (u, ui) in self.allnodeitems_in_layer(j) if not ui.outgoing]
            for u in prnnodes:
                self._remove_node(u)

    def compile_top_down(self, numLayers, domainFunc, trFunc, costFunc, rootState, isFeas):
        """Compile the MDD top-down according to a DP formulation.

        Perform a top-down compilation of the MDD according to a dynamic
        programming (DP) formulation. For the DP-specifying functions
        domainFunc, trFunc, and costFunc, layers should be numbered
        0, 1, ..., numLayers-1.

        Args:
            numLayers (int): number of layers (i.e., variables)
            domainFunc (int -> list of ints): domainFunc(j) returns the domain
                of layer 'j'
            trFunc ((object, int, int) -> object): trFunc(s,d,j) returns the
                state transitioned to when domain value 'd' is selected at
                current state 's', current layer 'j'
            costFunc ((object, int, int) -> float): costFunc(s,d,j) returns the
                cost of selecting domain value 'd' at current state 's',
                current layer 'j'
            rootState (object): state of the root node
            isFeas ((object, int) -> bool): isFeas(s,j) returns True if node
                state 's' in layer j is feasible and False otherwise
        """
        # First, clear the MDD
        self._clear()
        # Create first layer, containing only the root
        self._append_new_layer()
        self._add_node(MDDNode(0, rootState))
        for j in range(numLayers):
            # Create the next layer of nodes
            self._append_new_layer()
            # For each node in the current layer and each possible assignment...
            for u in self.allnodes_in_layer(j):
                for d in domainFunc(j):
                    # Apply transition function
                    vstate = trFunc(u.state, d, j)
                    # Add if node is feasible
                    if isFeas(vstate, j):
                        v = MDDNode(j+1, vstate)
                        # Check if equivalent node exists
                        if v not in self.allnodes_in_layer(j+1):
                            self._add_node(v)
                        # Add appropriate arc
                        self._add_arc(MDDArc(d, costFunc(u.state, d, j), u, v))

    def reduce_bottom_up(self, mergeFunc, adjFunc):
        """Reduce the MDD bottom-up by merging equivalent nodes.
        
        Merge all equivalent nodes in the MDD, i.e., nodes which have the same
        suffix set.  The state of the new node is determined by mergeFunc and
        the weight of arcs coming into the new node is adjusted by adjFunc.

        Args:
            mergeFunc ((list(object), int) -> object): mergeFunc(slist,j) returns
                the node state resulting from merging node states in 'slist'
                in layer j
            adjFunc ((float, object, object, int) -> float): adjFunc(w,os,ms,j)
                returns the adjusted weight of an arc with weight w, tail node
                state os, and head node (i.e., merged supernode in layer j)
                state ms
        """
        # Merge from bottom up
        for j in range(self.numLayers, 0, -1):
            # Cluster nodes by their outNeighbors
            outDict = dict()
            for v in self.allnodes_in_layer(j):
                outNeighbors = frozenset([(a.label, a.head, a.weight) for a in self.nodes[j][v].outgoing])
                if outNeighbors not in outDict.keys():
                    outDict[outNeighbors] = []
                outDict[outNeighbors].append(v)

            # Nodes that have the same outNeighbors can be merged together
            for mnodes in outDict.values():
                self.merge_nodes(mnodes, j, lambda slist: mergeFunc(slist,j), lambda w,os,ms: adjFunc(w,os,ms,j))

    # Find an 'optimal' root-terminal path in the MDD.
    # Args:
    #     longest (bool): True/False if computing longest/shortest path resp
    def _find_opt_path(self, longest):
        if longest:
            (limVal, limCmp) = (float('-inf'), lambda x,y: x < y)
        else:
            (limVal, limCmp) = (float('inf'), lambda x,y: x > y)
        # Remember there are self.numLayers ARC layers, and therefore
        # self.numLayers + 1 NODE layers
        for j in range(self.numLayers+1):
            for (u, ui) in self.allnodeitems_in_layer(j):
                if j == 0:
                    ui._tmp = (0, None)
                else:
                    ui._tmp = (limVal, None)
        # Compute the optimal path, layer by layer
        for j in range(self.numLayers):
            for (u, ui) in self.allnodeitems_in_layer(j):
                for a in ui.outgoing:
                    if limCmp(self.nodes[j+1][a.head]._tmp[0], ui._tmp[0] + a.weight):
                        self.nodes[j+1][a.head]._tmp = (ui._tmp[0] + a.weight, a)
        # Identify the optimal path, layer by layer
        optVal = limVal
        optNode = None
        for (u, ui) in self.allnodeitems_in_layer(self.numLayers):
            if limCmp(optVal, ui._tmp[0]):
                optVal = ui._tmp[0]
                optNode = u
        lpath = []
        for j in range(self.numLayers, 0, -1):
            optArc = self.nodes[j][optNode]._tmp[1]
            lpath.append(optArc.label)
            optNode = optArc.tail

        # Reset tmp attribute
        for j in range(self.numLayers+1):
            for (u, ui) in self.allnodeitems_in_layer(j):
                ui._tmp = None

        return (list(reversed(lpath)), optVal)

    def find_longest_path(self):
        """Find a longest root-terminal path in the MDD."""
        return self._find_opt_path(True)

    def find_shortest_path(self):
        """Find a shortest root-terminal path in the MDD."""
        return self._find_opt_path(False)

    # Default functions/args for GraphViz output
    @staticmethod
    def _default_nodedotfunc(state):
        return '[label="' + str(state) + '"];'

    @staticmethod
    def _default_arcdotfunc(label, weight):
        if label == 0:
            return '[style=dotted,label="' + str(weight) + '"];'
        else:
            return '[label="' + str(weight)  + '"];'

    _default_arcsortargs =  {'key': lambda a: a.label}
    _default_nodesortargs = {'key': lambda v: v.state, 'reverse': True}

    def output_to_dot(self, nodeDotFunc=None, arcDotFunc=None, arcSortArgs=None, nodeSortArgs=None):
        """Write the graphical structure of the MDD to a file.

        Write the graphical structure of the MDD to a file (<MDDName>.gv) in
        the DOT language.  The MDD can then be visualized with GraphViz.

        Args:
            nodeDotFunc (object -> str): nodeDotFunc(s) returns a string with
                the DOT options to use given node state 's'; if None (default),
                a sensible default is used
            arcDotFunc (object, float -> str): arcDotFunc(l,w) returns a string
                with the DOT options to use given arc label 'l' and arc weight
                'w'; if None (default), a sensible default is used
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
        """

        # Use default output functions if unspecified
        if nodeDotFunc is None:
            nodeDotFunc = self._default_nodedotfunc
        if arcDotFunc is None:
            arcDotFunc = self._default_arcdotfunc

        outf = open(self.name + '.gv', 'w')
        outf.write('digraph ' + self.name + ' {\n')
        if arcSortArgs is not None:
            outf.write('ordering=out;\n')
        for v in self.allnodes():
            outf.write(str(hash(v)) + nodeDotFunc(v.state) + '\n')
        for j in range(self.numLayers):
            for (u, ui) in self.allnodeitems_in_layer(j):
                if arcSortArgs is not None:
                    arcsinlayer = [a for a in ui.outgoing]
                    arcsinlayer.sort(**arcSortArgs)
                    for arc in arcsinlayer:
                        outf.write(str(hash(arc.tail)) + ' -> ' + str(hash(arc.head)) + arcDotFunc(arc.label, arc.weight) + '\n')
                else:
                    for arc in ui.outgoing:
                        outf.write(str(hash(arc.tail)) + ' -> ' + str(hash(arc.head)) + arcDotFunc(arc.label, arc.weight) + '\n')
        if nodeSortArgs is not None:
            for j in range(self.numLayers+1):
                nodesinlayer = [v for v in self.allnodes_in_layer(j)]
                if len(nodesinlayer) > 1:
                    nodesinlayer.sort(**nodeSortArgs)
                    for i in range(len(nodesinlayer) - 1):
                        outf.write(str(hash(nodesinlayer[i])) + ' -> ' + str(hash(nodesinlayer[i+1])) + '[style=invis];' + '\n')
                    outf.write('{rank=same')
                    for v in nodesinlayer:
                        outf.write(';' + str(hash(v)))
                    outf.write('}\n')
        outf.write('}')
        outf.close()
