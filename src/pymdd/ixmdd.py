from pymdd import mdd
class IxTuple(object):
    """IxTuple represents info about one type of -ix at one MDDNode.

    IxTuple represents information about one type of prefix or suffix
    for a single MDDNode.

    Args:
        node (MDDNode): associated node
        ixtype (str): type of ix
        weight (float): weight of -ix
        nextarc (MDDArc): next arc of optimal -ix (default: None)
        need_check(bool): True if optimal -ix may have changed (default: True)

    """
    def __init__(self, node, ixtype, weight, nextarc=None, need_check=True):
        """Construct a new IxTuple object."""
        #if not need_check:
        #    print('\tSet {}.{} to ({}, {}, {})'.format(node, ixtype, weight, nextarc, need_check))
        self.node = node
        self.ixtype = ixtype
        self._weight = weight
        self._nextarc = nextarc
        self._need_check = need_check

    @property
    def weight(self):
        return self._weight

    @weight.setter
    def weight(self, value):
        #print('\tSet {}.weight to {}'.format(self.ixtype, value))
        self._weight = value

    @property
    def nextarc(self):
        return self._nextarc

    @nextarc.setter
    def nextarc(self, value):
        #print('\tSet {}.nextarc to {}'.format(self.ixtype, value))
        self._nextarc = value

    @property
    def need_check(self):
        return self._need_check

    @need_check.setter
    def need_check(self, value):
        #if not value:
        #    print('\tConfirm {}.{} is ({}, {}, {})'.format(self.node, self.ixtype, self.weight, self.nextarc, value))
        self._need_check = value

    def __str__(self):
        s = '(' + str(self.weight) + ',' + str(self.nextarc) + ')'
        if self.need_check:
            s += '*'
        return s
    def __repr__(self):
        return 'IxTuple(' + repr(self.node) + ', ' + repr(self.ixtype) + ', ' + repr(self.weight) + ', ' + repr(self.nextarc) + ', ' + repr(self.need_check) + ')'

class IxInfo(object):
    """IxInfo represents -ix information about a single MDDNode.

    IxInfo reprsents prefix and suffix information about a single MDDNode.

    Args:
        node (MDDNode): associated node
        min_suffix (IxTuple): minimum weight suffix
        max_suffix (IxTuple): maximum weight suffix
        min_prefix (IxTuple): minimum weight prefix
        max_prefix (IxTuple): maximum weight prefix
    """

    def __init__(self, node, min_suffix=None, max_suffix=None, min_prefix=None, max_prefix=None):
        self.node = node
        self.min_suffix = min_suffix
        if self.min_suffix is None:
            self.min_suffix = IxTuple(node, 'min_suffix', float('inf'))

        self.max_suffix = max_suffix
        if self.max_suffix is None:
            self.max_suffix = IxTuple(node, 'max_suffix', float('-inf'))

        self.min_prefix = min_prefix
        if self.min_prefix is None:
            self.min_prefix = IxTuple(node, 'min_prefix', float('inf'))

        self.max_prefix = max_prefix
        if self.max_prefix is None:
            self.max_prefix = IxTuple(node, 'max_prefix', float('-inf'))

class IxParam(object):
    """IxParam stores parameters specific to each IxInfo.

    IxParam stores parameters specific to each IxInfo.

    Args:
        ixType (str): type of ixinfo
        numNodeLayers: current number of node layers in MDD
    """
    all_ix_types = ('min_suffix', 'max_suffix', 'min_prefix', 'max_prefix')
    def __init__(self, ixType, numNodeLayers):
        if ixType not in self.all_ix_types:
            raise ValueError('unknown ixType: %s' % ixType)
        isSuffix = (ixType[-6:] == 'suffix')
        isMin = (ixType[:3] == 'min')

        self.type = ixType
        self.otherEnd = 'head' if isSuffix else 'tail'
        self.nextArcs = 'outgoing' if isSuffix else 'incoming'
        self.lastNodeLayer = numNodeLayers-1 if isSuffix else 0
        self.iterDir = -1 if isSuffix else 1
        self.limVal = float('inf') if isMin else float('-inf')
        self.limCmp = (lambda a,b: a < b) if isMin else (lambda a,b: a > b)
        self.recEnd = 'tail' if isSuffix else 'head'
        self.recArcs = 'incoming' if isSuffix else 'outgoing'

class IxMDD(mdd.MDD):
    """IxMDD is an MDD that maintains the optimal -ix at every node.

    IxMDD is an extension of the MDD class that maintains the minimum and
    maximum prefix and suffix at every node.

    Args:
        name (str): name of IxMDD (default: 'ixmdd')
        nodes (List[Dict[MDDNode, MDDNodeInfo]]): nodes of MDD;
            if None (default), set to empty list
    """
    def __init__(self, name='ixmdd', nodes=None):
        """Construct a new 'IXMDD' object."""
        #super().__init__(name, nodes)
        mdd.MDD.__init__(self, name, nodes)
        self.ixinfo = {v: IxInfo(v) for v in self.allnodes()}
        self.update_ixinfo_all()

    # Fundamental getters
    def _opt_ixweight(self, ixType, node):
        """Return weight of optimal ix associated with node.

        Return weight of optimal ix associated with node. Note
        this function does NOT verify if optimal ix needs to be
        checked!
        Args:
            ixType (str): type of ix
            node (MDDNode): node being queried

        Returns:
            float: next arc of optimal ix
        """
        return getattr(self.ixinfo[node], ixType).weight
    def _opt_ixarc(self, ixType, node):
        """Return next arc of optimal ix associated with node.
        
        Return next arc of optimal ix associated with node, or None
        if there is no such arc. Make sure any function that uses this
        takes into account that ixarc may be None!!!

        Args:
            ixType (str): type of ix
            node (MDDNode): node being queried

        Returns:
            MDDArc: next arc of optimal ix (or None)
        """
        return getattr(self.ixinfo[node], ixType).nextarc
    def _opt_ixcheck(self, ixType, node):
        """Return if optimal ix associated with node needs to be checked.

        Return if optimal ix associated with node needs to be checked.
        If True, the optimal ix associated with this node may not be accurate,
        and should be found by MDD traversal instead of lookup.

        Args:
            ixType (str): type of ix
            node (MDDNode): node being queried

        Returns:
            bool: True if optimal ix needs to be checked
        """
        return getattr(self.ixinfo[node], ixType).need_check
    def _opt_ixpath(self, ixType, node):
        """Return optimal ix associate with node.

        NOTE: This function returns an error if the total weight of the
        traversed path does not match the optimal ix weight of node. It also
        does not verify whether node ixes need to be checked.

        Args:
            ixType (str): type of ix
            node (MDDNode): node being queried

        Returns:
            Tuple[float, List[object]]: optimal weight and ix path
        """
        ixp = IxParam(ixType, self.numNodeLayers)
        path = []
        w = 0
        currLayer = node.layer
        nextarc = getattr(self.ixinfo[node], ixp.type).nextarc
        while nextarc is not None:
            nextNode = getattr(nextarc, ixp.otherEnd)
            path.append(nextarc.label)
            w += nextarc.weight
            currLayer = nextNode.layer
            nextarc = self._opt_ixarc(ixp.type, nextNode)
        if currLayer != ixp.lastNodeLayer:
            (w, path) = (ixp.limVal, [])
        if w != self._opt_ixweight(ixp.type, node):
            raise RuntimeError('Internal error: weight of path is not optimal weight!')
        return (w, path)

    def opt_ixweight(self, ixType, node):
        """Return weight of optimal ix associated with node.

        Return weight of optimal ix associated with node.
        Args:
            ixType (str): type of ix
            node (MDDNode): node being queried

        Returns:
            float: next arc of optimal ix
        """
        if ixType not in IxParam.all_ix_types:
            raise ValueError('unknown ixType: %s' % ixType)
        if self._opt_ixcheck(ixType, node):
            raise RuntimeWarning('%s needs check: reverting to explicit search' % str(node))
            return self._find_optix(ixType, node)[0]
        else:
            return self._opt_ixweight(ixType, node)


    def opt_ixpath(self, ixType, node):
        """Return optimal ix associate with node.

        NOTE: This function returns an error if the total weight of the
        traversed path does not match the optimal ix weight of node. It also
        does not verify that the ixes are checked.

        Args:
            ixType (str): type of ix
            node (MDDNode): node being queried

        Returns:
            Tuple[float, List[object]]: optimal weight and ix path
        """
        ixp = IxParam(ixType, self.numNodeLayers)
        if self._opt_ixcheck(ixType, node):
            raise RuntimeWarning('%s needs check: reverting to explicit search' % str(node))
            return self._find_optix(ixType, node)
        path = []
        w = 0
        nextarc = getattr(self.ixinfo[node], ixp.type).nextarc
        while nextarc is not None:
            nextNode = getattr(nextarc, ixp.otherEnd)
            if self._opt_ixcheck(ixp.type, nextNode):
                raise RuntimeWarning('%s needs check: reverting to explicit search' % str(nextNode))
                return self._find_optix(ixType, node)
            path.append(nextarc.label)
            w += nextarc.weight
            nextarc = getattr(self.ixinfo[nextNode], ixp.type).nextarc
        return (w, path)

    def _find_optix(self, ixType, node):
        """Internal function that calls _find_optimal_ix of MDD."""
        suffixes = (ixType[-6:] == 'suffix')
        longest = (ixType[:3] == 'max')
        #return super()._find_optimal_ix(node, suffixes, longest)
        return mdd.MDD._find_optimal_ix(self, node, suffixes, longest)

    def _find_optimal_ix(self, node, suffixes, longest):
        """Replaces _find_optimal_ix of MDD to look at cached value first."""
        ixType = ('max' if longest else 'min') + ('_suffix' if suffixes else '_prefix')
        return self.opt_ixpath(ixType, node)

    ###########################
    # IxInfo Update Functions #
    ###########################

    def _reset_ixinfo(self):
        """Clear IxInfo for all nodes in IxMDD."""
        self.ixinfo = {v: IxInfo(v) for v in self.allnodes()}

    def update_ixinfo_all(self):
        """Update IxInfo for all nodes via bottom-up/top-down passes."""
        # Reset ixinfo
        self._reset_ixinfo()
        # Update each ixinfo
        for ixType in IxParam.all_ix_types:
            self._update_ixinfo_bylayer(ixType)

    def _update_ixinfo_bylayer(self, ixType, termLayer=-1):
        """Update IxInfo layer by layer.

        Update IxInfo, by layer up to termLayer.

        Args:
            ixType (str): type of ix
            termLayer (int): one layer beyond last layer to update ixinfo;
                if termLayer < 0 or termlayer > self.numNodeLayers,
                update all layers (default: -1)
        """
        ixp = IxParam(ixType, self.numNodeLayers)
        terminal_layer = self.numArcLayers-ixp.lastNodeLayer+ixp.iterDir if termLayer < 0 or termLayer >= self.numNodeLayers else termLayer + ixp.iterDir

        for j in range(ixp.lastNodeLayer, terminal_layer, ixp.iterDir):
            for (u, ui) in self.allnodeitems_in_layer(j):
                if j == ixp.lastNodeLayer:
                    setattr(self.ixinfo[u], ixp.type, IxTuple(u, ixp.type, 0, None, False))
                    #assert 0 == self._find_optix(ixp.type, u)[0]
                else:
                    for a in getattr(ui, ixp.nextArcs):
                        new_weight = a.weight + self._opt_ixweight(ixp.type, getattr(a, ixp.otherEnd))
                        if ixp.limCmp(new_weight, self._opt_ixweight(ixp.type, u)):
                            setattr(self.ixinfo[u], ixp.type, IxTuple(u, ixp.type, new_weight, a, True))
                    setattr(getattr(self.ixinfo[u], ixp.type), 'need_check', False)
                    #assert self._opt_ixweight(ixp.type, u) == self._find_optix(ixp.type, u)[0]

    def _update_ixinfo_bynode(self, ixType, nodeList):
        """Update IxInfo recursively from list of nodes.

        Update IxInfo recursively from nodes in 'nodeList'. If IxInfo is
        updated, all nodes in the subtree of the optimal ixtree must also
        be updated; this is done recursively. Note that 'nodeList' can only
        contain nodes from a single layer.

        Args:
            ixType (str): type of ix
            nodeList (List[MDDNode]): initial list of nodes to update

        Raises:
            ValueError: 'nodeList' cannot contain nodes from different layers
        """
        # Check all nodes in nodeList are on same layer
        if len(set(v.layer for v in nodeList)) > 1:
            raise ValueError('nodeList cannot contain nodes from different layers: %s' % str(nodeList))

        # Adjust path weights of nodes in optimal subtree
        ixp = IxParam(ixType, self.numNodeLayers)
        updateNodes = set(nodeList)
        for u in updateNodes:
            setattr(getattr(self.ixinfo[u], ixp.type), 'need_check', True)
        while len(updateNodes) > 0:
            nextUpdateNodes = set()
            for u in updateNodes:
                # If node requires a check, peform a full update
                if self._opt_ixcheck(ixp.type, u):
                    old_weight = self._opt_ixweight(ixp.type, u)
                    old_ixtuple = getattr(self.ixinfo[u], ixp.type)
                    setattr(self.ixinfo[u], ixp.type, IxTuple(u, ixp.type, ixp.limVal, None, True))
                    if u.layer == ixp.lastNodeLayer:
                        setattr(self.ixinfo[u], ixp.type, IxTuple(u, ixp.type, 0, None, False))
                        #assert 0 == self._find_optix(ixp.type, u)[0]
                    else:
                        currNextArcs = getattr(self._get_node_info(u), ixp.nextArcs)
                        for a in currNextArcs:
                            nextNode = getattr(a, ixp.otherEnd)
                            if self._opt_ixcheck(ixp.type, nextNode):
                                raise RuntimeWarning('node in closer layer not checked; performing layer update')
                                self._update_ixinfo_bylayer(ixp.type, nextNode.layer)
                            new_weight = a.weight + self._opt_ixweight(ixp.type, nextNode)
                            if ixp.limCmp(new_weight, self._opt_ixweight(ixp.type, u)):
                                setattr(self.ixinfo[u], ixp.type, IxTuple(u, ixp.type, new_weight, a, True))
                        setattr(getattr(self.ixinfo[u], ixp.type), 'need_check', False)
                        #assert self._opt_ixweight(ixp.type, u) == self._find_optix(ixp.type, u)[0]
                # Recursing down the optimal subtree, add nodes that need to be updated
                for a in getattr(self._get_node_info(u), ixp.recArcs):
                    # Inefficient but correct
                    #recNode = getattr(a, ixp.recEnd)
                    #setattr(getattr(self.ixinfo[recNode], ixp.type), 'need_check', True)
                    #nextUpdateNodes.add(recNode)

                    # More efficient, but is it correct?
                    recNode = getattr(a, ixp.recEnd)
                    curr_weight = self._opt_ixweight(ixp.type, recNode)
                    new_weight = a.weight + self._opt_ixweight(ixp.type, u)
                    if self._opt_ixcheck(ixp.type, recNode):
                        # If a node needs to be checked, add it to be fully updated
                        nextUpdateNodes.add(recNode)
                    elif ixp.limCmp(new_weight, curr_weight):
                        # If a node is checked AND new_weight is better, we can do a simple update
                        setattr(self.ixinfo[recNode], ixp.type, IxTuple(recNode, ixp.type, new_weight, a, False))
                        #assert self._opt_ixweight(ixp.type, recNode) == self._find_optix(ixp.type, recNode)[0]
                        #if self._opt_ixweight(ixp.type, recNode) == ixp.limVal:
                        #    print('\t\t!!!Node {}: Incoming = {}, Outgoing = {}'.format(recNode, self._get_node_info(recNode).incoming, self._get_node_info(recNode).outgoing))
                        # This assert won't always be true, but when the iteration ends it will be true.
                        ##setattr(getattr(self.ixinfo[recNode], ixp.type), 'need_check', True)
                        nextUpdateNodes.add(recNode)
                    elif ixp.limCmp(curr_weight, new_weight) and self._opt_ixarc(ixp.type, recNode) == a:
                        # If a node is checked AND new_weight is worse, we only need to update if its optimal arc is a
                        # If ixarc == None, something is wrong: since we're only checking nodes that have optimal arcs, ixarc should never be None here.
                        setattr(getattr(self.ixinfo[recNode], ixp.type), 'need_check', True)
                        nextUpdateNodes.add(recNode)

            updateNodes = nextUpdateNodes

    ################################
    # Incremental update functions #
    ################################
    # The following set of functions are all augmented versions of functions
    # in MDD, in order to maintain the IxInfo invariant every time the MDD
    # is changed.

    def _add_arc(self, newarc):
        #super()._add_arc(newarc)
        mdd.MDD._add_arc(self, newarc)
        #print('Adding arc {}'.format(newarc))
        for ixType in IxParam.all_ix_types:
            ixp = IxParam(ixType, self.numNodeLayers)
            newWeight = newarc.weight + self._opt_ixweight(ixp.type, getattr(newarc, ixp.otherEnd))
            currWeight = self._opt_ixweight(ixp.type, getattr(newarc, ixp.recEnd))
            if ixp.limCmp(newWeight, currWeight):
                self._update_ixinfo_bynode(ixp.type, [getattr(newarc, ixp.recEnd)])

    def _remove_arc(self, rmvarc):
        toProcess = dict()
        for ixType in IxParam.all_ix_types:
            ixp = IxParam(ixType, self.numNodeLayers)
            toProcess[ixType] = []
            optArc = self._opt_ixarc(ixp.type, getattr(rmvarc, ixp.recEnd))
            if optArc is not None and optArc == rmvarc:
                toProcess[ixType].append(getattr(rmvarc, ixp.recEnd))
        #print('Removing arc {}'.format(rmvarc))
        #super()._remove_arc(rmvarc)
        mdd.MDD._remove_arc(self, rmvarc)
        for ixType in IxParam.all_ix_types:
            self._update_ixinfo_bynode(ixType, toProcess[ixType])

    def _add_node(self, newnode):
        #super()._add_node(newnode)
        mdd.MDD._add_node(self, newnode)
        #print('Adding node {}'.format(newnode))
        self.ixinfo[newnode] = IxInfo(newnode)
        for ixType in IxParam.all_ix_types:
            ixp = IxParam(ixType, self.numNodeLayers)
            if newnode.layer == ixp.lastNodeLayer:
                setattr(self.ixinfo[newnode], ixp.type, IxTuple(newnode, ixp.type, 0, None, False))
                #assert 0 == self._find_optix(ixp.type, newnode)[0]
            else:
                setattr(getattr(self.ixinfo[newnode], ixType), 'need_check', False)
                #assert self._opt_ixweight(ixType, newnode) == self._find_optix(ixType, newnode)[0]

    def _remove_node(self, rmvnode):
        toProcess = dict()
        for ixType in IxParam.all_ix_types:
            ixp = IxParam(ixType, self.numNodeLayers)
            toProcess[ixType] = []
            # Look at all relevant arcs to be removed...
            for arc in getattr(self._get_node_info(rmvnode), ixp.recArcs):
                # If any of these arcs are optimal ix arcs for the node
                # on the opposite side, that node needs to be updated.
                optArc = self._opt_ixarc(ixp.type, getattr(arc, ixp.recEnd))
                if optArc is not None and optArc == arc:
                    #self._update_ixinfo_bynode(ixp.type, [getattr(arc, ixp.recEnd)])
                    toProcess[ixType].append(getattr(arc, ixp.recEnd))
        #print('Removing node {}'.format(rmvnode))
        #super()._remove_node(rmvnode)
        mdd.MDD._remove_node(self, rmvnode)
        for ixType in IxParam.all_ix_types:
            self._update_ixinfo_bynode(ixType, toProcess[ixType])


    def _append_new_layer(self):
        #super()._append_new_layer()
        mdd.MDD._append_new_layer(self)
        for ixType in ('min_suffix', 'max_suffix'):
            ixp = IxParam(ixType, self.numNodeLayers)
            for u in self.allnodes_in_layer(self.numArcLayers-1):
                setattr(self.ixinfo[u], ixp.type, IxTuple(u, ixp.type, ixp.limVal, None, False))
                #assert ixp.limVal == self._find_optix(ixp.type, u)[0]
            self._update_ixinfo_bynode(ixType, [u for u in self.allnodes_in_layer(self.numArcLayers-1)])

    def _clear(self):
        #super()._clear()
        mdd.MDD._clear(self)
        self.ixinfo.clear()

#    def _remove_nodes(self, rmvnodes):
#    def add_node(self, newnode):
#    def remove_node(self, rmvnode):
#    def add_arc(self, newarc):
#    def remove_arc(self, rmvarc):
#    def _merge_nodes_internal(self, mnodes, mlayer, nodefun, inarcfun=None, outarcfun=None):
#    def _merge_nodes(self, mnodes, mlayer, nsfun, awinfun=None, awoutfun=None):
#    def merge_nodes(self, mnodes, nsfun, awinfun=None, awoutfun=None):
#    def _redirect_incoming_arcs(self, old_head, new_head):
#    def prune_all(self):
#    def prune_recursive(self, origNodeList):
#    def compile_top_down(self, numLayers, domainFunc, trFunc, costFunc, rootState, isFeas, maxWidth=None, mergeFunc=None, adjFunc=None, nodeSelFunc=None):
#    def compile_trivial(self, numLayers, domainFunc, costFunc, nodeStateFunc):
#    def filter_and_refine_constraint(self, trFunc, rootState, isFeas, nodeStateFunc, maxWidth=None):
#    def compile_pathlist(self, pathList, minArcCostFunc=None, nodeStateFunc=None):
#    def reduce_bottom_up(self, mergeFunc, adjInFunc=None, adjOutFunc=None, ignoreLastLayer=False):
