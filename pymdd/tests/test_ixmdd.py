# Using py.test-3
from pymdd.mdd import MDDNode, MDDArc
from pymdd.ixmdd import IxMDD, IxParam

def _compile_topdown_misp_ixmdd():
    # MISP instance
    G = [[(1,2), (1,3)], [(2,1), (2,3), (2,4)], [(3,1), (3,2), (3,4)], [(4,2), (4,3), (4,5)], [(5,4)]]
    neighbors = [ [v for (u,v) in G[j]] for j in range(5) ]
    weight = [3, 4, 2, 2, 7]

    # DP model
    numLayers = 5
    domain = lambda j: (0,1)
    def trFunc(s,d,j):
        if d == 1:
            if j+1 in s:
                return s - {j+1} - set(neighbors[j])
            else:
                return None
        else:
            return s - {j+1}
    costFunc = lambda s,d,j,ns: d*weight[j]
    rootState = frozenset([1,2,3,4,5])
    isFeas = lambda s,j: s is not None
    name = 'misp'

    # Compile MDD
    mymdd = IxMDD(name)
    mymdd.compile_top_down(numLayers, domain, trFunc, costFunc, rootState, isFeas)

    return mymdd

def test_ixmdd_initial():
    mymdd = IxMDD()
    assert len(mymdd.ixinfo) == 0

    mymdd._append_new_layer()
    rootNode = MDDNode(0, 0)
    mymdd._add_node(rootNode)
    opt_weight = { rootNode: {'min_suffix': 0, 'max_suffix': 0, 'min_prefix': 0, 'max_prefix': 0} }
    for ixType in IxParam.all_ix_types:
        assert not mymdd._opt_ixcheck(ixType, rootNode)
        assert mymdd._opt_ixweight(ixType, rootNode) == opt_weight[rootNode][ixType]

    mymdd._append_new_layer()
    opt_weight = { rootNode: {'min_suffix': float('inf'), 'max_suffix': float('-inf'), 'min_prefix': 0, 'max_prefix': 0} }
    for ixType in IxParam.all_ix_types:
        assert not mymdd._opt_ixcheck(ixType, rootNode)
        assert mymdd._opt_ixweight(ixType, rootNode) == opt_weight[rootNode][ixType]

    secondNode = MDDNode(1, 1)
    mymdd._add_node(secondNode)
    assert len(mymdd.ixinfo) == 2
    opt_weight = { rootNode: {'min_suffix': float('inf'), 'max_suffix': float('-inf'), 'min_prefix': 0, 'max_prefix': 0}, secondNode: {'min_suffix': 0, 'max_suffix': 0, 'min_prefix': float('inf'), 'max_prefix': float('-inf')} }
    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert not mymdd._opt_ixcheck(ixType, v)
            assert mymdd._opt_ixweight(ixType, v) == opt_weight[v][ixType]

    mymdd._add_arc(MDDArc(1, -2, rootNode, secondNode))
    mymdd._add_arc(MDDArc(2, 5, rootNode, secondNode))
    opt_weight = { rootNode: {'min_suffix': -2, 'max_suffix': 5, 'min_prefix': 0, 'max_prefix': 0}, secondNode: {'min_suffix': 0, 'max_suffix': 0, 'min_prefix': -2, 'max_prefix': 5} }
    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert not mymdd._opt_ixcheck(ixType, v)
            assert mymdd._opt_ixweight(ixType, v) == opt_weight[v][ixType]

def test_ixmdd_basic_json():
    mymdd = IxMDD()
    mymdd.loadJSON('data/misp.json')

    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert not mymdd._opt_ixcheck(ixType, v)

def test_ixmdd_basic_compile_topdown():
    mymdd = _compile_topdown_misp_ixmdd()

    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert not mymdd._opt_ixcheck(ixType, v)

def test_ixmdd_basic_filter_refine():
    # MISP instance
    G = [[(i,j) for j in range(1,6) if i != j] for i in range(1,6)]
    #G = [[(1,2), (1,3)], [(2,1), (2,3), (2,4)], [(3,1), (3,2), (3,4)], [(4,2), (4,3), (4,5)], [(5,4)]]
    neighbors = [ [v for (u,v) in G[j]] for j in range(5) ]
    #weight = [3, 4, 2, 2, 7]
    weight = [1, 1, 1, 1, 1]

    # DP model
    numLayers = 5
    domainFunc = lambda j: (0,1)
    def trFunc(s,d,j):
        if d == 1:
            if j+1 in s:
                return s - {j+1} - set(neighbors[j])
            else:
                return None
        else:
            return s - {j+1}
    costFunc = lambda d,j: d*weight[j]
    rootState = frozenset([1,2,3,4,5])
    isFeas = lambda s,j: s is not None
    name = 'smisp'
    def hack_counter(k = [-1]):
        k[0] += 1
        return k[0]

    mymdd = IxMDD(name)
    mymdd.compile_trivial(numLayers, domainFunc, costFunc, lambda j: hack_counter())
    mymdd.filter_and_refine_constraint(trFunc, rootState, isFeas, lambda s, j: hack_counter())
    mymdd.reduce_bottom_up(lambda slist, j: slist[0])

    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert not mymdd._opt_ixcheck(ixType, v)

def test_ixmdd_update_all():
    mymdd = _compile_topdown_misp_ixmdd()
    mymdd.update_ixinfo_all()

    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert not mymdd._opt_ixcheck(ixType, v)

def test_find_optix():
    allNodes = [MDDNode(0, frozenset([1,2,3,4,5])), MDDNode(1, frozenset([2,3,4,5])), MDDNode(1, frozenset([4,5])), MDDNode(2, frozenset([5])), MDDNode(2, frozenset([3,4,5])), MDDNode(2, frozenset([4,5])), MDDNode(3, frozenset([5])), MDDNode(3, frozenset([4,5])), MDDNode(4, frozenset([5])), MDDNode(4, frozenset()), MDDNode(5, frozenset())]
    minSufKey = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    maxSufKey = [11, 11, 7, 7, 9, 7, 7, 7, 7, 0, 0]
    minPreKey = [0, 0, 3, 4, 0, 3, 2, 0, 0, 2, 0]
    maxPreKey = [0, 0, 3, 4, 0, 3, 4, 3, 4, 5, 11]

    mymdd = _compile_topdown_misp_ixmdd()
    for i in range(len(allNodes)):
        assert mymdd._opt_ixweight('min_suffix', allNodes[i]) == minSufKey[i]
    for i in range(len(allNodes)):
        assert mymdd._opt_ixweight('max_suffix', allNodes[i]) == maxSufKey[i]
    for i in range(len(allNodes)):
        assert mymdd._opt_ixweight('min_prefix', allNodes[i]) == minPreKey[i]
    for i in range(len(allNodes)):
        assert mymdd._opt_ixweight('max_prefix', allNodes[i]) == maxPreKey[i]

def test_ixmdd_invariants():
    mymdd = IxMDD()
    mymdd.loadJSON('data/misp.json')
    mymdd.update_ixinfo_all()

    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert mymdd._opt_ixweight(ixType, v) == mymdd._find_optix(ixType, v)[0]

    mymdd = IxMDD()
    mymdd.loadJSON('data/misp.json')
    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert mymdd._opt_ixweight(ixType, v) == mymdd._find_optix(ixType, v)[0]

    mymdd = _compile_topdown_misp_ixmdd()
    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert mymdd._opt_ixweight(ixType, v) == mymdd._find_optix(ixType, v)[0]

def test_ixmdd_invariants2():
    # MISP instance
    G = [[(i,j) for j in range(1,6) if i != j] for i in range(1,6)]
    #G = [[(1,2), (1,3)], [(2,1), (2,3), (2,4)], [(3,1), (3,2), (3,4)], [(4,2), (4,3), (4,5)], [(5,4)]]
    neighbors = [ [v for (u,v) in G[j]] for j in range(5) ]
    #weight = [3, 4, 2, 2, 7]
    weight = [1, 1, 1, 1, 1]

    # DP model
    numLayers = 5
    domainFunc = lambda j: (0,1)
    def trFunc(s,d,j):
        if d == 1:
            if j+1 in s:
                return s - {j+1} - set(neighbors[j])
            else:
                return None
        else:
            return s - {j+1}
    costFunc = lambda d,j: d*weight[j]
    rootState = frozenset([1,2,3,4,5])
    isFeas = lambda s,j: s is not None
    name = 'smisp'
    def hack_counter(k = [-1]):
        k[0] += 1
        return k[0]

    mymdd = IxMDD(name)
    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert mymdd._opt_ixweight(ixType, v) == mymdd._find_optix(ixType, v)[0]

    mymdd.compile_trivial(numLayers, domainFunc, costFunc, lambda j: hack_counter())
    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert mymdd._opt_ixweight(ixType, v) == mymdd._find_optix(ixType, v)[0]

    mymdd.filter_and_refine_constraint(trFunc, rootState, isFeas, lambda s, j: hack_counter())
    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert mymdd._opt_ixweight(ixType, v) == mymdd._find_optix(ixType, v)[0]

    mymdd.reduce_bottom_up(lambda slist, j: slist[0])
    for v in mymdd.allnodes():
        for ixType in IxParam.all_ix_types:
            assert mymdd._opt_ixweight(ixType, v) == mymdd._find_optix(ixType, v)[0]
