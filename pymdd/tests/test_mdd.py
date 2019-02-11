# Using py.test-3
from pymdd.mdd import MDD, MDDNode, MDDArc

def test_filter():
    G = [[(1,2), (1,3)], [(2,1), (2,3)], [(3,1), (3,2)]]
    neighbors = [ [v for (u,v) in G[j]] for j in range(3) ]
    weight = [1, 1, 1]
    numLayers = 3
    domain = lambda j: (0,1)

    rootState = frozenset([1,2,3])
    def trFunc(s,d,j):
        if d == 1:
            if j+1 in s:
                return s - {j+1} - set(neighbors[j])
            else:
                return None
        else:
            return s - {j+1}
    costFunc = lambda s,d,j,ns: d*weight[j]
    isFeas = lambda s,j: s is not None
    maxWidth = lambda j: 100
    name = 'smisp'
    def hack_counter(k = [-1]):
        k[0] += 1
        return k[0]

    # Construct the MDD
    mymdd = MDD(name)

    mymdd.compile_trivial(numLayers, domain, lambda d,j: costFunc(None,d,j,None), lambda j: hack_counter())
    assert mymdd.numArcLayers == numLayers

    mymdd.filter_and_refine_constraint(trFunc, rootState, isFeas, lambda s,j: hack_counter(), maxWidth)
    mymdd.reduce_bottom_up(lambda slist,j: slist[0])
    assert mymdd.widthList == [1,2,2,1]
    assert mymdd.find_shortest_path()[0] == 0
    lpath = mymdd.find_longest_path()
    assert lpath[0] == 1
    assert lpath[1] in ([1, 0, 0], [0, 1, 0], [0, 0, 1])

def test_pathlist():
    feasSolutions = [(0,1,0,1), (0,1,1,0), (0,1,1,1), (1,0,1,1), (1,1,0,0), (1,1,0,1), (1,1,1,0), (1,1,1,1)]
    cost = [6, 7, 8, 5, 6, 8, 7, 9]

    # MDD Parameters
    mergeFunc = lambda slist,i: min(slist)
    name = 'cac'

    # Construct the MDD
    mymdd = MDD(name)
    # Perform DP-based top-down compilation
    pathList = list(zip(cost, feasSolutions))
    mymdd.compile_pathlist(pathList)
    assert mymdd.widthList == [1,2,3,5,8]
    mymdd.reduce_bottom_up(mergeFunc)
    assert mymdd.widthList == [1,2,3,3,1]
    assert mymdd.find_shortest_path() == (5, [1,0,1,1])
    assert mymdd.find_longest_path() == (9, [1,1,1,1])
