# Using py.test-3
from pymdd.mdd import MDD, MDDNode
def _assert_clean_tmp(m):
    # Check tmp is clean afterwards
    for j in range(m.numNodeLayers):
        for (u, ui) in m.allnodeitems_in_layer(j):
            assert ui._tmp is None

def test_find_opt_rtpath():
    mymdd = MDD()
    mymdd.loadJSON('data/misp.json')
    assert mymdd._find_opt_rtpath(False) ==  (0.0, [0,0,0,0,0])
    _assert_clean_tmp(mymdd)
    assert mymdd._find_opt_rtpath(True) ==  (11.0, [0,1,0,0,1])
    _assert_clean_tmp(mymdd)

    mymdd.loadJSON('data/trivial.json')
    assert mymdd._find_opt_rtpath(False) ==  (0.0, [])
    _assert_clean_tmp(mymdd)
    mymdd._append_new_layer()
    mymdd.add_node(MDDNode(1, 'NEO'))
    assert mymdd.find_shortest_path() ==  (float('inf'), [])
    assert mymdd.find_longest_path() ==  (float('-inf'), [])

def test_find_optimal_ix():
    mymdd = MDD()
    mymdd.loadJSON('data/misp.json')

    testnode = MDDNode(3, frozenset([5]))
    assert mymdd._find_optimal_ix(testnode, False, False) == (2.0, [0,0,1])
    _assert_clean_tmp(mymdd)
    assert mymdd._find_optimal_ix(testnode, False, True) == (4.0, [0,1,0])
    _assert_clean_tmp(mymdd)
    assert mymdd._find_optimal_ix(testnode, True, False) == (0.0, [0,0])
    _assert_clean_tmp(mymdd)
    assert mymdd._find_optimal_ix(testnode, True, True) == (7.0, [0,1])
    _assert_clean_tmp(mymdd)

    assert mymdd._find_optimal_ix(MDDNode(0, frozenset([1,2,3,4,5])), False, True) == (0.0, [])
    _assert_clean_tmp(mymdd)
    assert mymdd._find_optimal_ix(MDDNode(5, frozenset()), True, False) == (0.0, [])
    _assert_clean_tmp(mymdd)

    mymdd.loadJSON('data/trivial.json')
    assert mymdd._find_optimal_ix(MDDNode(0, ''), False, True) == (0.0, [])
    _assert_clean_tmp(mymdd)
    assert mymdd._find_optimal_ix(MDDNode(0, ''), True, False) == (0.0, [])
    _assert_clean_tmp(mymdd)


def test_find_path():
    mymdd = MDD()
    mymdd.loadJSON('data/misp.json')

    srcList = list(mymdd.allnodes_in_layer(4))
    destList = list(mymdd.allnodes_in_layer(1))
    assert mymdd._find_optimal_path(srcList, destList, False) == (0.0, [0,0,0])
    _assert_clean_tmp(mymdd)
    assert mymdd._find_optimal_path(srcList, destList, True) == (4.0, [0,0,1])
    _assert_clean_tmp(mymdd)

    srcList = [MDDNode(0, frozenset([1,2,3,4,5]))]
    destList = [MDDNode(5, frozenset())]
    assert mymdd._find_optimal_path(srcList, destList, False) == (0.0, [0,0,0,0,0])
    _assert_clean_tmp(mymdd)
    assert mymdd._find_optimal_path(srcList, destList, True) == (11.0, [0,1,0,0,1])
    _assert_clean_tmp(mymdd)

    srcList = [MDDNode(2, frozenset([5])), MDDNode(2, frozenset([3,4,5]))]
    destList = [MDDNode(2, frozenset([4,5]))]
    assert mymdd._find_optimal_path(srcList, destList, True) == (float('-inf'), [])
    _assert_clean_tmp(mymdd)
    destList.append(MDDNode(2, frozenset([3,4,5])))
    assert mymdd._find_optimal_path(srcList, destList, True) == (0.0, [])
    _assert_clean_tmp(mymdd)
