# Using py.test-3
from pymdd.mdd import MDD, MDDNode
def _assert_clean_tmp(m):
    # Check tmp is clean afterwards
    for j in range(m.numNodeLayers):
        for (u, ui) in m.allnodeitems_in_layer(j):
            assert ui._tmp is None

def test_enumerate_all_paths():
    mymdd = MDD()
    mymdd.loadJSON('data/misp.json')
    assert len(mymdd.enumerate_all_paths()) == 10

def test_enumerate_paths():
    mymdd = MDD()
    mymdd.loadJSON('data/misp.json')

    testnode = MDDNode(2, frozenset([3,4,5]))
    assert sorted(mymdd.enumerate_all_suffixes(testnode)) == [(0, [0,0,0]), (2, [0,1,0]), (2, [1,0,0]), (7, [0,0,1]), (9, [1,0,1])]
    _assert_clean_tmp(mymdd)
    testnode = MDDNode(4, frozenset([5]))
    assert sorted(mymdd.enumerate_all_prefixes(testnode)) == [(0, [0,0,0,0]), (2, [0,0,1,0]), (3, [1,0,0,0]), (4, [0,1,0,0])]
    _assert_clean_tmp(mymdd)
    
    mymdd.loadJSON('data/trivial.json')
    assert sorted(mymdd.enumerate_all_suffixes(MDDNode(0, ''))) == []
    _assert_clean_tmp(mymdd)
