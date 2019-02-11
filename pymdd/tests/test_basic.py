# Using py.test-3
from pymdd.mdd import MDD, MDDNode

def test_initial():
    mymdd = MDD()
    assert mymdd.numNodeLayers == 0
    assert mymdd.numArcLayers == -1
    assert list(mymdd.allnodes()) == []

    mymdd._append_new_layer()
    mymdd._add_node(MDDNode(0, 0))
    mymdd._append_new_layer()
    assert mymdd.widthList == [1,0]
    mymdd._add_node(MDDNode(1, 1))
    mymdd._add_node(MDDNode(1, 2))
    assert mymdd.numNodeLayers == 2
    assert mymdd.widthList == [1,2]

def test_properties():
    mymdd = MDD()
    mymdd.loadJSON('data/misp.json')

    assert mymdd.numNodeLayers == 6
    assert mymdd.numArcLayers == 5
    assert mymdd.widthList == [1, 2, 3, 2, 2, 1]
    assert mymdd.maxWidth == 3

    mymdd.loadJSON('data/trivial.json')
    assert mymdd.numNodeLayers == 1
    assert mymdd.numArcLayers == 0
    assert mymdd.widthList == [1]
    assert mymdd.maxWidth == 1
