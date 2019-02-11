from abc import ABCMeta, abstractmethod

class DPModel(object):
    __metaclass__ = ABCMeta

    @property
    @abstractmethod
    def number_of_layers(self):
        pass

    @property
    @abstractmethod
    def root_state(self):
        pass

    @abstractmethod
    def domain_function(self, layer):
        pass

    @abstractmethod
    def transition_function(self, state, domain_value, layer):
        pass

    @abstractmethod
    def weight_function(self, state, domain_value, layer, next_state):
        pass

    @abstractmethod
    def is_feasible(self, state, layer):
        pass

    def node_selection_function(self, node_list, layer):
        return NotImplemented

    def state_merge_function(self, state_list, layer):
        return NotImplemented

    def adjusted_weight_function(self, old_weight, old_head_state, merged_head_state, layer):
        return old_weight

    def buffer_root_state(self):
        return NotImplemented

    def buffer_function(self, state, transition_list, layer):
        return NotImplemented

class KnapsackDP(DPModel):

    def __init__(self, numItems, capacity, profit, weight):
        self.numItems = numItems
        self.capacity = capacity
        self.profit = profit
        self.weight = weight

    @property
    def number_of_layers(self):
        return self.numItems

    @property
    def root_state(self):
        return self.capacity

    def domain_function(self, i):
        return (0,1)

    def transition_function(self, s, d, i):
        return s - d*self.weight[i]

    def weight_function(self, s, d, i, ns):
        return d*self.profit[i]

    def is_feasible(self, s, i):
        return s >= 0

    def state_merge_function(self, slist, i):
        return min(slist)

class MaxIndepSetDP(DPModel):

    def __init__(self, adjacencyMatrix, nodeWeight):
        self.numNodes = len(nodeWeight)
        self.weight = nodeWeight
        self.neighbors = [ [k for k in range(self.numNodes) if adjacencyMatrix[j][k] > 0] for j in range(self.numNodes) ]

    @property
    def number_of_layers(self):
        return self.numNodes

    @property
    def root_state(self):
        return frozenset(range(self.numNodes))

    def domain_function(self, j):
        return (0,1)

    def transition_function(self, s, d, j):
        if d == 1:
            if j in s:
                return s - {j} - set(self.neighbors[j])
            else:
                return None
        else:
            return s- {j}

    def weight_function(self, s, d, j, ns):
        return d*self.weight[j]

    def is_feasible(self, s, j):
        return s is not None

    def node_selection_function(self, vlist, j):
        return [min(vlist), max(vlist)]

    def state_merge_function(self, slist, j):
        return frozenset(set.union(*(set(s) for s in slist)))

    @property
    def buffer_root_state(self):
        return 0

    def buffer_function(self, state, transition_list, layer):
        return min(pst + a.weight for (pst, a) in transition_list)

class MaxCutDP(DPModel):
    
    def __init__(self, numVertices, arcWeight):
        self.numVertices = numVertices
        self.weight = arcWeight
        self.rootValue = sum(min(arcWeight[i][j], 0) for i in range(numVertices) for j in range(numVertices) if i < j)

    @property
    def number_of_layers(self):
        return self.numVertices

    @property
    def root_state(self):
        return tuple(0 for j in range(self.numVertices))

    def domain_function(self, k):
        if k == 0:
            return ['S']
        else:
            return ['S', 'T']

    def transition_function(self, s, d, k):
        if d == 'S':
            return tuple(s[l] + self.weight[k][l] if l > k else 0 for l in range(self.numVertices))
        else:
            return tuple(s[l] - self.weight[k][l] if l > k else 0 for l in range(self.numVertices))

    def weight_function(self, s, d, k, ns):
        if k == 0:
            return self.rootValue
        else:
            if d == 'S':
                return max(-s[k], 0) + sum(min(abs(s[l]), abs(self.weight[k][l])) for l in range(k+1, self.numVertices) if s[l]*self.weight[k][l] <= 0)
            else:
                return max(s[k], 0) + sum(min(abs(s[l]), abs(self.weight[k][l])) for l in range(k+1, self.numVertices) if s[l]*self.weight[k][l] >= 0)

    def is_feasible(self, s, k):
        return True

    def state_merge_function(self, slist, k):
        newstate = []
        for l in range(self.numVertices):
            if l < k:
                newstate.append(0)
            elif all(u[l] >= 0 for u in slist):
                newstate.append(min( u[l] for u in slist ))
            elif all(u[l] <= 0 for u in slist):
                newstate.append(-min( abs(u[l]) for u in slist ))
            else:
                newstate.append(0)
        return tuple(newstate)

    def adjusted_weight_function(self, w, os, ms, k):
        return w + sum(abs(os[l]) - abs(ms[l]) for l in range(k, self.numVertices))
