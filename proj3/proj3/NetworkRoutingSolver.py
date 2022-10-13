#!/usr/bin/python3
import math
from copy import copy

from CS4412Graph import *
import time
from PriorityHeap import PriorityHeap
from PriorityArray import PriorityArray


class NetworkRoutingSolver:
    def __init__( self):
        pass

    def initializeNetwork( self, network ):
        assert( type(network) == CS4412Graph )

        self.network = network
        self.dist = [math.inf]*len(self.network.nodes)
        self.prev = [None]*len(self.network.nodes)


    def getShortestPath( self, destIndex ):
        self.dest = destIndex
        # TODO: RETURN THE SHORTEST PATH FOR destIndex
        #       INSTEAD OF THE DUMMY SET OF EDGES BELOW
        #       IT'S JUST AN EXAMPLE OF THE FORMAT YOU'LL 
        #       NEED TO USE
        path_edges = []
        total_length = 0

        #destination
        nodeIndex = self.dest
        #iterate until source node is reached
        while self.prev[nodeIndex] != None:
            #previous node closest to source
            prevNodeIndex = self.prev[nodeIndex]
            prevNode = self.network.nodes[prevNodeIndex]
            src_edge = prevNode.neighbors[0]
            #find the edge from the neighbor to the current node
            for edge in prevNode.neighbors:
                if edge.dest.node_id == nodeIndex:
                    src_edge = edge
                    break

            path_edges.insert(0, (src_edge.src.loc, src_edge.dest.loc, '{:.0f}'.format(src_edge.length)) )
            total_length += src_edge.length
            nodeIndex = self.prev[nodeIndex]
        return {'cost':total_length, 'path':path_edges}

    def computeShortestPaths( self, srcIndex, use_heap=False ):
        self.source = srcIndex
        t1 = time.time()
        # TODO: RUN DIJKSTRA'S TO DETERMINE SHORTEST PATHS.
        #       ALSO, STORE THE RESULTS FOR THE SUBSEQUENT
        #       CALL TO getShortestPath(dest_index)
        #initialize source distance to 0
        self.dist[self.source] = 0
        pq_size = 0
        if use_heap:
            pq = PriorityHeap(copy(self.network.nodes), copy(self.dist), self.source)
            pq_size = pq.heapSize
        elif not use_heap:
            pq = PriorityArray(copy(self.network.nodes), copy(self.dist))
            pq_size = pq.size

        while pq_size != 0:
            #deleteMin gets called |V| times, DecreaseKey is log|V| for binaryHeap and is called
            #|E| times, and Insert is coupled together with the constructor of both the array and
            #the heap implementations. 


            #delete the minimum from the pq
            u = pq.deleteMin()
            if not use_heap:
                #array structure is different
                min = self.network.nodes[u[0]]
            else: min = u[0]
            #decrement size
            pq_size -=1
            #lowest cost
            lowestDist = u[1]
            for edge in min.neighbors:  #iterate through the minimum node neighbors
                v = edge.dest.node_id
                if self.dist[v] > lowestDist + edge.length:
                    self.dist[v] = lowestDist + edge.length
                    self.prev[v] = min.node_id
                    if not use_heap: index = v
                    else: index = pq.pointer[v]
                    pq.decreaseKey(index, self.dist[v])

        t2 = time.time()
        return (t2-t1)

