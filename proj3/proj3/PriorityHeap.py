import math


class PriorityHeap:
    heapSize = 0
    pointer = None
    bHeap = []
    priorities = []

    def __init__(self, vertices, priorities, src_node):
        # makeQueue method, and the loop in dijkstra is included here. Complexity is therefore O(V)
        self.pointer = [-1] * (len(vertices) + 1)  # size of pointer list should be the length of the greatest
        # element in the vertex set
        self.src_node = src_node
        self.bHeap = vertices
        self.heapSize = len(self.bHeap)
        self.priorities = priorities
        i = 0
        for i in range(len(self.bHeap)):
            #if src node is encountered, move it to root
            if self.priorities[i] == 0:
                temp = self.bHeap[i]
                self.bHeap[i] = self.bHeap[0]
                self.pointer[self.bHeap[0].node_id] = i
                self.bHeap[0] = temp
                self.pointer[temp.node_id] = 0

                temp = self.priorities[i]
                self.priorities[i] = self.priorities[0]
                self.priorities[0] = temp

            self.pointer[self.bHeap[i].node_id] = i

    def decreaseKey(self, currentNode, newPriority):
        # All of the nodes and costs are organized by a tree, so it takes a maximum of log_2 (|V|) to bubble up
        self.priorities[currentNode] = newPriority  # set the node with the new cost
        bubble = True
        parent = 0
        while bubble and currentNode > 0:  # while the cost of the parent node is less than the current
            parent = math.floor((currentNode - 1) / 2)    #parent node
            bubble = self.bubbleDown(currentNode, parent)   #bubble the child node up
            currentNode = parent

    def deleteMin(self):
        #This method is also log|V| because we're only looking
        #at one branch of the tree at each stage, which means that the most elements that are looked at is the depth of
        #the tree, or log|V|.

        min = self.bHeap[0] #root
        lowestPriority = self.priorities[0]
        if self.heapSize > 1:
            self.bHeap[0] = self.bHeap.pop(self.heapSize - 1)   #move the last node to the root
            self.heapSize -= 1  #decrement heap size
            self.priorities[0] = self.priorities.pop(self.heapSize - 1) #move the last priorities element to the first
            self.pointer[self.bHeap[0].node_id] = 0     #update root pointer
            self.pointer[0] = -1  #update old root pointer
            bubble = True
            i = 0
            while i < len(self.bHeap) and bubble:
                child1 = 2 * i + 1   #children
                child2 = (2 * i + 1) + 1

                #compare the two children to find which is smallest

                #Make sure that the children don't go out of bounds of the tree
                if child1 < len(self.priorities) - 1 and child2 <= len(self.priorities) - 1:
                    #Case 1: cost of child1 < cost of child2
                    if self.priorities[child1] < self.priorities[child2]:
                        #Switch between child1 and the new root
                        bubble = self.bubbleDown(child1, i)
                        i = child1
                    else:
                        #switch between child2 and the new root
                        bubble = self.bubbleDown(child2, i)
                        i = child2
                #If there's only one child, switch with that one
                elif child1 < len(self.priorities) and child2 >= len(self.priorities):
                    bubble = self.bubbleDown(child1, i)
                    i = child1
                else:
                    i = child1
        return min, lowestPriority

    def bubbleDown(self, switchNode, currentNode):
        #compare current node and node being switched. This method runs in constant time.
        if self.priorities[switchNode] < self.priorities[currentNode]:
            #update arrays
            temp = self.bHeap[switchNode]
            self.bHeap[switchNode] = self.bHeap[currentNode]
            self.pointer[self.bHeap[currentNode].node_id] = switchNode
            self.bHeap[currentNode] = temp
            self.pointer[temp.node_id] = currentNode

            temp = self.priorities[switchNode]
            self.priorities[switchNode] = self.priorities[currentNode]
            self.priorities[currentNode] = temp
            return True
        else:
            return False
