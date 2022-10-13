import math


class PriorityArray:

    def __init__(self,vertices, priorities):
        #insert method, O(1) complexity
        self.priority_queue = priorities
        self.size = len(vertices)
        self.visited = [False]*len(vertices)

    def decreaseKey(self, node, newCost):
        #constant complexity
        self.priority_queue[node] = newCost

    def deleteMin(self):
        #O(|V|) complexity
        i = 0
        min_index = 0
        min = 0
        for k in range(len(self.visited)):
            if  self.visited[k] == False:
                min = self.priority_queue[k]
                break
        for i in range(len(self.priority_queue)):
            if self.priority_queue[i] < min and not self.visited[i]:
                min_index = i
                min = self.priority_queue[i]
        self.visited[min_index] = True
        min_cost = self.priority_queue[min_index]
        self.size -= 1
        return min_index, min_cost
