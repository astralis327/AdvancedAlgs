#!/usr/bin/python3
import copy
import math
import queue

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
import heapq
import itertools



class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''
	
	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for 
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

	def greedy( self,time_allowance=60.0 ):
		results = {}
		start_time = time.time()
		cities = self._scenario.getCities()

		route = []
		unvisitedCities = copy.deepcopy(cities)

		fullTour = False
		nextCityIndex = 0
		numSolns = 0
		city = unvisitedCities.pop(0)
		startCity = city
		totalDist = 0

		while(not fullTour and time.time() - start_time < time_allowance):
			lowestCost = math.inf
			closestCity = None
			if len(unvisitedCities) == 0:
				disToStart = city.costTo(startCity)
				fullTour = True
				route.append(startCity)
				totalDist += disToStart
				continue
			for index in range(len(unvisitedCities)):

				distance = city.costTo(unvisitedCities[index])
				if distance < lowestCost:
					closestCity	= unvisitedCities[index]
					lowestCost = distance
					city = closestCity
					nextCityIndex = index			# save the index of closest city to delete from unvisited cities

			route.append(closestCity)
			del unvisitedCities[nextCityIndex]
			totalDist += lowestCost					# add the latest distance to the total
		foundTour = False
		if totalDist == math.inf:
			print("No Solution to GREEDY! Using a RandomTour.")			#if no soln is to be found with the given config of cities
			results = self.defaultRandomTour()
			return results

		bssf = TSPSolution(route)
		bssf.cost = totalDist

		numSolns += 1

		if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = numSolns
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results



	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''
		
	def bab_init( self, time_allowance=60.0 ):
		start_time = time.time()
		results = {}
		S = queue.PriorityQueue()					#priority queue to save the current active problems
		cities = self._scenario.getCities()
		ncities = len(cities)
		bssf = self.greedy()['soln']
		if bssf == False:
			bssf = self.defaultRandomTour()
		costMatrix = []
		for row in range(ncities): 				#make the initial cost matrix
			currentRow = []
			city = cities[row]
			for col in range(ncities):
				cityToBeEntered = cities[col]
				if row == col: cost = math.inf		# make the diagonal cost infinity
				else:
					cost = city.costTo(cityToBeEntered)
				currentRow.append(cost)
			costMatrix.append(currentRow)
		start_node = self.lowerbound(costMatrix, [], 0, None, start_time, time_allowance) 	#reduce the cost matrix and find initial lowerbound
		S.put((start_node[0],(start_node[1], [0])))											#put this initial matrix into the pq
		sol_route, sol_cost, count, num_pruned, max = \
			self.branch_and_bound(S,bssf, start_time, time_allowance)

		if sol_cost < bssf.cost:
			for i in range(len(sol_route)):

				bssf.route[i] = cities[sol_route[i]]

			bssf.route.append(cities[0])
			bssf.cost = sol_cost



		end_time = time.time()
		results['cost'] = bssf.cost
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = max
		results['total'] = max + num_pruned
		results['pruned'] = num_pruned
		return results


	def branch_and_bound(self, active_probs, bssf, start_time, time_allowance):
		route = bssf.route			#initial route found by greedy or random
		cost = bssf.cost
		max_num_subproblems = 0		#maximum number of subproblems to exist during runtime

		count_solns = 0				#num solns found
		num_pruned = 0				#num subproblems pruned
		first_pass = True
		while active_probs.qsize() > 0 	and time.time() - start_time < time_allowance:
			if active_probs.qsize() > max_num_subproblems: max_num_subproblems = active_probs.qsize()
			active_node = active_probs.get()
			lb = active_node[0]						#lb of the current subproblem
			cost_Matrix = active_node[1][0]
			tour = active_node[1][1]				#tour so far in the current subproblem

			if len(tour) > 0:
				paths = cost_Matrix[tour[-1]] 		#possible cities to go next
			else: paths = cost_Matrix[0]
			for col in range(len(paths)):

				if paths[col] == math.inf or (col == 0 and not first_pass):		#Skip over column if the cost is infinity or it's the first column, since we won't go back there until a valid tour's been found
					continue
				else:
					#make a new branch and reduece the resulting matrix
					new_lb, branch_cost = self.branch(col, copy.deepcopy(cost_Matrix), tour, copy.deepcopy(lb), start_time, time_allowance)

					#if the time is up return the current bssf
					if branch_cost == True: return bssf.route, bssf.cost, count_solns, num_pruned, max_num_subproblems
					if new_lb == math.inf or new_lb >= cost:		#if the newest lb is inf or larger than bssf, prune
						num_pruned += 1
						continue
					else:
						tour.append(col)		#add the city

						active_probs.put((new_lb,(branch_cost, copy.deepcopy(tour))))
						if len(tour) == len(cost_Matrix[0]): #if the length of the tour is the number of cities, valid soln
							route = tour
							cost = lb
							count_solns += 1
				del tour[-1]
			first_pass = False		#once the first city has been passed
		return route, cost, count_solns, num_pruned, max_num_subproblems

	def branch(self, column, branch_matrix, cities_visited, lb_parent, start_time, time_allowance):
		#set the row and column of the current city to inf

		dist = branch_matrix[cities_visited[-1]][column]
		for col2 in range(len(branch_matrix[cities_visited[-1]])):
			branch_matrix[cities_visited[-1]][col2] = math.inf
		branch_matrix[column][cities_visited[-1]] = math.inf
		for row in range(len(branch_matrix)):
			branch_matrix[row][column] = math.inf
		cost= self.calculate_cost(dist, branch_matrix, cities_visited, lb_parent, column, start_time, time_allowance)
		return cost


	def calculate_cost(self, distance, cost_matrix, cities_visited, lb_parent, current_city, start_time, time_allowance):
		residual_cost = self.lowerbound(cost_matrix, cities_visited, lb_parent, current_city, start_time, time_allowance)
		return residual_cost[0] + distance, residual_cost[1]

	def lowerbound(self, distMatrix, path, lb_parent, current_city, start_time, time_allowance):
		distMatrix, lb = self.reduceByRow(distMatrix, path, lb_parent,start_time, time_allowance)
		if lb == True: return 0, True	#if the time is up, return to the outer while loop
		distMatrix, lb = self.reduceByCol(distMatrix, path, lb, current_city,start_time, time_allowance)
		if lb == True: return 0, True
		return lb, distMatrix

	def reduceByRow(self, matrix, path, lb, start_time, time_allowance):
		for r in range(len(matrix)):
			if (time.time() - start_time) > time_allowance:
				return [0, True]
			minInRow = min(matrix[r])

			seen = self.city_already_visited(path, r)		#if city already visited
			if seen is not False:
				continue
			else:
				lb += minInRow								#update lb
				for c in range(0, len(matrix[r])):
						newCellVal = matrix[r][c] - minInRow	# update matrix
						matrix[r][c] = newCellVal
		return matrix, lb

	def reduceByCol(self, matrix, visited_cities, lb, current_city,start_time, time_allowance):
		for c in range(len(matrix)):
			if (time.time() - start_time) > time_allowance:
				return [0, True]
			minCol = math.inf				#min val in column
			for r in range(len(matrix)):
				if minCol == 0: break		#min isn't going any lower, so break
				if matrix[r][c] < minCol:
					minCol = matrix[r][c]
			seen = self.city_already_visited(visited_cities, c)	#if city already visited
			if ( seen or current_city == c ) and \
					(minCol == math.inf) or minCol ==0: continue	#skip over cities already visited are are currently being visited, or if the min col = 0 and so the values aren't going to be changed anyway
			else:
				for r in range(len(matrix)):
					matrix[r][c] -= minCol
				lb += minCol
		return matrix, lb

	def city_already_visited(self, visited_cities, city):
		return city in visited_cities

	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''
		
	def fancy( self,time_allowance=60.0 ):
		pass
		


