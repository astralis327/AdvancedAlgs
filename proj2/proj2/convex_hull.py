from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF, QObject
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF, QObject
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time

# Some global color constants that might be useful
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)

# Global variable that controls the speed of the recursion automation, in seconds
#
PAUSE = 0.25


#
# This is the class you have to complete.
#
class ConvexHullSolver(QObject):

    # Class constructor
    def __init__(self):
        super().__init__()
        self.pause = False

    # Some helper methods that make calls to the GUI, allowing us to send updates
    # to be displayed.

    def showTangent(self, line, color):
        self.view.addLines(line, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseTangent(self, line):
        self.view.clearLines(line)

    def blinkTangent(self, line, color):
        self.showTangent(line, color)
        self.eraseTangent(line)

    def showHull(self, polygon, color):
        self.view.addLines(polygon, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseHull(self, polygon):
        self.view.clearLines(polygon)

    def showText(self, text):
        self.view.displayStatusText(text)

    # This is the method that gets called by the GUI and actually executes
    # the finding of the hull
    def compute_hull(self, points, pause, view):
        self.pause = pause
        self.view = view
        assert (type(points) == list and type(points[0]) == QPointF)

        t1 = time.time()
        # TODO: SORT THE POINTS BY INCREASING X-VALUE

        points_x = self.quick_sort1(points, 0, len(points) - 1)

        t2 = time.time()

        t3 = time.time()
        # this is a dummy polygon of the first 3 unsorted points
        hull_points = self.divide_and_conquer(points_x)
        polygon = [QLineF(hull_points[i], points[(i+1) % len(hull_points)]) for i in range(len(hull_points))]
        # TODO: REPLACE THE LINE ABOVE WITH A CALL TO YOUR DIVIDE-AND-CONQUER CONVEX HULL SOLVER
        t4 = time.time()

        # when passing lines to the display, pass a list of QLineF objects.  Each QLineF
        # object can be created with two QPointF objects corresponding to the endpoints
        self.showHull(polygon, RED)
        self.showText('Time Elapsed (Convex Hull): {:3.3f} sec'.format(t4 - t3))

    #This is the driver for the algorithm. It takes the series of points and divides
    #them up into 2 subproblems, of size n/2. According to the master theorem, a = 2, b=2,
    #and d =2. Therefore the time complexity of the algorithm is n^2.
    def divide_and_conquer(self, points):
        convex_hull = []
        if len(points) >= 2:
            left_hull = self.divide_and_conquer(points[:int(len(points) / 2)])
            right_hull = self.divide_and_conquer(points[int(len(points) / 2):])
        else:
            return points               #base case of one point
        maximum = self.largest(left_hull)  # rightmost point of left
        minimum = self.smallest(right_hull)
        ul, ur = self.find_upper(left_hull, right_hull, minimum, maximum)

        ll, lr = self.find_lower(left_hull, right_hull, minimum, maximum)

        convex_hull = self.make_new_hull(ul, ll, left_hull,convex_hull)
        return self.make_new_hull(lr,ur, right_hull,convex_hull)


    #This method makes a side of a new hull. It starts at the lower left
    # and looks for the upper left node to make the left side of the new hull.
    #The right side is then passed in and it starts at the lower right node and looks
    #for the upper right node. It is linear in time complexity.
    def make_new_hull(self, end, start, side, new_hull):
        encountered = False
        i = side.index(start)
        while encountered == False and i <= len(side):
            if side[i] == end: encountered = True
            new_hull.append(side[i])
            i = i+1
            if encountered == False and i >= len(side):
                i = 0                                       #if the end of the hull array has been reached and the ul/lr point has not been found, restart
        encountered = False

        return new_hull


    #This method attempts to find the upper tangent points. It starts with the points with
    #the minimum and maximum x-values and attempts to iterate counterclockwise over the
    #left hull and clockwise over the right. This works much of the time, but there are
    #always edge cases that wreak havoc on further calculations. The time complexity of this
    #method is n^2.
    def find_upper(self, left_hull, right_hull, minimum, maximum):
        upper_left = left_hull[maximum]               #maximum x value of left hull
        upper_right = right_hull[minimum]              #minimum x value of right hull
        # leftmost point of right #slope of current two points
        slope0 = self.find_slope(upper_right, upper_left)

        i = maximum-1

        while i != maximum-1:
            if i == 0 and len(left_hull) != 1:
                i = len(left_hull)-1
                # once past rightmost point, go to back of the array and
            if(len(left_hull)==1): break
            point = left_hull[i]
            i = i - 1  # go counterclockwise


            slope1 = self.find_slope(upper_right, point)
            if slope1 <= slope0:
                for point0 in right_hull:
                    slope12 = self.find_slope(upper_left, point0)
                    if slope12 < slope0: break
                    slope0 = slope12
                    upper_right = point0
            upper_left = point

        return upper_left, upper_right

    #This method attempts to find the lower tangent points of the polygon. As with find_upper,
    #it starts with the points with the minimum and maximum x-values and attempts to iterate
    #clockwise over the left hull and counterclockwise over the right hull. The complexity of this
    #method is n^2.
    def find_lower(self, left_hull, right_hull,minimum, maximum):

        lower_left = left_hull[maximum]                #maximum x_value of left hull
        lower_right = right_hull[minimum]              #minimum x_value of right hull

        slope0 = self.find_slope(lower_right, lower_left)
        #slope of current two points
        #
        #
        j = maximum - 1

        # for index in (range(0, len(left_hull))):
        while j >= 0:

            point = left_hull[j]
            slope1 = self.find_slope(lower_right, point)
            j = j - 1

            if slope1 >= slope0 or left_hull[maximum] == point:
                i = minimum - 1
                while i != minimum:
                   if i == minimum -1 and len(right_hull) != 1: i = len(right_hull)
                   elif len(right_hull) == 1: break
                   #once past leftmost point, go to back of the array and
                   i = i-1                          #go counterclockwise


                   point0 = right_hull[i]
                   slope12 = self.find_slope(lower_left, point0)
                   if slope12 > slope0: break
                   slope0 = slope12
                   lower_right = point0

                lower_left = point

        return lower_left, lower_right

    # Python3 implementation of QuickSort

    # Function to find the partition position
    def partition(self, array, low, high):

        # Choose the rightmost element as pivot
        pivot = QPointF.x(array[high])

        # Pointer for greater element
        i = low - 1

        # Traverse through all elements
        # compare each element with pivot
        for j in range(low, high):
            if QPointF.x(array[j]) <= pivot:
                # If element smaller than pivot is found
                # if array[j] <= pivot:
                # swap it with the greater element pointed by i
                i = i + 1

                # Swapping element at i with element at j
                (array[i], array[j]) = (array[j], array[i])

        # Swap the pivot element with the greater element specified by i
        (array[i + 1], array[high]) = (array[high], array[i + 1])

        # Return the position from where partition is done
        return i + 1

    # Function to perform quicksort
    def quick_sort1(self, array, low, high):
        if low < high:
            # Find pivot element such that
            # element smaller than pivot are on the left
            # element greater than pivot are on the right
            pi = self.partition(array, low, high)

            # Recursive call on the left of pivot
            self.quick_sort1(array, low, pi - 1)

            # Recursive call on the right of pivot
            self.quick_sort1(array, pi + 1, high)
        return array

    #Calculates slope. This method is constant in time complexity.
    def find_slope(self, point1, point2):
        return (QPointF.y(point1) - QPointF.y(point2)) / (QPointF.x(point1) - QPointF.x(point2))

    #finds the smallest element in the hull array. Linear in time complexity.
    def largest(self, arr):
        # Initialize maximum element
        max = arr[0]
        max_index = 0
        # Traverse array elements from second
        # and compare every element with
        # current max
        for i in range(1, len(arr)):
            if QPointF.x(arr[i]) > QPointF.x(max):
                max = arr[i]
                max_index = i
        return max_index
    #finds the point with the smallest x-value in the hull array. Linear in time
    #complexity.
    def smallest(self, arr):
        # Initialize min element
        min = arr[0]
        min_index = 0
        # Traverse array elements from second
        # and compare every element with
        # current max
        for i in range(1, len(arr)):
            if QPointF.x(arr[i]) < QPointF.x(min):
                min = arr[i]
                min_index = i
        return min_index

