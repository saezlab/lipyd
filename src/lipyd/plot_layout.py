import math

class Rectangle(object):
    
    def __init__ (
            self,
            x1 = .0,    # lower left corner;
            y1 = .0, 
            x2 = .0,    # upper right corner;
            y2 = .0,
            **kwargs,
        ):
            self.x1 = x1 + 0.0 # lower left corner
            self.y1 = y1 + 0.0
            self.x2 = x2 + 0.0 # upper right corner
            self.y2 = y2 + 0.0
            self.orig = (self.x1, self.y1, self.x2, self.y2) 

    def set_values(self, t):
            self.x1 = t[0] + 0.0 # lower left corner
            self.y1 = t[1] + 0.0
            self.x2 = t[2] + 0.0 # upper right corner
            self.y2 = t[3] + 0.0

    def set_orig(self, t):
        self.orig = (t[0], t[1], t[2], t[3])
            
    def get_orig(self):
        return self.orig

    def get_values(self):
        return (self.x1, self.y1, self.x2, self.y2)

    def is_intersect_rectangles(self, r): # r is rectangle;
        r_in_self = False
        
        if r.x1 <= self.x2 and r.x2 >= self.x1 and \
            r.y1 <= self.y2 and r.y2 >= self.y1 :
                r_in_self = True
        
        return r_in_self


class Layout_shape(object):
    def __init__(
            self,
            plot_x_min = 60, #60
            plot_x_max = 800, #800
            plot_y_min = 20, #20
            plot_y_max = 100, #100
            annot_list = [], #structure element of list: (mz, intens, annot)
            **kwargs,
        ):
        self.plot_x_min = plot_x_min
        self.plot_x_max = plot_x_max
        self.plot_y_min = plot_y_min
        self.plot_y_max = plot_y_max
        
        self.cell_x = 40 #50 
        self.cell_y = 3 #7
        self.step_x = 1 #1
        self.step_y = 1 #1
        self.search_x_min = 50 #100
        self.search_x_max = 300 #500
        self.search_y_min = 20 #20
        self.search_y_max = 80 #80
        self.search_area = (0,0,0,0) # area of search;
        self.no_next_point = (0, 0)
        
        self.shape_list = []    # [(x1, y1, x2, y2), ...]
        
        for e in annot_list:      # e[0] is x; e[1] is y;
            self.shape_list.append( (e[0], 0, e[0], e[1]) )

    def set_search_area(self, x, y): # set size of search area for peak;
        x_min = x - self.search_x_min \
            if x - self.search_x_min >= self.plot_x_min else self.plot_x_min
        x_max = x + self.search_x_max \
            if x + self.search_x_max <= self.plot_x_max else self.plot_x_max
        y_min = y + self.search_y_min \
            if y + self.search_y_min <= .8 * self.plot_y_max else .8 * y
        y_max = y + self.search_y_max \
            if y + self.search_y_max <= .95 * self.plot_y_max else .95 * y
        
        self.search_area = (x_min, x_max, y_min, y_max)

    def is_point_in_search_area(self, x, y):
        return True if self.search_area[0] <= x and \
                        self.search_area[1] >= x and \
                        self.search_area[2] <= y and \
                        self.search_area[3] >= y \
                        else False

    def add_cellule(self, x1, y1, x2, y2 ):
        self.shape_list.append( (x1, y1, x2, y2) )

    def make_cell(self, x1, y1, x2, y2):
        xa = x1 - self.cell_x \
                if (x1 - self.cell_x) > self.plot_x_min else self.plot_x_min
        ya = y1 - self.cell_y \
                if (y1 - self.cell_y) > self.plot_y_min else self.plot_y_min
        xb = x2 + self.cell_x \
                if (x2 + self.cell_x) < self.plot_x_max else self.plot_x_max
        yb = y2 + self.cell_y \
                if (y2 + self.cell_y) < self.plot_y_max else self.plot_y_max
        return (xa, ya, xb, yb)

    def make_cell2(self, t): # t is tuple;
        return self.make_cell( t[0], t[1], t[2], t[3] )

    def is_cellule_free(self, rc):
        p = Rectangle()
        p.set_values( self.make_cell2( rc ) )
        q = Rectangle()

        for c in self.shape_list:
            q.set_values( self.make_cell2( c ) )
            if q.is_intersect_rectangles( p ):
                return False
        return True

    def get_next_point(self, x1, y1 ):
        x = x1
        y = y1
        x += self.step_x
        if not self.is_point_in_search_area(x, y):
            x = self.search_area[0]
            y += self.step_y
            if not self.is_point_in_search_area(x, y):
                x, y = self.no_next_point   # no next point
        return (x, y)

    def get_free_place(self, x, y):
        self.set_search_area(x, y)
        
        q = self.search_area[0], self.search_area[2]    # x of search area, y of search area
        cell = Rectangle()  # current cellule;
        cell.set_values( self.make_cell(q[0], q[1], q[0], q[1]) )
        cell.set_orig( (q[0], q[1], q[0], q[1]) )
        
        while True:
            orig = cell.get_orig()
            if self.is_cellule_free( orig ):
                self.add_cellule( orig[0], orig[1], orig[2], orig[3] )
                return orig[0], orig[1]

            q = self.get_next_point(q[0], q[1]) # get next point;
            if q == self.no_next_point:
                self.add_cellule( x, y, x, y )
                return x, y   #no free points;

            cell.set_orig( (q[0], q[1], q[0], q[1]) )
            cell.set_values( self.make_cell(q[0],q[1],q[0],q[1]) )
        # end of get_free_place;

