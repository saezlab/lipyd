class Rectangle(object):
    def __init__ (
            self,
            x1 = .0,    # lower left corner;
            y1 = .0, 
            x2 = .0,    # upper right corner;
            y2 = .0,
            x = .0,     # original value;
            y = .0,
            **kwargs,
        ):
            self.x1 = x1 + 0.0 # lower left corner
            self.y1 = y1 + 0.0
            self.x2 = x2 + 0.0 # upper right corner
            self.y2 = y2 + 0.0
            self.x = x + 0.0
            self.y = y + 0.0

    def set_points(self,
                    x1 = .0,
                    y1 = .0, 
                    x2 = .0,
                    y2 = .0,
                    x = .0,
                    y = .0,
        ):
            self.x1 = x1 + 0.0 # lower left corner
            self.y1 = y1 + 0.0
            self.x2 = x2 + 0.0 # upper right corner
            self.y2 = y2 + 0.0
            self.x = x + 0.0
            self.y = y + 0.0

    def get_points(self):
        return (self.x1, self.y1, self.x2, self.y2, self.x, self.y)

    def get_center(self):   # return (x,y) of center of rectangle;
        return (self.x1+abs(self.x2 - self.x1)/2, self.y1+abs(self.y2 - self.y1)/2)

    def is_rectangle_intersect_rectangle(self, r): # r is rectangle;
        r_in_self = False
        #print(80*"-")
        #print("intersect: r:",r.x1,r.y1,r.x2,r.y2)
        #print("intersect: s:",self.x1,self.y1,self.x2,self.y2)
        
        if r.x1 <= self.x2 and r.x2 >= self.x1 and \
            r.y1 <= self.y2 and r.y2 >= self.y1 :
                r_in_self = True
        
        return r_in_self


class Layout_shape(object):
    pad_x = 5
    pad_y = 5
    cellule_x = 10
    cellule_y = 10
    number_steps = 50

    def __init__(
            self,
            plot_x_min = 0,
            plot_x_max = 400,
            plot_y_min = 0,
            plot_y_max = 100,
            annot_list = [], #structure element of list: (mz, intens, annot)
            **kwargs,
        ):
        self.plot_x_min = plot_x_min
        self.plot_x_max = plot_x_max
        self.plot_y_min = plot_y_min
        self.plot_y_max = plot_y_max
        self.shape_list = []
        self.cellule_list = []
        self.trial_cellule = Rectangle()
        
        for e in annot_list:      # e[0] is x; e[1] is y;
            self.shape_list.append( self.create_rectangle(e[0], 0, e[0], e[1], e[0], e[1]) )

    def create_rectangle(self, x1, y1, x2, y2, x, y):
        '''
        x_1 = x1 - Layout_shape.pad_x \
                if (x1 - Layout_shape.pad_x) > self.plot_x_min else self.plot_x_min

        y_1 = y1 - Layout_shape.pad_y \
                if (y1 - Layout_shape.pad_y) > self.plot_y_min else self.plot_y_min

        x_2 = x2 + Layout_shape.pad_x \
                if (x2 + Layout_shape.pad_x) < self.plot_x_max else self.plot_x_max

        y_2 = y2 + Layout_shape.pad_y \
                if (y2 + Layout_shape.pad_y) < self.plot_y_max else self.plot_y_max
        '''

        return Rectangle(x1, y1, x2, y2, x, y)

    def create_cellule(self, x, y):
        x_1 = x - Layout_shape.cellule_x \
                if (x - Layout_shape.cellule_x) > self.plot_x_min else self.plot_x_min

        y_1 = y - Layout_shape.cellule_y \
                if (y - Layout_shape.cellule_y) > self.plot_y_min else self.plot_y_min

        x_2 = x + Layout_shape.cellule_x \
                if (x + Layout_shape.cellule_x) < self.plot_x_max else self.plot_x_max

        y_2 = y + Layout_shape.cellule_y \
                if (y + Layout_shape.cellule_y) < self.plot_y_max else self.plot_y_max

        return Rectangle(x_1, y_1, x_2, y_2, x, y)


    def add_cellule(self, x, y):
        self.cellule_list.append( self.create_cellule(x, y) )

    def is_cellule_free(self, cellule): # cellule is object of Rectangle class;
        """
        True  - the place is free;
        False - the place is allocate;
        """
        for c in self.cellule_list:
            if c.is_rectangle_intersect_rectangle(cellule):
                return False

        return True


    def is_shape_free(self, rectangle):
        """
        True  - the place is free;
        False - the place is allocate;
        """
        for c in self.shape_list:
            if c.is_rectangle_intersect_rectangle(rectangle):
                return False

        return True


    def get_free_place(self, mz):
        for e in self.shape_list:
            x1, y1, x2, y2, x, y = e.get_points()
            #print("get_free_place: ", x1, y1, x2, y2, x, y)
            if mz != x:
                continue

            shift_y = y

            for s in range(Layout_shape.number_steps):
                if shift_y >= self.plot_y_max:
                    return (x, y)

                self.trial_cellule = self.create_cellule(x, shift_y)
                if self.is_shape_free(self.trial_cellule) and \
                    self.is_cellule_free(self.trial_cellule):
                        self.add_cellule(x, y)
                        return self.trial_cellule.get_center()

                shift_y += Layout_shape.cellule_y

            return (0, 0)

