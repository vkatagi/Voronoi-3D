import numpy as np

class Vertex:
    def __init__(self, p):
        self.p = np.array(p)
        self.adj = set() # adjacent triangles

class Triangle:
    def __init__(self, p0, p1, p2):
        self.vert = [p0, p1, p2]
        for vert in self.vert:
            vert.adj.add(self)
        self.updateCircum();

    def updateCircum(self):
        ax = self.vert[0].p[0]; ay = self.vert[0].p[1]
        bx = self.vert[1].p[0]; by = self.vert[1].p[1]
        cx = self.vert[2].p[0]; cy = self.vert[2].p[1]

        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
        uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
        self.circum = [ux, uy]
        self.circum_radius = np.linalg.norm(self.vert[0] - self.circum)

    def sharesEdgeWith(self, other):
        count = 0
        for vert in self.vert:
            if vert in other.vert:
                count+=1
        return count == 2

    def vertexInCircum(self, vert):
        return np.linalg.norm(vert.p - self.circum) < self.circum_radius

class Edge:
    def __init__(self, p0, p1):
        self.vert = [p0, p1]

    def __eq__(self, other):
        # Assumes p0 != p1 always.
        if self.vert[0] == other.vert[0] or self.vert[0] == other.vert[1]:
            if self.vert[1] == other.vert[0] or self.vert[1] == other.vert[1]:
                return True
        return False


class Voronoi:
    def __init__(self, points):
        self.points = points
        



##
## Driver
##

# imported later. these are only used for plotting & drawing
import scipy as sp
import scipy.spatial
# import pylab
import matplotlib.pyplot as plt
import matplotlib.collections

def test():
    tri = Triangle([-1,0], [0,1], [1, 0])


    circle = plt.Circle(tri.circumcenter(), tri.circumradius())

    fig, ax = plt.subplots()
    
    ax.add_artist(circle)
    ax.set_aspect('equal','box')
    plt.show()


def main():
    test()
    return

    points = [[0, 0], [3.2, 1.4], [3.1, 5], [2.7, 4.1], [2.9, 1], [1, 3], [5, 5], [4.6, 2.1]]
    v = Voronoi(points)
    
    lines = []
    
    for seg in v.segments:
        lines.append([(seg.start.x, seg.start.y), (seg.end.x, seg.end.y)])

    lines.append([(0, 0), (1, 1)])
    
    lc = matplotlib.collections.LineCollection(lines, linewidths=2)
    fig, ax = plt.subplots()
    ax.add_collection(lc)
    x,y = np.array(points).T
    plt.xlim([-0.5,5.5])
    plt.ylim([-0.5,5.5])
    plt.scatter(x, y)



    vor = sp.spatial.Voronoi(points)
    fig = sp.spatial.voronoi_plot_2d(vor)

    plt.show()


if __name__ == "__main__":
    main()