import numpy as np

class Vertex:
    def __init__(self, p):
        self.p = np.array(p)
        self.adj = set() # adjacent tetrahedrons

class Tetrahedron:
    def __init__(self, p0, p1, p2, p3):
        self.vert = [p0, p1, p2, p3]
        for vert in self.vert:
            vert.adj.add(self)
        self.updateCircum();

    def updateCircum(self):
        points = [vtx.p for vtx in self.vert]
        pts = points[1:] - points[0]

        (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = pts

        l1 = x1 * x1 + y1 * y1 + z1 * z1
        l2 = x2 * x2 + y2 * y2 + z2 * z2
        l3 = x3 * x3 + y3 * y3 + z3 * z3

        # Compute determinants:
        dx = +l1 * (y2 * z3 - z2 * y3) - l2 * (y1 * z3 - z1 * y3) + l3 * (y1 * z2 - z1 * y2)
        dy = +l1 * (x2 * z3 - z2 * x3) - l2 * (x1 * z3 - z1 * x3) + l3 * (x1 * z2 - z1 * x2)
        dz = +l1 * (x2 * y3 - y2 * x3) - l2 * (x1 * y3 - y1 * x3) + l3 * (x1 * y2 - y1 * x2)
        aa = +x1 * (y2 * z3 - z2 * y3) - x2 * (y1 * z3 - z1 * y3) + x3 * (y1 * z2 - z1 * y2)
        a = 2 * aa

        center = (dx / a, -dy / a, dz / a)
        self.circum = np.array([
            center[0] + points[0][0],
            center[1] + points[0][1],
            center[2] + points[0][2],
        ])
        self.circum_radius = np.linalg.norm(self.vert[0].p - self.circum)

    def sharesFaceWith(self, other):
        # TODO:
        count = 0
        for vert in self.vert:
            if vert in other.vert:
                count+=1
        return count == 2

    def vertexInCircum(self, vert):
        return np.linalg.norm(vert.p - self.circum) < self.circum_radius

    # Returns neighboring triangles set. (All triangles that share an edge with this one)
    def neighbors(self):
        # TODO: 3d
        # PERF: can be optimized
        result = set()
        for vtx in self.vert:
            for adj in vtx.adj:
                if adj is not self and self.sharesEdgeWith(adj):
                    result.add(adj)
        return result;

class Edge:
    def __init__(self, p0, p1):
        self.vert = [p0, p1]

    def __eq__(self, other):
        # Assumes p0 != p1 always.
        if self.vert[0] == other.vert[0] or self.vert[0] == other.vert[1]:
            if self.vert[1] == other.vert[0] or self.vert[1] == other.vert[1]:
                return True
        return False



def findPolygon(badTriangles):
    # PERF: can be optimized a lot
    all_edges = []
    for tri in badTriangles:
        all_edges.append(Edge(tri.vert[0], tri.vert[1]))
        all_edges.append(Edge(tri.vert[1], tri.vert[2]))
        all_edges.append(Edge(tri.vert[2], tri.vert[0]))
    
    edges = []
    for edge in all_edges:
        found = 0
        for tested in all_edges:
            if tested == edge:
                found += 1
                if found > 1:
                    break
        if found == 1:
            edges.append(edge)

    return edges


def triangulate(points):
    verts = []
    for p in points:
        verts.append(Vertex(p))
    
    # super triangle
    verts.append(Vertex([1000, 0]))
    verts.append(Vertex([0, 1000]))
    verts.append(Vertex([0, 0]))

    triangles = set()
    triangles.add(Triangle(verts[-3], verts[-2], verts[-1]))

    for vtx in verts:
        # Find bad triangles, (all the triangles that vtx is inside of the circumcircle)
        bad = set(filter(lambda tri: (tri.vertexInCircum(vtx)), triangles))

        # Find polygonal hole from bad triangles
        polygon = findPolygon(bad)

        # Remove bad triangles from triangulation
        for tri in bad:
            for triVtx in tri.vert:
                triVtx.adj.remove(tri)
            triangles.remove(tri)

        for edge in polygon:
            # assumes all edges are valid (non points)
            triangles.add(Triangle(vtx, edge.vert[0], edge.vert[1]))
    
    return triangles
        
def voronoi(triangulation):
    edges = []
    for tri in triangulation:
        for neighbor in tri.neighbors():
            edges.append(Edge(tri.circum, neighbor.circum))
    return edges 



##
## Driver
##

# imported later. these are only used for plotting & drawing
import scipy as sp
import scipy.spatial
import matplotlib.pyplot as plt
import matplotlib.collections
from mpl_toolkits.mplot3d import Axes3D 


def main():
    tetr = Tetrahedron(
        Vertex([0, 0, 0]), 
        Vertex([2, 0, 0]), 
        Vertex([1, 1, 0]), 
        Vertex([1, 0, 1])
    )


    return
    points = np.array([[0, 0, 0], [3.2, 1.4, 0.1] , [3.1, 5, 2], [2.7, 4.1, 1], [2.9, 1, 1.1], [1, 3, 3.4], [5, 5, 1], [4.6, 2.1, 4]])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    print(points[:,0]);

    ax.scatter(points[:,0], points[:,1], points[:,2])
    plt.show()


    # triangles = triangulate(points)
    return

    # voro_edges = voronoi(triangles)
    
    # lines = []
    # for edge in voro_edges:
    #     lines.append(edge.vert)
    # lc = matplotlib.collections.LineCollection(lines, linewidths=2)
    # fig, ax = plt.subplots()
    # ax.add_collection(lc)

    # plt.xlim([-0.5,5.5])
    # plt.ylim([-0.5,5.5])

    # vor = sp.spatial.Voronoi(points)
    # fig = sp.spatial.voronoi_plot_2d(vor)

    # plt.show()
    return


if __name__ == "__main__":
    main()