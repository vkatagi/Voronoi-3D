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

    def vertexInCircum(self, vert):
        return np.linalg.norm(vert.p - self.circum) < self.circum_radius
    
    def sharesFaceWith(self, other):
        return len(set(self.vert).intersection(other.vert)) >= 3

    # Returns neighboring triangles set. (All triangles that share an edge with this one)
    # def neighbors(self):
    #     # TODO: 3d
    #     # PERF: can be optimized
    #     result = set()
    #     for vtx in self.vert:
    #         for adj in vtx.adj:
    #             if adj is not self and self.sharesEdgeWith(adj):
    #                 result.add(adj)
    #     return result;

class Face:
    def __init__(self, p0, p1, p2):
        self.vert = [p0, p1, p2]

    def __eq__(self, other):
        return len(set(self.vert).intersection(other.vert)) == 3
        
def findPolyhedron(badTetras):
    # PERF: can be optimized a lot
    all_faces = []
    for poly in badTetras:
        all_faces.append(Face(poly.vert[0], poly.vert[1], poly.vert[2]))
        all_faces.append(Face(poly.vert[1], poly.vert[2], poly.vert[3]))
        all_faces.append(Face(poly.vert[2], poly.vert[3], poly.vert[0]))
        all_faces.append(Face(poly.vert[3], poly.vert[0], poly.vert[1]))
    
    faces = []
    for face in all_faces:
        found = 0
        for tested in all_faces:
            if tested == face:
                found += 1
                if found > 1:
                    break
        if found == 1:
            faces.append(face)

    return faces

def triangulate(points):
    verts = []
    for p in points:
        verts.append(Vertex(p))
    
    # super polyhedron
    d = 100
    verts.append(Vertex([d, 0, 0]))
    verts.append(Vertex([0, d, 0]))
    verts.append(Vertex([0, 0, d]))
    verts.append(Vertex([0, 0, 0]))

    delaunay = set()
    delaunay.add(Tetrahedron(verts[-4], verts[-3], verts[-2], verts[-1]))

    for vtx in verts:
        # Find bad tetrahedrons
        bad = set(filter(lambda tetr: (tetr.vertexInCircum(vtx)), delaunay))

        # Find polygonal hole from bad tetrahedrons
        polyhedron = findPolyhedron(bad)

        # Remove bad tetrahedrons from triangulation
        for tetra in bad:
            for tetraVtx in tetra.vert:
                tetraVtx.adj.remove(tetra)
            delaunay.remove(tetra)

        for face in polyhedron:
            delaunay.add(Tetrahedron(face.vert[0], face.vert[1], face.vert[2], vtx))

    superTetras = set()
    for vert in verts[-4:]:
        superTetras = superTetras.union(vert.adj)
    
    for tetra in superTetras:
        delaunay.remove(tetra)

    return delaunay
        
def voronoi(triangulation):
    edges = []
    for tri in triangulation:
        for neighbor in tri.neighbors():
            edges.append(Edge(tri.circum, neighbor.circum))
    return edges 

# imports. these are only used for plotting & drawing
import scipy as sp
import scipy.spatial
import matplotlib.pyplot as plt
import matplotlib.collections
from mpl_toolkits.mplot3d import Axes3D 
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from itertools import chain

##
## Driver
##

def main():

    # tetr = Tetrahedron(
    #     Vertex([0, 0, 0]), 
    #     Vertex([2, 0, 0]), 
    #     Vertex([1, 1, 0]), 
    #     Vertex([1, 0, 1])
    # )

    points = np.array([[0.1, 0.1, 0.1], [3.2, 1.4, 0.1] , [3.1, 5, 2], [2.7, 4.1, 1], [2.9, 1, 1.1], [1, 3, 3.4], [5, 5, 1], [4.6, 2.1, 4]])

    delu = triangulate(points)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:,0], points[:,1], points[:,2])

    lines = []
    for tetra in delu:
        lines.append((tetra.vert[0].p, tetra.vert[1].p))
        lines.append((tetra.vert[0].p, tetra.vert[2].p))
        lines.append((tetra.vert[0].p, tetra.vert[3].p))
        lines.append((tetra.vert[1].p, tetra.vert[3].p))
        lines.append((tetra.vert[1].p, tetra.vert[2].p))
        lines.append((tetra.vert[2].p, tetra.vert[3].p))
        
    for a, b in lines:
        ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]])

    plt.show()

    return

    



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