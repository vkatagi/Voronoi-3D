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

        if a == 0:
            print("A was 0")

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

def triangulate(points, superSize = 50, removeSuper = False):
    verts = []
    for p in points:
        verts.append(Vertex(p))
    
    # super polyhedron
    d = superSize
    verts.append(Vertex([d, 0, 0]))
    verts.append(Vertex([0, d, 0]))
    verts.append(Vertex([0, 0, d]))
    verts.append(Vertex([0, 0, 0]))

    superTetra = Tetrahedron(verts[-4], verts[-3], verts[-2], verts[-1])
    delaunay = set()
    delaunay.add(superTetra)

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

    if removeSuper:
        superTetras = set()
        for vert in verts[-4:]:
            superTetras = superTetras.union(vert.adj)
        
        for tetra in superTetras:
            delaunay.remove(tetra)

    return (delaunay, superTetra)

#
# Voronoi
#

class Polyface:
    def __init__(self):
        self.points = []
        self.lines = []
    
    def add(self, p):
        self.points.append(p);
        if len(self.points) > 1:
            self.lines.append([self.points[-2], self.points[-1]])

class Edge:
    def __init__(self, v0, v1):
        if v0.p[0] > v1.p[0]: # Have explicit ordering to determine uniqueness
            v0, v1 = v1, v0
        self.vtx = [v0, v1]

    def __eq__(self, other):
        return self.vtx[0] == other.vtx[0] and self.vtx[1] == other.vtx[1]

    def __hash__(self):
        return hash((self.vtx[0], self.vtx[1]))

def voronoi(delaunay, superTetra):
    edges = collections.defaultdict(Edge) # Unique delaunay edges
    
    for tetra in delaunay:
        if tetra.sharesFaceWith(superTetra):
            continue
        if m.isnan(tetra.circum_radius):
            continue
        edges.setdefault(Edge(tetra.vert[0], tetra.vert[1]), set()).add(tetra)
        edges.setdefault(Edge(tetra.vert[0], tetra.vert[2]), set()).add(tetra)
        edges.setdefault(Edge(tetra.vert[0], tetra.vert[3]), set()).add(tetra)
        edges.setdefault(Edge(tetra.vert[1], tetra.vert[3]), set()).add(tetra)
        edges.setdefault(Edge(tetra.vert[1], tetra.vert[2]), set()).add(tetra)
        edges.setdefault(Edge(tetra.vert[2], tetra.vert[3]), set()).add(tetra)

    # make a poly face for each edge
    voro = []
    for edge, tetras in edges.items():
        face = Polyface()
        
        # iterate all tetrahedrons around this edge and add their circumcenters as points
        # we check for adjacency to form the proper polygon
        current = tetras.pop()
        face.add(current.circum)

        while len(tetras) > 0:
            for tetra in tetras:
                if tetra.sharesFaceWith(current):
                    current = tetra
                    face.add(current.circum)
                    break;
            try: tetras.remove(current)
            except: break
            
        voro.append(face)



    return voro 

# imports. these are only used for plotting & drawing
import scipy as sp
import scipy.spatial
import matplotlib.pyplot as plt
import matplotlib.collections
from mpl_toolkits.mplot3d import Axes3D 
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from itertools import chain
import collections
import math as m

##
## Driver
##

def plot_delu(delu, ax):
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

def plot_voro_lines(voro, ax):
    for poly in voro:
        for a, b in poly.lines:
            ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]])

    
    
def main():

    # tetr = Tetrahedron(
    #     Vertex([0, 0, 0]), 
    #     Vertex([2, 0, 0]), 
    #     Vertex([1, 1, 0]), 
    #     Vertex([1, 0, 1])
    # )

    #points = np.array([[0.1, 0.1, 0.1], [3.2, 1.4, 0.1] , [3.1, 5, 2], [2.7, 4.1, 1], [2.9, 1, 1.1], [1, 3, 3.4], [5, 5, 1], [4.6, 2.1, 4]])
    points = np.array([[0.1, 0.1, 0.1], [4.2, 2.2, 4.2], [1, 3.0, 4.2], [2.7, 4.1, 1]])

    delu, superTetra = triangulate(points, removeSuper = False, superSize = 10)

    voro = voronoi(delu, superTetra)


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:,0], points[:,1], points[:,2])

    #plot_delu(delu, ax)

    plot_voro_lines(voro, ax)
    for poly in voro:
        for a, b in poly.lines:
            ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]])
        
    ax.set_xlim(min(points[:,0]), max(points[:,0]))
    ax.set_ylim(min(points[:,1]), max(points[:,1]))
    ax.set_zlim(min(points[:,2]), max(points[:,2]))


    plt.show()

    return



if __name__ == "__main__":
    main()