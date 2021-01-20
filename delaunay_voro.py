import sys
if sys.version_info[0] < 3:
    raise Exception("delaunay_voro.py requires Python 3.")

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
        self.updateCircum()

    def updateCircum(self):
        ax = self.vert[0].p[0]; ay = self.vert[0].p[1]
        bx = self.vert[1].p[0]; by = self.vert[1].p[1]
        cx = self.vert[2].p[0]; cy = self.vert[2].p[1]

        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
        uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
        self.circum = [ux, uy]
        self.circum_radius = np.linalg.norm(self.vert[0].p - self.circum)

    def sharesEdgeWith(self, other):
        count = 0
        for vert in self.vert:
            if vert in other.vert:
                count+=1
        return count == 2

    def vertexInCircum(self, vert):
        return np.linalg.norm(vert.p - self.circum) < self.circum_radius

    # Returns neighboring triangles set. (All triangles that share an edge with this one)
    def neighbors(self):
        # PERF: can be optimized
        result = set()
        for vtx in self.vert:
            for adj in vtx.adj:
                if adj is not self and self.sharesEdgeWith(adj):
                    result.add(adj)
        return result

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
    verts.append(Vertex([-500, -500]))
    verts.append(Vertex([1000, -500]))
    verts.append(Vertex([-500, 1000]))

    triangles = set()
    triangles.add(Triangle(verts[-3], verts[-2], verts[-1]))

    for vtx in verts[:-3]:
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

# these are only used for plotting & drawing
import scipy as sp
import scipy.spatial
import matplotlib.pyplot as plt
import matplotlib.collections

def get_cmap(n, name='hsv'):
    return plt.cm.get_cmap(name, n)

def plot_tris(triangles):
    cmap = get_cmap(len(triangles))
    
    for i, tri in enumerate(triangles):
        pts = [tri.vert[0].p, tri.vert[1].p, tri.vert[2].p]
        t = plt.Polygon(pts, True, color=cmap(i))
        plt.gca().add_patch(t)

def main():
    # Points expected to be in [0, 1) range
    NUM_POINTS = 10
    points = np.random.rand(NUM_POINTS, 2)
    

    triangles = triangulate(points)

    voro_edges = voronoi(triangles)
    
    lines = []
    for edge in voro_edges:
        lines.append(edge.vert)
    lc = matplotlib.collections.LineCollection(lines, linewidths=2)
    fig, ax = plt.subplots()
    ax.add_collection(lc)

    ax.scatter(*zip(*points))

    plt.xlim([-0.1,1.1])
    plt.ylim([-0.1,1.1])
    fig.canvas.set_window_title('Our Voronoi Implementation (using delaunay O(N log N))')



    vor = sp.spatial.Voronoi(points)
    # NOTE: this will not exactly match our plot limits so the result will not look exactly the same
    fig = sp.spatial.voronoi_plot_2d(vor)
    fig.canvas.set_window_title('sp.spatial.Voronoi() for reference')

    plt.show()
    return


if __name__ == "__main__":
    main()