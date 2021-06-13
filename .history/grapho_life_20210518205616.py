# %% codecell
import graph_tool.all as gt
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
import os
from gi.repository import Gtk, Gdk, GdkPixbuf, GObject, GLib
import argparse

# %%

parser = argparse.ArgumentParser(
    description="Open an interactive GTK window running a Graph of Life."
)
parser.add_argument(
    "is_offscreen",
    action="store_true",
    help="if True, will dump the frames to the disk instead of running an interactive GTK window",
)
parser.add_argument(
    "niter",
    metavar="n_iter",
    action="store_const",
    const = 60,
    default = 60,
    help="number of iterations to be performed. If zero, runs until convergence or a predefined maximum number of iterations",
)

parser.add_argument(
    "-dims",
    "--dimensions",
    action="store_const",
    default=[50, 50],
    nargs = 2,
    help="-dims n m where n, m are the size of the underlying space used for the generation of the base random geometric graph. random ",
)

parser.add_argument(
    "-dens",
    "--density",
    action="store_const",
    const = 0.65
    default=0.65,
    help="parameter of the binomial distribution whence the states come; thus the higher the value, the denser the number of 'live' states",
)

# %%

args = parser.parse_args()
print(args.accumulate(args.integers))

n = 75
m = 75
n_iter = 10 * 25
plt_scale = 13
density = 0.65
num_states = 2
N = n * m // 10
print(N)
space_scale = 2.9
P = 3
Q = 2 * P
npr.RandomState(2)
# %% codecell
def gen_ints(N, p, n=1):
    # n, p = 1, .5  number of trials, probability of each trial
    s = np.random.binomial(n, p, N)
    return s


# %% codecell
def make_rpGA(dim, num_states, density=0.5):
    points = npr.random((N, 2)) * space_scale
    g = gt.Graph(directed=False)
    g.vp.pos = g.new_vertex_property("vector<double>")
    g, g.vp.pos = gt.geometric_graph(points, 10 * space_scale / (n + m))
    g.gp.age = g.new_graph_property("int")

    g.vp.state = g.new_vertex_property("int")
    g.vp.color = g.new_vertex_property("string")
    g.ep.weight = g.new_edge_property("int")
    nedg = g.num_edges()
    g.ep.weight.a = np.array([1 for i in range(nedg)])

    states_vals = gen_ints(dim, density, n=1)
    # print(sum(states_vals))
    # print(g.list_properties())
    for i, v in enumerate(g.vertices()):
        g.gp.age = 0
        g.vp.state[v] = states_vals[i]
        if g.vp.state[v] == 0:
            g.vp.color[v] = "black"
        else:
            g.vp.color[v] = "white"
    return g


# %% codecell
g = make_rpGA(N, num_states, density=density)
# %% codecell
def update_state(
    g, v
):  # updates the state of a vertex with rules analogous to the original GoL

    Nv = g.get_all_neighbors(v)
    sum_nbstates = sum([g.vp.state[u] for u in Nv])
    nb_size = len(Nv)
    acoef = sum_nbstates / nb_size
    # print(acoef)
    if g.vp.state[v] == 1:
        if 2 / nb_size <= acoef <= 3 / nb_size:
            g.vp.state[v] = 1
            g.vp.color[v] = "white"
        else:
            g.vp.state[v] = 0
            g.vp.color[v] = "black"

    elif g.vp.state[v] == 0:
        # if (acoef == 3/nb_size):
        if 2 / nb_size <= acoef <= 4 / nb_size:
            g.vp.state[v] = 1
            g.vp.color[v] = "white"
        else:
            g.vp.state[v] = 0
            g.vp.color[v] = "black"
    return g.vp.state[v], g.vp.color[v]


# %% codecell
def update_configuration(g, rule):  # applies the update_state function to each vertex
    G = g.copy()
    to_remove = []
    for v in g.vertices():
        if len(g.get_all_edges(v)) == 0:
            to_remove.append(v)
        else:
            G.vp.state[v], G.vp.color[v] = rule(G, v)

    G.remove_vertex(to_remove)
    return G


# %% codecell
def update_topology(
    g,
):  # updates the graph topology, based on the dynamics i.e., the vertices states
    G = g.copy()
    for v in g.vertices():
        neighbors = g.get_all_neighbors(v)
        for u in neighbors:
            if g.vp.state[v] == g.vp.state[u]:

                e_ = G.edge(v, u, all_edges=True)

                if g.vp.state[v] == 1:

                    for e in e_:
                        G.ep.weight[e] = min([G.ep.weight[e] + 1, Q * G.gp.age])

                elif g.vp.state[v] == 0:
                    for e in e_:

                        if G.ep.weight[e] < 1:
                            G.remove_edge(e)

                        else:
                            if npr.rand() < 0.2:
                                G.ep.weight[e] -= 1

    G.gp.age += 1
    return G


# %% codecell
def make_gif(g, rule, n_iter=n_iter, outname="def", draw=True):
    fname_prefix = os.path.join("./frames", outname)
    if draw == True:
        gt.graph_draw(
            g,
            pos=g.vp.pos,
            vertex_shape="circle",
            vertex_fill_color=g.vp.color,
            vertex_size=9,
            output_size=(plt_scale * n, plt_scale * m),
            output=fname_prefix + "_draw0000.png",
        )

    Gs = [g]

    for i in range(1, n_iter + 1):
        G = update_configuration(Gs[-1], rule=rule)
        # print(len(G.vp.state.a), len(Gs[-1].vp.state.a))
        msize = min([len(G.vp.state.a), len(Gs[-1].vp.state.a)])

        if np.all(G.vp.state.a[-msize:-1] == Gs[-1].vp.state.a[-msize:-1]):
            print("early convergence")
            return Gs
        else:
            G = update_topology(G)
            npos = gt.sfdp_layout(
                G, pos=G.vp.pos, eweight=G.ep.weight, max_iter=1, init_step=0.01, K=0.5
            )
            G.vp.pos = npos

        if draw == True:
            fname_prefix = os.path.join("./frames", outname)
            gt.graph_draw(
                G,
                pos=G.vp.pos,
                vertex_shape="circle",
                vertex_fill_color=G.vp.color,
                vertex_size=9,
                output_size=(plt_scale * n, plt_scale * m),
                output=fname_prefix + "_draw%04d.png" % i,
            )
        fname_prefix + "_draw%04d.png"
        Gs.append(G)

    return Gs


# %% codecell
rule = update_state
# %% codecell
g = make_rpGA(N, num_states, density=density)
g.vp.pos = gt.sfdp_layout(g, pos=g.vp.pos, eweight=g.ep.weight, K=0.5)
gt.graph_draw(
    g,
    pos=g.vp.pos,
    vertex_shape="circle",
    vertex_fill_color=g.vp.color,
    vertex_size=9,
    output_size=(plt_scale * n, plt_scale * m),
)
# %% codecell

offscreen = sys.argv[1] == "offscreen" if len(sys.argv) > 1 else False
max_count = 5000
if offscreen and not os.path.exists("./frames"):
    os.mkdir("./frames")

# This creates a GTK+ window with the initial graph layout
if not offscreen:
    win = GraphWindow(g, pos, geometry=(500, 400))
else:
    win = Gtk.OffscreenWindow()
    win.set_default_size(500, 400)
    win.graph = GraphWidget(g, pos)
    win.add(win.graph)
# %% codecell

# %% codecell
