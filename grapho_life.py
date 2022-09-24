# %%

from defusedxml import DTDForbidden
import graph_tool.all as gt
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
import os, sys
from gi.repository import Gtk, Gdk, GdkPixbuf, GObject, GLib

# %%

n = 80
m = 80
n_iter = 10 * 25
plt_scale = 13
density = 0.65 #era 0.65
num_states = 2
N = n * m // 10
space_scale = 2.5 #era 2.9
P = 3
Q = 2 * P
npr.RandomState(2)

OFFSCREEN = sys.argv[1] == "offscreen" if len(sys.argv) > 1 else False
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
    g.ep.weight.a = np.full(nedg, 1)

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
pos = gt.sfdp_layout(g, pos=g.vp.pos, eweight=g.ep.weight, K=0.5)
g.vp.pos = pos
g_vertices = g.get_vertices()

max_count = 2000


count = 0

if OFFSCREEN and not os.path.exists("./frames"):
    os.mkdir("./frames")

# This creates a GTK+ window with the initial graph layout
if not OFFSCREEN:
    win = gt.GraphWindow(g, pos, geometry=(plt_scale*n,plt_scale*m),
    vertex_shape="circle",
    vertex_fill_color=g.vp.color,
    vertex_size=9,
    )
else:
    win = Gtk.OffscreenWindow()
    win.set_default_size(plt_scale*n,plt_scale*m)
    win.graph = gt.GraphWidget(g, pos,
    vertex_shape="circle",
    vertex_fill_color=g.vp.color,
    vertex_size=9,
    )
    win.add(win.graph)


# %%
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
    to_remove = []
    for v in g.vertices():
        if len(g.get_all_edges(v)) == 0:
            to_remove.append(v)
        else:
            g.vp.state[v], g.vp.color[v] = rule(g, v)

    g.remove_vertex(to_remove)

    return g


# %% codecell
def update_topology(
    g,
):  # updates the graph topology, based on the dynamics i.e., the vertices states
    to_add = []
    for v in g.vertices():
        
        neighbors = g.get_all_neighbors(v)
        state_click_counter = 0
        num_neighbors = len(neighbors)

        for u in neighbors:
            if g.vp.state[v] == g.vp.state[u]:
                
                e_ = g.edge(v, u, all_edges=True)

                if g.vp.state[v] == 1:

                    state_click_counter += 1

                    for e in e_:
                        g.ep.weight[e] = min([g.ep.weight[e] + 1, Q * g.gp.age])

                elif g.vp.state[v] == 0:
                    state_click_counter -= 1
                    for e in e_:

                        if g.ep.weight[e] < 1:
                            g.remove_edge(e)

                        else:
                            if npr.rand() < 0.2:
                                g.ep.weight[e] -= 1
        # print(state_click_counter > num_neighbors//2)                                
        if state_click_counter >= 2*num_neighbors//3:
            if npr.random() < 0.025:
                to_add.append([g.vp.pos[v], v])
    # print(len(to_add), g.num_vertices()) 
    nv_1 = g.num_vertices()
    # print(nv_1)       
    for (pos, father) in to_add:
        g.add_vertex()    
        g.vp.pos[g.num_vertices() - 1] = pos + npr.random(2) * 0.25
        g.vp.state[g.num_vertices() - 1] = 1
        g.vp.color[g.num_vertices() - 1] = "white"
        g.add_edge(g.num_vertices() - 1, father)

    print(nv_1, g.num_vertices())
    g.gp.age += 1
    return g


# %%
def run_simulation():

    global count, g

    update_configuration(g, update_state)
    update_topology(g)

    gt.sfdp_layout(
        g, pos=g.vp.pos, eweight=g.ep.weight, max_iter=1, init_step=0.01, K=0.5
         )

    count += 1

    win.graph.regenerate_surface()
    win.graph.queue_draw()
    
    if OFFSCREEN:
        pixbuf = win.get_pixbuf()
        pixbuf.savev(r'./frames/graphol%06d.png' % count, 'png', [], [])

    if count >= max_count:
        sys.exit(0)

    return True


# %%

cid = GLib.idle_add(run_simulation)
win.connect("delete_event", Gtk.main_quit)
win.show_all()
Gtk.main()
    
# %%

# %%
