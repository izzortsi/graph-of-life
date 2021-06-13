# %% codecell
from defusedxml import DTDForbidden
import graph_tool.all as gt
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
import os, sys
from gi.repository import Gtk, Gdk, GdkPixbuf, GObject, GLib


# %%

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


# %%
Nv = g.get_all_neighbors(1)

# %%

# %%

Nv
# %% codecell
def update_state(g, v):  # updates the state of a vertex with rules analogous to the original GoL

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
    states = []
    colors = []
    for v in g.iter_vertices():
        if len(g.get_all_edges(v)) == 0:
            to_remove.append(v)
        else:
            state, color = rule(g, v)
            states.append(state)
            colors.append(color)
    
    states = np.array(states)
    colors = np.array(colors)

    g.vp.state.a = states
    g.vp.color.a = colors
    g.remove_vertex(to_remove)

    return g


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
rule = update_state
# %% codecell
g = make_rpGA(N, num_states, density=density)
pos = gt.sfdp_layout(g, pos=g.vp.pos, eweight=g.ep.weight, K=0.5)

# %% codecell

#offscreen = sys.argv[1] == "offscreen" if len(sys.argv) > 1 else False

max_count = 250

win = gt.GraphWindow(g, pos, geometry=(500, 400),
    vertex_shape="circle",
    vertex_fill_color=g.vp.color,
    vertex_size=9,
    output_size=(plt_scale * n, plt_scale * m),
    )
count = 0
##
def run_simulation():

    global count
    
    to_remove = []
    states = []
    colors = []
    for v in g.iter_vertices():
        if len(g.get_all_edges(v)) == 0:
            to_remove.append(v)
        else:
            Nv = g.get_all_neighbors(v)
            sum_nbstates = sum([g.vp.state[u] for u in Nv])
            nb_size = len(Nv)
            acoef = sum_nbstates / nb_size
    # print(acoef)
        if g.vp.state[v] == 1:
            if 2 / nb_size <= acoef <= 3 / nb_size:
                state = 1
                color = "white"
            else:
                state = 0
                color = "black"

        elif g.vp.state[v] == 0:
            # if (acoef == 3/nb_size):
            if 2 / nb_size <= acoef <= 4 / nb_size:
                state = 1
                color = "white"
            else:
                state = 0
                color = "black"

            
            states.append(state)
            colors.append(color)
    
    states = np.array(states)
    colors = np.array(colors)

    g.vp.state.a = states
    g.vp.color.a = colors
    g.remove_vertex(to_remove)

    # print(len(G.vp.state.a), len(Gs[-1].vp.state.a))
    #g = update_topology(g)
    gt.sfdp_layout(
        g, pos=g.vp.pos, eweight=g.ep.weight, max_iter=1, init_step=0.01, K=0.5
        )

    count += 1

    win.graph.regenerate_surface()
    win.graph.queue_draw()

    if count >= max_count:
        sys.exit(0)

    return True

##

cid = GLib.idle_add(run_simulation)

# We will give the user the ability to stop the program by closing the window.
win.connect("delete_event", Gtk.main_quit)

# Actually show the window, and start the main loop.
win.show_all()
Gtk.main()
##
##

# %%

# %%
