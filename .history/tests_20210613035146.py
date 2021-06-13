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
g_vertices = g.get_vertices()
offscreen = False
max_count = 2000

win = gt.GraphWindow(g, pos, geometry=(plt_scale*n,plt_scale*m),
    vertex_shape="circle",
    vertex_fill_color=g.vp.color,
    vertex_size=9,
    )
count = 0



# %%
def run_simulation():

    global count, g, pos
    
    to_remove = []
    states = []
    colors = []
    alterations = []
    r = 0
    

    for v in g.get_vertices():
        if len(g.get_all_edges(v)) == 0:
            print(r)
            r+=1
            to_remove.append(v)
        else:
            Nv = g.get_all_neighbors(v)
            sum_nbstates = sum([g.vp.state[u] for u in Nv])
            nb_size = len(Nv)
            acoef = sum_nbstates / nb_size
    
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
            
            alterations.append([v, state, color])
            states.append(state)
            colors.append(color)

    ###some checks
    #print(vertices)
    diff = len(to_remove)
    if diff >= 1:
        print(f"diff {diff}, count {count}") #this eventually results in something ≠ 0, I still haven't find out why 
    ###
        
    for (v, new_state, new_color) in alterations:
        g.vp.state[v] = new_state
        g.vp.color[v] = new_color

    g.remove_vertex(to_remove)   
    
    # for v in g.get_vertices():
    #     #print(v)
    #     g.vp.state[v] = states[v]
    #     g.vp.color[v] = colors[v]


   
    
    #g.remove_vertex(to_remove, fast=True)
    pos = gt.sfdp_layout(
        g, pos=pos, eweight=g.ep.weight, max_iter=1, init_step=0.01, K=0.5
         )

    count += 1

    win.graph.regenerate_surface()
    win.graph.queue_draw()
    
    if count >= max_count:
        sys.exit(0)

    return True


# %%

cid = GLib.idle_add(run_simulation)
win.connect("delete_event", Gtk.main_quit)
win.show_all()
Gtk.main()
    
# %%
