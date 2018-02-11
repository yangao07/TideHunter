from graph_tool.all import *

def graph_visual():
    graph = Graph()
    
    for :
        graph.add_vertex()
    graph.add_vertex(2) # start and end

    graph.edge_properties["label"] = graph.new_edge_property("string")
    graphe = graph.edge_properties["label"]

    for :
        graph.add_edge(from, to)

    for i in :
        if graph.vertex(i).in_degree() == 0:
            graph.add_edge(graph.vertex(v_n), graph.vertex(i))
        if graph.vertex(i).out_degree() == 0:
            graph.add_edge(graph.vertex(i), graph.vertex(v_n + 1))

    # graph.save('POA_graph.xml.gz')
    graph_draw(graph, vertex_text=graph.vertex_index, output_size=(10000, 10000), vertex_color=[1, 1, 1, 0], vertex_size=1,
               vcmap=matplotlib.cm.gist_heat_r, edge_pen_width=1.2,
               output='POA_graph.pdf')
