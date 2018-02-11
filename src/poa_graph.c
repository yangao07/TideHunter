#include <stdio.h>
#include <stdlib.h>
#include "poa.h"
#include "poa_graph.h"
#include "poa_align.h"
#include "utils.h"

poa_node_t *poa_init_node(int n) {
    poa_node_t *node = (poa_node_t*)_err_calloc(n, sizeof(poa_node_t));
    return node;
}

void poa_set_graph_node(poa_graph_t *graph, int node_i) {
    graph->node[node_i].node_id = -1;
    graph->node[node_i].cumul_len = 0;
    graph->node[node_i].in_edge_n = 0;
    graph->node[node_i].in_edge_m = 0;
    graph->node[node_i].out_edge_n = 0;
    graph->node[node_i].out_edge_m = 0;
    graph->node[node_i].aligned_node_n = 0;
    graph->node[node_i].aligned_node_m = 0;

}

void poa_free_node(poa_node_t *node, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (node[i].in_edge_m > 0) free(node[i].in_id);
        if (node[i].out_edge_m > 0) free(node[i].out_id);
        if (node[i].aligned_node_m > 0) free(node[i].aligned_node_id);
    }
    free(node);
}

poa_edge_t *poa_init_edge(int n) {
    poa_edge_t *edge = (poa_edge_t*)_err_calloc(n, sizeof(poa_edge_t));
    return edge;
}

void poa_free_edge(poa_edge_t *edge, int n) {
    free(edge);
}

// 0: in_edge, 1: out_edge
poa_graph_t *poa_realloc_graph_edge(poa_graph_t *graph, int io, int id) {
    if (io == 0) {
        if (graph->node[id].in_edge_m <= 0) {
            graph->node[id].in_edge_m = 1;
            graph->node[id].in_id = (int*)_err_malloc(sizeof(int));
        }
        if (graph->node[id].in_edge_n == graph->node[id].in_edge_m) {
            graph->node[id].in_edge_m <<= 1;
            graph->node[id].in_id = (int*)_err_realloc(graph->node[id].in_id, graph->node[id].in_edge_m * sizeof(int));
        }
    } else {
        if (graph->node[id].out_edge_m <= 0) {
            graph->node[id].out_edge_m = 1;
            graph->node[id].out_id = (int*)_err_malloc(sizeof(int));
        }
        if (graph->node[id].out_edge_n == graph->node[id].out_edge_m) {
            graph->node[id].out_edge_m <<= 1;
            graph->node[id].out_id = (int*)_err_realloc(graph->node[id].out_id, graph->node[id].out_edge_m * sizeof(int));
        }
    }
    return graph;
}

poa_graph_t *poa_realloc_graph_node(poa_graph_t *graph) {
    if (graph->node_m <= 0) {
        graph->node_m = 1;
        graph->node = (poa_node_t*)_err_calloc(1, sizeof(poa_node_t));
    }
    if (graph->node_n == graph->node_m) {
        int i;
        graph->node_m <<= 1;
        graph->node = (poa_node_t*)_err_realloc(graph->node, graph->node_m * sizeof(poa_node_t));
        for (i = graph->node_m >> 1; i < graph->node_m; ++i) {
            poa_set_graph_node(graph, i);
        }
    }
    return graph;
}

poa_graph_t *poa_init_graph(int n) {
    poa_graph_t *graph = (poa_graph_t*)_err_calloc(n, sizeof(poa_graph_t));
    int i;
    for (i = 0; i < n; ++i) {
        graph[i].node_n = 2, graph[i].node_m = 2;
        graph[i].node = poa_init_node(2);
        graph[i].node[0].node_id = 0; graph[i].node[0].cumul_len = 0;
        graph[i].node[1].node_id = 1; graph[i].node[1].cumul_len = 0;
        graph[i].edge_n = 0; graph[i].edge_m = 0;
        graph[i].cons_seq = NULL;
    }
    return graph;
}

void poa_free_graph(poa_graph_t *graph, int n) {
    int i;
    for (i = 0; i < n; ++i) {
        if (graph[i].node_m > 0) poa_free_node(graph[i].node, graph[i].node_m);
        if (graph[i].edge_m > 0) poa_free_edge(graph[i].edge, graph[i].edge_m);
        if (graph[i].cons_seq) free(graph[i].cons_seq);
    }
    free(graph);
}

int poa_graph_node_id_to_rank(poa_graph_t *graph, int node_id) {
    if (node_id == 0) return 0;
    else return node_id-1;
}

int poa_graph_rank_to_node_id(poa_graph_t *graph, int rank_i) {
    if (rank_i == 0) return 0;
    else return rank_i+1;
}

//TODO 
//1. node_rank
//2. node.cumul_len
int poa_topological_sort(poa_graph_t *graph) {
    return 0;
}


int poa_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar) {
    if (graph->node_n <= 2 || qlen <= 0) { // empty graph or seq
        err_func_format_printf(__func__, "graph node: %d\tquery: %d\n", graph->node_n, qlen);
        return -1;
    }
#ifdef __DEBUG__
    int i;
    for (i = 0; i < qlen; ++i)
        printf("%c", "ACGTN"[query[i]]); 
    printf("\n");
#endif
    int score = -1;
    if (ppt->align_mode == POA_GLOBAL_FLAG) {
        score = poa_global_align_sequence_with_graph(graph, query, qlen, ppt, n_cigar, graph_cigar);
    } else if (ppt->align_mode == POA_LOCAL_FLAG) {
        score = poa_local_align_sequence_with_graph(graph, query, qlen, ppt, n_cigar, graph_cigar);
    } else if (ppt->align_mode == POA_EXTEND_FLAG) {
        score = poa_extend_align_sequence_with_graph(graph, query, qlen, ppt, n_cigar, graph_cigar);
    } else {
        err_fatal(__func__, "Unknown align mode: %d\n", ppt->align_mode);
    }
    return score;
}

int poa_add_graph_node(poa_graph_t *graph, uint8_t base) {
    int node_id = graph->node_n;
    graph = poa_realloc_graph_node(graph);
    // add node
    graph->node[node_id].node_id = node_id;
    graph->node[node_id].base = base;

    ++graph->node_n;
    return node_id;
}

int poa_add_graph_edge(poa_graph_t *graph, int from_id, int to_id, int check_edge) {
    if (from_id < 0 || from_id >= graph->node_n || to_id < 0 || to_id >= graph->node_n) err_fatal(__func__, "node_n: %d\tfrom_id: %d\tto_id: %d\n", graph->node_n, from_id, to_id);

    if (check_edge) {
        int i;
        for (i = 0; i < graph->node[from_id].out_edge_n; ++i) {
            if (graph->node[from_id].out_id[i] == to_id) // edge exists
                return 0;
        }
    }
    // add edge
    graph = poa_realloc_graph_edge(graph, 0, to_id); 
    graph->node[to_id].in_id[graph->node[to_id].in_edge_n] = from_id;
    ++graph->node[to_id].in_edge_n;
    graph = poa_realloc_graph_edge(graph, 1, from_id); 
    graph->node[from_id].out_id[graph->node[from_id].out_edge_n] = to_id;
    ++graph->node[from_id].out_edge_n;
    return 0;
}

int poa_add_graph_sequence(poa_graph_t *graph, uint8_t *seq, int seq_l, int start, int end) {
    if (seq_l <= 0 || start >= seq_l || end <= start) err_fatal(__func__, "seq_l: %d\tstart: %d\tend: %d\n", seq_l, start, end);
    if (start < 0) start = 0; if (end > seq_l) end = seq_l;
    int node_id = poa_add_graph_node(graph, seq[start]);
    poa_add_graph_edge(graph, POA_SRC_NODE_ID, node_id, 0);
    int i; 
    for (i = start+1; i < end; ++i) {
        node_id = poa_add_graph_node(graph, seq[i]);
        poa_add_graph_edge(graph, node_id-1, node_id, 0);
    }
    poa_add_graph_edge(graph, node_id, POA_SINK_NODE_ID, 0);
    
    return 0;
}

// fusion stratergy :
// 1. Match: merge to one node
// 2. Mismatch: check if B is identical to A' aligned nodes, then merge to node; if not, add node
// 3. Insertion: add node
// 4. Deletion: nothing
// 5. Clipping: add node
// 6. For all first/last node, link to virtual start/end node
int poa_add_graph_alignment(poa_graph_t *graph, uint8_t *seq, int seq_l, int n_cigar, poa_cigar_t *poa_cigar) {
    if (graph->node_n == 2) { // empty graph FIXME : when 
        return poa_add_graph_sequence(graph, seq, seq_l, 0, seq_l);
    } else {
        if (graph->node_n < 2) {
            err_fatal(__func__, "Graph node: %d.\n", graph->node_n);
        } else if (n_cigar == 0) {
            err_fatal(__func__, "Empty graph cigar.\n");
        }
    }
    // normal graph, normal graph_cigar
    int i, j; int op, len, node_id, query_id, last_new = 0, last_id = POA_SRC_NODE_ID, new_id;//, aligned_id;
    for (i = 0; i < n_cigar; ++i) {
        op = poa_cigar[i] & 0xf;
        if (op == POA_CMATCH) {
            node_id = (poa_cigar[i] >> 34) & 0x3fffffff;
            poa_add_graph_edge(graph, last_id, node_id, 1-last_new);
            last_id = node_id; last_new = 0;
        } else if (op == POA_CDIFF) {
            node_id = (poa_cigar[i] >> 34) & 0x3fffffff;
            query_id = (poa_cigar[i] >> 4) & 0x3fffffff;
            // check if query base is identical to node_id's aligned node XXX ??? XXX
            //if ((aligned_id = poa_get_aligned_id(graph, node_id, seq[query_id])) >= 0) {
            //    poa_add_graph_edge(graph, last_id, aligned_id);
            //    last_id = aligned_id;
            //}
            new_id = poa_add_graph_node(graph, seq[query_id]);
            poa_add_graph_edge(graph, last_id, new_id, 0);
            // add new_id to node_id's aligned node
            last_id = new_id; last_new = 1;
        } else if (op == POA_CINS || op == POA_CSOFT_CLIP || op == POA_CHARD_CLIP) {
            query_id = (poa_cigar[i] >> 34) & 0x3fffffff;
            len = (poa_cigar[i] >> 4) & 0x3fffffff;
            for (j = len-1; j >= 0; --j) { // XXX use dynamic id, instead of static query_id
                new_id = poa_add_graph_node(graph, seq[query_id-j]);
                poa_add_graph_edge(graph, last_id, new_id, 0);
                last_id = new_id; last_new = 1;
            }
        } else if (op == POA_CDEL) {
            // nothing;
            continue;
        }
    } poa_add_graph_edge(graph, last_id, POA_SINK_NODE_ID, 1-last_new);
    poa_topological_sort(graph);
    return 0;
}

uint8_t *poa_generate_consensus(poa_graph_t *graph, int *cons_seq_l) {
    if (graph->node_n <= 2) {
        err_func_format_printf(__func__, "graph node: %d\n", graph->node_n);
        *cons_seq_l = 0;
        return NULL;
    }
    uint8_t *cons_seq = (uint8_t*)_err_malloc(graph->node_n * sizeof(uint8_t));

    *cons_seq_l = 1;
    return cons_seq;
}
