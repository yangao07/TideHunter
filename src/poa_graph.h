#ifndef POA_GRAPH_H
#define POA_GRAPH_H

#include <stdint.h>
#include "poa.h"

#ifdef __cplusplus
extern "C" {
#endif

//#define CIGAR_STR "MIDNSHP=XB"
//#define POA_GRAPH_CIGAR_STR "=XIDNSH"
//#define POA_GRAPH_CEQUAL 0
//#define POA_GRAPH_CMISMATCH 1
//#define POA_GRAPH_CINS 2
//#define POA_GRAPH_CDEL 3
//#define POA_GRAPH_CREF_SKIP 4
//#define POA_GRAPH_CCLIP 5

#define POA_SRC_NODE_ID 0
#define POA_SINK_NODE_ID 1

typedef struct {
    int node_id, rank, level, cumul_len;
    int in_edge_n, in_edge_m, *in_id;
    int out_edge_n, out_edge_m, *out_id;
    int aligned_node_n, aligned_node_m, *aligned_node_id; // mismatch; aligned node will have same level
    uint8_t base; // 0~m
    // ID, pos ???
} poa_node_t;

typedef struct {
    int from_id, to_id;
    int weight;
} poa_edge_t;

typedef struct {
    int node_n, node_m; poa_node_t *node; 
    int *rank_to_node_id, *node_id_to_rank;
    int *level_to_node_id, *node_id_to_level, *level_index, level_n; // level: for anti-diagnal band width
    int edge_n, edge_m; poa_edge_t *edge;
    int seq_n;
    uint8_t *cons_seq;
} poa_graph_t;

// XXX max of in_edge is pow(2,30)
// for MATCH/MISMATCH: node_id << 34  | query_id << 4 | op
// for INSERTION:      query_id << 34 | op_len << 4   | op
// for DELETION:       node_id << 34  | op_len << 4   | op // op_len is always equal to 1
// for CLIP            query_id << 34 | op_len << 4   | op 
#define poa_cigar_t int64_t 
                       
//typedef struct {
//    uint8_t op; // 0:match, 1:insertion, 2:deletion, 3:mismatch
//    int32_t len, m, *node_id, *query_id;
//} poa_graph_cigar_t; // alignment result of mapping a sequence to a graph


poa_node_t *poa_init_node(int n);
void poa_free_node(poa_node_t *node, int n);
poa_edge_t *poa_init_edge(int n);
void poa_free_edge(poa_edge_t *edge, int n);
poa_graph_t *poa_init_graph(int n);
void poa_free_graph(poa_graph_t *graph, int n);
int poa_align_sequence_with_graph(poa_graph_t *graph, uint8_t *query, int qlen, poa_para_t *ppt, int *n_cigar, poa_cigar_t **graph_cigar);
int poa_add_graph_alignment(poa_graph_t *graph, uint8_t *query, int qlen, int n_cigar, poa_cigar_t *poa_cigar);
int poa_graph_node_id_to_rank(poa_graph_t *graph, int node_id);
int poa_graph_rank_to_node_id(poa_graph_t *graph, int rank_i);
uint8_t *poa_generate_consensus(poa_graph_t *graph, int *cons_seq_l);
#ifdef __cplusplus
}
#endif

#endif
