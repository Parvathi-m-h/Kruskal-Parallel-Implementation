#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>

typedef enum { FALSE, TRUE } boolean;

void parse_input();
void parse_edge_list_input();
int compare_edges(const void *a, const void *b);
boolean should_send();
int get_merge_partner_rank(boolean);
void merge();
void send_recieve_local_mst();
void merge_mst();
void print_received_mst_edges();
void start_sort_measure();
void end_sort_measure();
void start_communication_measure();
void end_communication_measure();

int number_of_processors;  
int rank;                  
int root = 0;              
MPI_Datatype mpi_edge;     

FILE * f;                  

int number_of_vertices;    
int vertices_per_process;  
int number_of_edges;      

int merge_iteration;       /* Current merge iteration */
int max_merge_iteration;   /* Max merge iteration (if there is 2^n processors
                              there will be n merge iterations) */
boolean has_sent;         

typedef struct edge_t {    
    int v;
    int u;
    int weight;
} edge_s;

edge_s * edges;            /* Array of edges belonging to this process */
edge_s * local_mst_edges;  /* Array of edges which form local MST */
edge_s * recv_mst_edges;   /* Array of edges received from other process */
edge_s * merged_mst_edges; /* Array of merged local and received edges */
int local_mst_edge_count;  /* Number of edges in array of local MST edges */
int recv_mst_edge_count;   /* Number of edges in array of received edges */
int merged_mst_edge_count; /* Number of edges in array of merged edges */

typedef struct node_t {    /* Union-Find data structure */
    struct node_t * parent;
    int depth;
} u_node;

u_node * uf_set;           /* Array indicating which set does vertex belong to
                              (used in Union-Find algorithm) */

void uf_make();
u_node * uf_find(u_node * a);
void uf_union(u_node * a, u_node * b);

double comp_start_time;    /* Start time of computation */
double comp_end_time;      /* End time of computation */

double sort_start_time;    /* Start time of qsort-increasing order  */
double sort_end_time;      /* End time of qsort-increasing order */

double comm_start_time;    /* Start time of communication */
double comm_end_time;      /* End time of communication */
double total_comm_time;    /* Total time of communication */

double parse_start_time;   /* Start time of parsing edges*/
double parse_end_time;     /* End time of parsing edges*/

void initialisation(int argc, char** argv) {

	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processors);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	/* Parsing graph size */
	f = fopen(argv[1], "r");
	//fread(&number_of_vertices, sizeof(number_of_vertices), 1, f);
          fscanf(f,"%d ",&number_of_vertices);
	/* Create 'edge' type (u, v, distance) */
	MPI_Type_contiguous(3, MPI_INT, &mpi_edge);
	MPI_Type_commit(&mpi_edge);

	vertices_per_process = number_of_vertices / number_of_processors;

	total_comm_time = 0;

	if (rank == root) {
		printf("Number of vertices: %d. Each processor will get %d vertices.\n", number_of_vertices, vertices_per_process);
	}

}

void check_conditions() {

	/* Check for error conditions */
	unsigned char isError = 0;
	if (number_of_vertices < number_of_processors) {
		printf("Number of vertices (%d) is smaller than number of processors (%d)\n", number_of_vertices, number_of_processors);
		isError = 1;
	}

	if (number_of_vertices % number_of_processors != 0) {
		puts("Number of vertices % number of processors must be 0\n");
		isError = 1;
	}

	if (!((number_of_processors != 0) && !(number_of_processors & (number_of_processors - 1)))) {
		puts("Number of processors must be power of 2\n");
		isError = 1;
	}

	if (isError == 1) {
		fclose(f);
		MPI_Type_free(&mpi_edge);
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

}

void parse_input() {

	/* Start and end index of vertex belonging to this process */
	int rank_range_start = rank * vertices_per_process;
	int rank_range_end = (rank + 1) * vertices_per_process;

	number_of_edges = 0;

	int total_edges, u, v, edge_weight;
	//fread(&total_edges, sizeof(int), 1, f);
	fscanf(f,"%d\n",&total_edges);
	/* Create arrays for holding edge structures */
	edges = (edge_s *) malloc(vertices_per_process * number_of_vertices * sizeof(edge_s));
	local_mst_edges = (edge_s *) malloc((number_of_vertices - 1) * sizeof(edge_s));
	merged_mst_edges = (edge_s *) malloc(2 * (number_of_vertices - 1) * sizeof(edge_s));

	int i;int k=0;
	int f_edge[3];
	
	for (i=0; i < total_edges; i++) {
		
		for(k=0;k<3;k++)
		fscanf(f,"%d ",&f_edge[k]);
		//fread(&f_edge, sizeof(int), 3, f);
		
		v = f_edge[0]; u = f_edge[1]; edge_weight = f_edge[2];
		
		if (((v >= rank_range_start) && (v < rank_range_end))
		|| ((u >= rank_range_start) && (u < rank_range_end))) {
			edge_s e = { v, u, edge_weight };
			edges[number_of_edges++] = e;
		}
	}
	
	fclose(f);

	MPI_Barrier(MPI_COMM_WORLD);
}

void finalize() {

	/* Free memory used by dynamic data structures */
	free(uf_set);
	free(edges);
	free(local_mst_edges);
	if (rank % 2 == 0) {
		free(recv_mst_edges);
	}
	free(merged_mst_edges);

	/* Shut down MPI */
	MPI_Type_free(&mpi_edge);
	MPI_Finalize();

	//printf("Proc %d: All done. Terminating.\n", rank);
}

int compare_edges(const void *a, const void *b) {
	edge_s *ea = (edge_s *)a;
	edge_s *eb = (edge_s *)b;

	if (ea->weight < eb->weight) {
		return -1;
	} else if (ea->weight > eb->weight) {
		return 1;
	}

	// Never happens if weights are distinct
	return 0;
}

void find_local_mst() {

	local_mst_edge_count = 0;
	merged_mst_edge_count = 0;

	/* Sort edges belonging to each process */
	start_sort_measure();
	qsort(edges, number_of_edges, sizeof(edge_s), compare_edges);
	end_sort_measure();

	uf_make();

	int used_edge_index = 0;
	int i;
	for(i = 1;i<=number_of_edges;i++) {

		edge_s * min_edge = &edges[used_edge_index++];

		int v = (*min_edge).v;
		int u = (*min_edge).u;
		u_node * v_node = uf_find(uf_set + v);
		u_node * u_node = uf_find(uf_set + u);

		/* Add edge to MST if it doesn't form a cycle */
		if(v_node != u_node) {
			local_mst_edges[local_mst_edge_count++] = *min_edge;
			merged_mst_edges[merged_mst_edge_count++] = *min_edge;

			uf_union(v_node, u_node);
		}
	}

}

void merge_mst() {

	max_merge_iteration = 0;

	int temp = number_of_processors;
	while(temp > 1) {
		temp = temp / 2;
		max_merge_iteration++;
	}

	/* Skip odd ranks, as they will never receive, everyone else will receive at least once */
	if (rank % 2 == 0) {
		recv_mst_edges = (edge_s *) malloc((number_of_vertices - 1) * sizeof(edge_s));
	}

	has_sent = FALSE;
	merge_iteration = 0;
	while(merge_iteration < max_merge_iteration) {
		send_recieve_local_mst();
		merge_iteration++;
		if (has_sent) { // Process has sent it's data and can terminate now
			break;
		}
	}

}

void send_recieve_local_mst() {

	MPI_Status status;
	boolean sending = should_send();
	int merge_partner_rank = get_merge_partner_rank(sending);
	if (sending) { // should send
	        char processor_name[MPI_MAX_PROCESSOR_NAME];
	        int name_len;
                MPI_Get_processor_name(processor_name, &name_len);
		printf("Proc name %s with rank %d: Sending %d edges to processor #%d\n",processor_name, rank, local_mst_edge_count, merge_partner_rank);
		MPI_Send(local_mst_edges, local_mst_edge_count, mpi_edge, merge_partner_rank, 0, MPI_COMM_WORLD);
		has_sent = TRUE;
	} else { // should receive
		start_communication_measure();
		MPI_Recv(recv_mst_edges, number_of_vertices - 1, mpi_edge, merge_partner_rank, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, mpi_edge, &recv_mst_edge_count);
		end_communication_measure();
		if (recv_mst_edge_count != MPI_UNDEFINED) {
			char processor_name[MPI_MAX_PROCESSOR_NAME];
	                int name_len;
                        MPI_Get_processor_name(processor_name, &name_len);
			printf("Proc name %s with rank %d: Received %d edges from processor #%d\n", processor_name, rank, recv_mst_edge_count, merge_partner_rank);
			merge();
		} else {
			MPI_Abort(MPI_COMM_WORLD, 10); // Abort all
		}
	}

}

boolean should_send() {
	if ((rank / (int)pow(2, merge_iteration)) % 2 == 0) {
		return FALSE;
	} else {
		return TRUE;
	}
}

int get_merge_partner_rank(boolean should_send) {
	// Merge partner is rank +/- 2^merge_iteration
	if (should_send) {
		return rank - pow(2, merge_iteration); // rank of destination
	} else {
		return rank + pow(2, merge_iteration); // rank of source
	}
}

void merge() {

	int i;
	for (i = 0; i < recv_mst_edge_count; i++) {
		edge_s received_edge = recv_mst_edges[i];
		merged_mst_edges[merged_mst_edge_count++] = received_edge;
	}

	/* Sort local and received edges */
	qsort(merged_mst_edges, merged_mst_edge_count, sizeof(edge_s), compare_edges);

	uf_make();

	int used_edge_index = 0;
	local_mst_edge_count = 0;

	for (i = 0; i < merged_mst_edge_count; i++) {

		edge_s * min_edge = &merged_mst_edges[used_edge_index++];

		int v = (*min_edge).v;
		int u = (*min_edge).u;
		u_node * v_node = uf_find(uf_set + v);
		u_node * u_node = uf_find(uf_set + u);

		if(v_node != u_node) {
			local_mst_edges[local_mst_edge_count++] = *min_edge;

			uf_union(v_node, u_node);
		}

	}

	/* Transfer new local MST edges to merged edges */
	merged_mst_edge_count = 0;
	for (i = 0; i < local_mst_edge_count; i++) {
		edge_s merged_edge = local_mst_edges[i];
		merged_mst_edges[merged_mst_edge_count++] = merged_edge;
	}

}

void uf_make() {
	int size = sizeof(u_node) * number_of_vertices;// + sizeof(int)*(number_of_vertices - 1);
	uf_set = malloc(size);
	memset(uf_set, 0, size);
}

u_node * uf_find(u_node * a) {
    if (a->parent==NULL) return a;
    else return (a->parent = uf_find(a->parent));
}

void uf_union(u_node * a, u_node * b) {
    if (a->depth > b->depth) {
        b->parent = a;
    } else if (a->depth<b->depth) {
        a->parent = b;
    } else {
        a->parent = b;
        a->depth += 1;
    }
}

void print_local_edges() {
	int i;
	printf("Proc %d local edges:\n", rank);
	for (i = 0; i < number_of_edges; ++i) {
		printf("(%d,%d) = %d\n", edges[i].v, edges[i].u, edges[i].weight);
	}
}

void print_local_mst_edges() {
	int i;
	printf("Proc %d local MST edges:\n", rank);
	for(i = 0; i < local_mst_edge_count; i++) {
		printf("(%d,%d) = %d\n",local_mst_edges[i].v, local_mst_edges[i].u, local_mst_edges[i].weight);
	}
}

void print_received_mst_edges() {
	int i;
	printf("Proc %d received MST edges:\n", rank);
	for(i = 0; i < recv_mst_edge_count; i++) {
		printf("(%d,%d) = %d\n",recv_mst_edges[i].v, recv_mst_edges[i].u, recv_mst_edges[i].weight);
	}
}

void print_merged_mst_edges() {
	int i;
	printf("Proc %d merged Final MST edges:\n", rank);
	for(i = 0; i < merged_mst_edge_count; i++) {
		printf("(%d,%d) = %d\n",merged_mst_edges[i].v, merged_mst_edges[i].u, merged_mst_edges[i].weight);
	}
}

void start_computation_measure() {
	comp_start_time = MPI_Wtime();
}

void end_computation_measure() {
	comp_end_time = MPI_Wtime();
}

void start_sort_measure() {
	sort_start_time = MPI_Wtime();
}

void end_sort_measure() {
	sort_end_time = MPI_Wtime();
}

void start_communication_measure() {
	comm_start_time = MPI_Wtime();
}

void end_communication_measure() {
	comm_end_time = MPI_Wtime();
	total_comm_time += comm_end_time - comm_start_time;
}

void start_parse_measure() {
	parse_start_time = MPI_Wtime();
}

void end_parse_measure() {
	parse_end_time = MPI_Wtime();
}

void print_measurement_times() {
	if (rank == root) {
		printf("Total: %f, Sort: %f, Comm: %f, Parse: %f\n",
		comp_end_time - comp_start_time,
		sort_end_time - sort_start_time,
		total_comm_time,
		parse_end_time - parse_start_time);
	}
}

int main(int argc, char** argv) {

	if (argc < 2) {
		printf("usage: %s filename\n", argv[0]);
		return 1;
	}

	initialisation(argc, argv);

	check_conditions();

	start_parse_measure();
	parse_input();
	end_parse_measure();

	start_computation_measure();
	find_local_mst();
	merge_mst();
	end_computation_measure();

	print_measurement_times();

	if (rank == root) {
		print_merged_mst_edges();
	}
	
	else
	{
	   print_local_mst_edges();
	}

	finalize();

	return 0;
}

