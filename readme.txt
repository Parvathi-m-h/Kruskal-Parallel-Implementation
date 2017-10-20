mpicc kruskal.c -o krus

mpirun -np 8 --hosts master,client2 ./krus input.txt

//input.txt file format:
n m // n->no of vertices m->no of edges
v1 v2 e //v1->first vertex v2->second vertex  e->edge weight between the 2 vertices