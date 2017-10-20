mpicc kruskal.c -o krus

mpirun -np 8 --hosts master,client2 ./krus input.txt
