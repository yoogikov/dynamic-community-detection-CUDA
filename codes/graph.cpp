#include <fstream>
#include <iostream>

#ifdef __CUDACC__
#define CUDA_ __host__ __device__
#else
#define CUDA_
#endif 


class graph{
public:
    CUDA_ int n_v, n_e, m;
    CUDA_ int * offset;
    CUDA_ int * edges;
    CUDA_ int * weights;
    CUDA_ void display(){
        for(int i=0;i<n_v+1;i++)std::cout << offset[i] << " ";std::cout << std::endl;
        for(int i=0;i<n_e;i++)std::cout << edges[i] << " ";std::cout << std::endl;
        for(int i=0;i<n_e;i++)std::cout << weights[i] << " ";std::cout << std::endl;

    }
};

struct graph * input(char * inp){
    std::fstream fin;
    fin.open(inp, std::ios::in);
    int n_v, n_e;
    fin >> n_v >> n_e;
    
    graph * G = new graph();
    G->n_v = n_v;
    G->n_e = n_e;
    G->offset = (int *)malloc(sizeof(int)*(n_v+1));
    G->edges = (int *)malloc(sizeof(int)*n_e);
    G->weights = (int *)malloc(sizeof(int)*n_e);

    for(int i=0;i<n_v+1;i++){
        fin >> G->offset[i];
    }
    for(int i=0;i<n_e;i++){
        fin >> G->edges[i];
    }
    for(int i=0;i<n_e;i++){
        fin >> G->weights[i];
        G->m += G->weights[i];
    }
    
    G->m = G->m/2;

    return G;
};
