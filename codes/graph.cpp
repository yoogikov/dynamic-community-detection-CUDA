#include <fstream>
#include <iostream>
#include <vector>

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
        std::cout <<n_v << " " << n_e << " " << m << std::endl;
        for(int i=0;i<n_v+1;i++)std::cout << offset[i] << " ";std::cout << std::endl;
        for(int i=0;i<n_e;i++)std::cout << edges[i] << " ";std::cout << std::endl;
        for(int i=0;i<n_e;i++)std::cout << weights[i] << " ";std::cout << std::endl;

    }
    CUDA_ void initialize(int n_v, int n_e){
        this->n_v = n_v;
        this->n_e = n_e;
        m = 0;
        offset = (int*)malloc(sizeof(int)*(n_v+1));
        offset[0]=0;
        edges = (int*)malloc(sizeof(int)*n_e);
        weights = (int*)malloc(sizeof(int)*n_e);
    }

    CUDA_ void getm(){
        m = 0;
        for(int i=0;i<n_v;i++){
            for(int j=offset[i];j<offset[i+1];j++){
                m+=weights[j];
            }
        }
        m/=2;

        std::cout << "M found to be " << m << std::endl; 

    }
};

struct graph * input(char * inp){
    std::fstream fin;
    fin.open(inp, std::ios::in);
    int n_v, n_e;
    fin >> n_v >> n_e;
    
    graph * G = new graph();
    G->initialize(n_v, n_e);

    for(int i=0;i<n_v+1;i++){
        fin >> G->offset[i];
    }
    for(int i=0;i<n_e;i++){
        fin >> G->edges[i];
    }
    for(int i=0;i<n_e;i++){
        fin >> G->weights[i];
    }

    /*if(G->offset[0]==0 and G->offset[1]==0){
        for(int i=0;i<G->n_v;i++)G->offset[i] = G->offset[i+1];
        G->n_v-=1;
        for(int i=0;i<G->n_e;i++)G->edges[i]-=1;
    }*/


    for(int i=0;i<n_v;i++){
        for(int j=G->offset[i];j<G->offset[i+1];j++){
            if(i<=G->edges[j])G->m+=G->weights[j];
        }
    }
    

    return G;
};

struct graph * convert(std::vector<std::vector<std::pair<int, int>>> &adlist){
    graph * G = new graph();
    int n_v = adlist.size();
    int n_e = 0;
    for(int i=0;i<n_v;i++){
        n_e += adlist[i].size();
    }

    G->initialize(n_v, n_e);
        
    int offset = 0;
    G->offset[0] = offset;
    for(int i=0;i<n_v;i++){
        for(int j=0;j<adlist[i].size();j++){
            G->edges[offset] = adlist[i][j].first;
            G->weights[offset++] = adlist[i][j].second;
        }
        G->offset[i+1] = offset;
    }
    G->getm();
    //std::cout << G->m <<"yes"<< std::endl;
    return G;
}
