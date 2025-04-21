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
    CUDA_ int * community;
    CUDA_ int * in;
    CUDA_ int * tot;
    CUDA_ int * self_loops;
    CUDA_ bool weighted;
    CUDA_ void display(){
        std::cout <<n_v << " " << n_e << " " << m << std::endl;
        for(int i=0;i<n_v+1;i++)std::cout << offset[i] << " ";std::cout << std::endl;
        for(int i=0;i<n_e;i++)std::cout << edges[i] << " ";std::cout << std::endl;
        for(int i=0;i<n_e;i++)std::cout << weights[i] << " ";std::cout << std::endl;
        /*for(int i=0;i<n_v;i++)std::cout << community[i] << " ";std::cout << std::endl;*/
        /*for(int i=0;i<n_v;i++)std::cout << in[i] << " ";std::cout << std::endl;*/
        /*for(int i=0;i<n_v;i++)std::cout << tot[i] << " ";std::cout << std::endl;*/
        /*for(int i=0;i<n_v;i++)std::cout << self_loops[i] << " ";std::cout << std::endl;*/

    }
    CUDA_ void initialize(int n_v, int n_e){
        weighted = false;

        this->n_v = n_v;
        this->n_e = n_e;
        m = 0;
        offset = (int*)malloc(sizeof(int)*(n_v+1));
        offset[0]=0;
        edges = (int*)malloc(sizeof(int)*n_e);
        weights = (int*)malloc(sizeof(int)*n_e);
        community = (int*)malloc(sizeof(int)*n_v);
        in = (int*)malloc(sizeof(int)*n_v);
        tot = (int*)malloc(sizeof(int)*n_v);
        self_loops = (int*)malloc(sizeof(int)*n_v);
    }

    CUDA_ void getm(){
        m = 0;
        for(int i=0;i<n_v;i++){
            for(int j=offset[i];j<offset[i+1];j++){
                m+=weights[j];
            }
        }
        m/=2;


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

    for(int i=0;i<G->n_v;i++){
        G->community[i] = i;
        G->tot[i] = 0;
        G->in[i] = 0;
        for(int j=G->offset[i];j<G->offset[i+1];j++){
            G->tot[i] += G->weights[j];
            if(G->edges[j]==i)G->in[i] += G->weights[j];
        }
        G->self_loops[i] = G->in[i];
    }

    return G;
};

struct graph * convert(std::vector<std::vector<std::pair<int, int>>> &adlist){
    graph * G = new graph();


    int number_of_edges = 0;
    for(int i=0;i<adlist.size();i++){
        number_of_edges += adlist[i].size();
    }

    G->n_v = adlist.size();
    G->n_e = number_of_edges;
    G->m = 0;
    G->offset = (int*)malloc(sizeof(int)*(G->n_v+1));
    G->offset[0]=0;

    int n_v = G->n_v;

    G->edges = (int*)malloc(sizeof(int)*G->n_e);
    G->weights = (int*)malloc(sizeof(int)*G->n_e);
    G->community = (int*)malloc(sizeof(int)*G->n_v);
    G->in = (int*)malloc(sizeof(int)*n_v);
    G->tot = (int*)malloc(sizeof(int)*n_v);
    G->self_loops = (int*)malloc(sizeof(int)*n_v);

    int off = 0;
    for(int i=0;i<n_v;i++){
        G->community[i] = i;
        G->self_loops[i] = 0;
        G->in[i] = 0;
        G->tot[i] = 0;
        for(std::pair<int, int> p:adlist[i]){
            G->edges[off] = p.first;
            G->weights[off] = p.second;
            G->tot[i]+=p.second;
            if(p.first==i)G->in[i]+=p.second;
            if(p.first==i)G->self_loops[i]+=p.second;
            off++;
        }
        G->offset[i+1]=off;
    }
    G->weighted=true;
    return G;
}
