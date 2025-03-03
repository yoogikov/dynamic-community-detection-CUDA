#include "graph.cpp"
#include <set>
#include <vector>

#define thresh1 0.000001

double modularity(graph * G, std::vector<int> &colour){
    int m = 0;
    for(int i=0;i<G->n_e;i++){
        m+=G->weights[i];
    }
    m/=2;
    std::vector<int>degrees(G->n_v);

    for(int i=0;i<G->n_v;i++){
        for(int j = G->offset[i];j<G->offset[i+1];j++){
            degrees[i] += G->weights[j];
        }
    }
    
    double sub = 0;
    for(int i=0;i<G->n_v;i++){
        for(int j=0;j<G->n_v;j++){
           if(colour[i]==colour[j])sub += (double)(degrees[i]*degrees[j]/(2.0*(double)m));
        }
    }

    double add = 0;
    for(int i=0;i<G->n_v;i++){
        for(int j=G->offset[i];j<G->offset[i+1];j++){
            if(colour[i]==colour[G->edges[j]])add += (double)(G->weights[j]);
        }
    }

    return (add-sub)/((double)2.0*m);
};

double modularity_colour(graph * G, std::vector<int> &colour, int col){
    double add=0, sub=0;

    for(int i=0;i<G->n_v;i++){
        if(colour[i]!=col) continue;
        for(int j=G->offset[i];j<G->offset[i+1];j++){
            sub += G->weights[j];
            if(colour[G->edges[j]]==col)add += G->weights[j];
        }
    }

    sub = sub/((double)2*(G->m));
    add = add/((double)2*(G->m));
    return add - sub*sub;
}

double summed_modularity(graph * G, std::vector<int> &colour){
    std::set<int> cols;
    for(auto i:colour)cols.insert(i);
    double Q = 0;
    for(auto i:cols)Q+=modularity_colour(G, colour, i);
    return Q;
}

void initialize(graph * G, std::vector<int> &colour, std::vector<int> &tot_degree, std::vector<int> &degree){
    for(int src=0;src<G->n_v;src++){
        int start = G->offset[src];
        int end = G->offset[src+1];
        colour[src] = src;
        for(int i=start;i<end;i++){
            tot_degree[src]+=G->weights[i];
            degree[src]+=G->weights[i];
        }
    }
}


double optimize(graph * G, std::vector<int> &colour, std::vector<int> &tot_degree, std::vector<int> &degree){
    double d = 0;
    for(int i=0;i<G->n_v;i++){
        int best_new_colour = colour[i];
        double max_delta_Q = 0;

        int kx=0, wa=0, ea=0;
        std::vector<int> col_degrees(G->n_v);

        for(int j=G->offset[i];j<G->offset[i+1];j++){
            kx += G->weights[j];
            col_degrees[colour[G->edges[j]]] += G->weights[j];
        }

        for(int j=G->offset[i];j<G->offset[i+1];j++){
            int s_col = colour[i];
            int d_col = colour[G->edges[j]];

            int wa=col_degrees[s_col];
            int wb=col_degrees[d_col];


            double delta_Q1 = static_cast<double>(wb-wa)/G->m - (static_cast<double>(degree[i])/(2*G->m))*(static_cast<double>(2*tot_degree[d_col] - 2*tot_degree[s_col] + 2*degree[i] - ((s_col==d_col)?(2*kx):0))/(2*G->m));

            if(delta_Q1>max_delta_Q){
                max_delta_Q = delta_Q1;
                best_new_colour = d_col;
            }
        }
        d += max_delta_Q;

        int s_col = colour[i];
        int d_col = best_new_colour;

        colour[i] = best_new_colour;

        tot_degree[s_col]-=degree[i];
        tot_degree[d_col]+=degree[i];
    }
    return d;
}

double moving_phase(graph * G, std::vector<int> &colour, std::vector<int> &tot_degree, std::vector<int> &degree){
    double delta = 0;
    double d = 100;
    while(d>=thresh1){
        d = optimize(G, colour, tot_degree, degree);
        delta+=d;
        std::set<int> col;
        for(auto i:colour)col.insert(i);
        std::cout << col.size() << std::endl;
    }
    return delta;
}





int main(int argc, char *argv[]){
    graph * G = input(argv[1]);
//  G->display();
    std::vector<int> colour(G->n_v);
    std::vector<int> tot_degree(G->n_v);
    std::vector<int> degree(G->n_v);


    initialize(G, colour, tot_degree, degree);

    std::cout << G->n_v << std::endl;

    double moving_phase_delta = moving_phase(G, colour, tot_degree, degree);    
    


}
