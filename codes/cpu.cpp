#include "graph.cpp"
#include <map>
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

            //std::cout << s_col << " " << d_col << std::endl;

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
        //std::cout << col.size() << std::endl;
    }
    return delta;
}


graph* aggregate_phase(graph*G, std::vector<int> &colour, std::vector<int> &tot_degree, std::vector<int> &degree){
    
    std::vector<int> renumber(G->n_v, -1);

    int new_col = 0;
    for(int i=0;i<G->n_v;i++){
        if(renumber[colour[i]]!=-1)continue;
        renumber[colour[i]]=new_col++;
    }
    
    std::vector<std::vector<std::pair<int, int>>> new_adlist(new_col);
    for(int i=0;i<G->n_v;i++){
        for(int j=G->offset[i];j<G->offset[i+1];j++){
            int src = i;
            int end = G->edges[j];
            int weight = G->weights[j];
            new_adlist[renumber[colour[src]]].push_back(std::make_pair(renumber[colour[end]], weight));
        }
    }


    std::vector<std::vector<std::pair<int, int>>> compressed_adlist(new_col);
    for(int comm=0;comm<new_col;comm++){
        std::map<int, int> count;
        for(int i=0;i<new_adlist[comm].size();i++){
            count[new_adlist[comm][i].first] += new_adlist[comm][i].second;
        }
        for(auto p:count){
            compressed_adlist[comm].push_back(p);
        }
    }

    graph * G_ = convert(compressed_adlist);
    return G_;

    
}

double w_degree(graph * G, int n){
    double d = 0;
    for(int i=G->offset[n];i<G->offset[n+1];i++)d += G->weights[i];
    return d;
}

double tot_degree(graph * G, std::vector<int> &colour, int col){
    double d = 0;
    for(int i=0;i<G->n_v;i++){
        if(colour[i]==col)d+=w_degree(G, i);
    }
    return d;
}

double comp_weight(graph * G, std::vector<int> &colour, int n, int col){
    double d = 0;
    for(int i=G->offset[n];i<G->offset[n+1];i++)if(colour[G->edges[i]]==col)d+=G->weights[i];
    return d;
}

bool one_level(graph * G, std::vector<int> &colour, std::vector<int> &tot_degrees, std::vector<int> &degrees){
    //missing randomization

    double new_mod = modularity(G, colour);
    //std::cout << "Modularity is " << new_mod << std::endl;
    double cur_mod = new_mod;
    bool improvement = false;
    int moves = 0;
    do {
        double cur_mod = new_mod;
        moves = 0;
        for(int src = 0;src<G->n_v;src++){
            int node = src;
            int node_comm = colour[node];
            double w_deg = w_degree(G, node);

            colour[node] = -1;
            tot_degrees[node_comm] -= w_deg;

            int best_comm = node_comm;
            double best_increase = -1e9;
            int neigh_col = best_comm;

            unsigned int i = G->offset[node];
            do{
                
                double totc = tot_degrees[neigh_col];
                double dnc = comp_weight(G, colour, node, neigh_col);
                double degc = w_deg;
                double m2 = 2 * (double)G->m;



                //std::cout << totc << " " << degc << " " << m2 << " " << dnc << std::endl;
                double increase = dnc - totc*degc/m2;

                //std::cout << node << " " << neigh_col << " " << increase << std::endl;
                if(increase>best_increase){
                    best_comm = neigh_col;
                    best_increase = increase;
                }

                if(i==G->offset[node+1])break;
                neigh_col = colour[G->edges[i]];
                if(neigh_col==-1)neigh_col=node_comm;

                i++;

            }while(i<=G->offset[node+1]);

            //std::cout << node << " " << best_comm << std::endl;

            colour[node] = best_comm;
            tot_degrees[colour[node]]+=w_deg;
            
            if(best_comm!=node_comm)moves++;
            if(moves>0)improvement=true;
            
            
        }
    } while(moves>0);
    return improvement;
}



int main(int argc, char *argv[]){
    graph * G = input(argv[1]);
    std::vector<int> colouring(G->n_v);
    for(int i=0;i<G->n_v;i++)colouring[i] = i;  
//  G->display();
    for(int i=0;i<5;i++){
        std::vector<int> colour(G->n_v);
        std::vector<int> tot_degree(G->n_v);
        std::vector<int> degree(G->n_v);


        initialize(G, colour, tot_degree, degree);

        for(int i=0;i<G->n_v;i++)colour[i]=i; 
        //double moving_phase_delta = moving_phase(G, colour, tot_degree, degree);    
        one_level(G, colour, tot_degree, degree);
        //for(int i=0;i<G->n_v;i++)std::cout << i << " " << colour[i] << std::endl; 

        graph* G_ = aggregate_phase(G, colour, tot_degree, degree);

        //G_->display();

        G = G_;
        G->display();

        G->getm();

        std::cout << modularity(G, colour) << std::endl;

        std::cout << "--------------" << std::endl;
    }



//    G_->display();


}
