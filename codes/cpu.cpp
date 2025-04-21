#include "graph.cpp"
#include <set>
#include <vector>

#define thresh1 0.000001

double modularity(graph * G, std::vector<int> &colour){

    double q = 0.;
    double m2 = (double)G->m;
    m2*=2;
    
    for(int i=0;i<G->n_v;i++){
        if(G->tot[i]>0)
            q += (double)G->in[i]/m2 - ((double)G->tot[i]/m2)*((double)G->tot[i]/m2);
    }

    return q;
};

int uniqueness(graph * G){
    long long int x = 1;
    for(int i=0;i<G->n_v;i++){
        int tot_degree=1;
        for(int j=G->offset[i]; j<G->offset[i+1];j++){
            tot_degree+=G->weights[j];
        }
        x*=tot_degree;
        x = x%((long long int)1e9+7);
    }
    return x;
}

double single_modularity(graph *G){
    double m = (double)G->m;
    double q = 0.;

    for(int i=0;i<G->n_v;i++){
        int in_degree = 0;
        int tot_degree = 0;
        for(int j=G->offset[i]; j<G->offset[i+1];j++){
            tot_degree += G->weights[j];
            if(G->edges[j]==i)in_degree += G->weights[j];
        }
        q += (double)in_degree - (double)(tot_degree*tot_degree)/(2*m);
    }
    return q/(2*m);
}

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


graph* aggregate_phase(graph*G, std::vector<int> &colour, std::vector<int> &tot_degree, std::vector<int> &degree, std::vector<int> &colouring){

    
    std::vector<int> renumber(G->n_v, -1);

    for(int i=0;i<G->n_v;i++){
        renumber[G->community[i]]++;
    }


    int new_col = 0;
    for(int i=0;i<G->n_v;i++){
        if(renumber[i]!=-1)
            renumber[i]=new_col++;
    }

    std::vector<std::vector<int>> communities(new_col);
    for(int i=0;i<G->n_v;i++){
        communities[renumber[G->community[i]]].push_back(i);
    }

    std::vector<std::vector<std::pair<int, int>>> compressed_adlist(new_col);



    std::vector<int> n_cols(new_col);
    std::vector<int> n_weights(new_col, -1);
    int n_last=0;
    for(int comm=0;comm<communities.size();comm++){
        for(auto node:communities[comm]){
            for(int j = G->offset[node]; j<G->offset[node+1];j++){
                int dst = G->edges[j];
                int wt = G->weights[j];
                if(n_weights[renumber[G->community[dst]]]==-1){
                    n_cols[n_last++] = renumber[G->community[dst]];
                    n_weights[renumber[G->community[dst]]] = 0;
                }
                n_weights[renumber[G->community[dst]]] += wt;
            }
        }
        for(int i=0;i<n_last;i++){
            compressed_adlist[comm].push_back(std::make_pair(n_cols[i], n_weights[n_cols[i]]));
            n_weights[n_cols[i]] = -1;
        }
        n_last = 0;
    }


    graph * G_ = convert(compressed_adlist);
    std::cout << G_->n_v << std::endl;
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

    bool improvement = false;
    int moves = 0;
    double new_mod = modularity(G, colour);

    std::vector<int> n_weights(G->n_v, -1);
    std::vector<int> n_cols(G->n_v);
    int n_last = 0;
    double cur_mod = new_mod;

    int best_nblinks = 0;

    

    do {
        cur_mod = new_mod;
        moves = 0;


        for(int src = 0;src<G->n_v;src++){




            int node = src;
            int node_comm = G->community[node];
            double w_deg = degrees[node];


            for(int i=0;i<n_last;i++){
                n_weights[n_cols[i]]=-1;
            }
            n_last=0;

            n_cols[0] = node_comm;
            n_weights[node_comm]=0;
            n_last=1;

            std::pair<int*, int*> p = std::make_pair(G->edges+G->offset[src], G->weights + G->offset[src]);

            int * ptr = G->community;

            int deg = G->offset[src+1] - G->offset[src];

            for(int i=0;i<deg;i++){
                int dst = *(p.first + i);
                int wt = G->weighted?*(p.second+i):1;
                int n_colour = *(ptr+dst);

                if(dst!=node){
                    if(n_weights[n_colour]==-1){
                        n_weights[n_colour] = 0;
                        n_cols[n_last++] = n_colour;
                    }
                    n_weights[n_colour]+=wt;
                }
            }

            G->tot[node_comm] -= degrees[node];
            G->in[node_comm] -= 2*n_weights[node_comm] + G->self_loops[node]; 
            G->community[node] = -1;

            int best_comm = node_comm;
            double best_increase = 0.;

            for(int i=0;i<n_last;i++){
                int neigh_col = n_cols[i];

                double totc = G->tot[n_cols[i]];
                double dnc = n_weights[n_cols[i]];
                double degc = w_deg;
                double m2 = 2*(double)G->m;

                //std::cout << node << " " << neigh_col << " " << totc << ' ' << dnc << " " << degc << " " << m2 << "\n";

                double increase = dnc-totc*degc/m2;

                if(increase>best_increase){
                    best_comm = n_cols[i];
                    best_increase =increase;
                    best_nblinks = n_weights[n_cols[i]];
                }
            }


            G->tot[best_comm] += degrees[node];
            G->in[best_comm] += 2*best_nblinks + G->self_loops[node];
            G->community[node] = best_comm;
            colour[node] = best_comm;

           // unsigned int i = G->offset[node];
           // do{
           //     
           //     double totc = tot_degrees[neigh_col];
           //     double dnc = n_weights[neigh_col];
           //     double degc = w_deg;
           //     double m2 = 2 * (double)G->m;
           //     std::cout << m2 << "\n";



           //     double increase = dnc - totc*degc/m2;

           //     if(increase>best_increase){
           //         best_comm = neigh_col;
           //         best_increase = increase;
           //     }

           //     if(i==G->offset[node+1])break;
           //     neigh_col = colour[G->edges[i]];
           //     if(neigh_col==-1)neigh_col=node_comm;

           //     i++;

           // }while(i<=G->offset[node+1]);
            

            
            //std::cout << node << " " << best_comm << " " << best_increase << std::endl;
            
            if(best_comm!=node_comm)moves++;
            if(moves>0)improvement=true;
            
            
        }
        new_mod = modularity(G, colour);
    } while(moves>0 and new_mod-cur_mod>thresh1);
    return improvement;
}



int main(int argc, char *argv[]){
    /*char input_file[] = "../data/csr/karate_club.csr";*/
    /*graph * G = input(input_file);*/
    graph * G = input(argv[1]);

    std::cout << G->n_v << " " << G->n_e << " " << G->m << std::endl;
    std::cout << "Input taken\n";
    std::vector<int> colouring(G->n_v);
    for(int i=0;i<G->n_v;i++)colouring[i] = i;  
//  G->display();
    std::cout << modularity(G, colouring) << std::endl;
    for(int i=0;i<5;i++){
        std::vector<int> colour(G->n_v);
        std::vector<int> tot_degree(G->n_v);
        std::vector<int> degree(G->n_v);


        initialize(G, colour, tot_degree, degree);
        std::cout << "vectors initialized\n";

        //for(int i=0;i<G->n_v;i++)colour[i]=i; 
        //double moving_phase_delta = moving_phase(G, colour, tot_degree, degree);    
        one_level(G, colour, tot_degree, degree);
        std::cout << "one_level over\n";
        //for(int i=0;i<G->n_v;i++)std::cout << i << " " << colour[i] << std::endl; 
        
       // for(int j=0;j<G->n_v;j++){
       //     colouring[j] = colour[colouring[j]];
       // }

        graph* G_ = aggregate_phase(G, colour, tot_degree, degree, colouring);
        std::cout << "aggregate_phase over\n";



        //G_->display();

        G_->m = G->m;
        G = G_;
        std::cout << "graph copied\n";
        //G->display();


        /*std::set<int> st;*/
        /*for(auto i:colouring)st.insert(i);*/
        //std::cout << st.size() << std::endl;

        std::cout << "--------------" << std::endl;
    }




    //    G_->display();


}
