#include<bits/stdc++.h>
using namespace std;

int main(int argc, char* argv[]){
    // Creating a filestream for input file
    fstream fin;
    fin.open(argv[1], ios::in);

    fstream fout;
    fout.open(argv[2], ios::out);
    
    // Parsing the file and populating the graph
    map<int, set<pair<int,int>>> graph;
    string s;
    int num_vert=0, num_edge=0;
    while(getline(fin, s)){
        bool hasWeight = false;
        stringstream str(s);    
        string s1, s2, s3;
        getline(str, s1, ',');
        getline(str, s2, ',');
        if(getline(str, s3, ','))hasWeight = true;
        
        bool isEdge = s1.size()>0 and isdigit(s1[0]);
        // Checking if the string represents an integer or not
        for(auto c:s1){
            if(!isdigit(c)){isEdge=false;break;}
        }
        for(auto c:s2){
            if(!isdigit(c)){isEdge=false;break;}
        }
        if(isEdge==false)continue;

        // Adding the edge
        int a = stoi(s1);
        int b = stoi(s2);
        num_vert = max(num_vert, max(a, b));
        num_edge-=graph[a].size()+graph[b].size();
        graph[a].insert(make_pair(b, hasWeight?stoi(s3):1));
        graph[b].insert(make_pair(a, hasWeight?stoi(s3):1));// Comment this line if the graph is directed
        num_edge+=graph[a].size()+graph[b].size();
    }
    num_vert+=1;

    vector<int> offset(num_vert+1);
    vector<int> edges(num_edge);
    vector<int> weights(num_edge);

    for(int i=0;i<num_vert;i++){
        int j=offset[i];
        for(pair<int, int> p : graph[i]){
            edges[j] = p.first;
            weights[j] = p.second;
            j++;
        }
        offset[i+1]=j;
    }

    fout << num_vert << " " << num_edge << endl;
    for(int i=0;i<num_vert+1;i++)fout << offset[i] << " ";
    fout << endl;
    for(int i=0;i<num_edge;i++)fout << edges[i] << " ";
    fout << endl;
    for(int i=0;i<num_edge;i++)fout << weights[i] << " ";
    fout << endl;

}
