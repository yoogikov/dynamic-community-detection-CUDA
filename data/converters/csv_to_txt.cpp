#include<bits/stdc++.h>
using namespace std;

int main(int argc, char* argv[]){
    // Creating a filestream for input file
    fstream fin;
    fin.open(argv[1], ios::in);

    fstream fout;
    fout.open(argv[2], ios::out);
    
    // Parsing the file and populating the graph
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

        fout << s1 << " " << s2;
        if(hasWeight)fout << s3 << "\n";
        else fout << " 1" << "\n";
    }
}
