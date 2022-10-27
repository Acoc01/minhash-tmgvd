#include <vector>
#include <iostream>
#include <limits.h>
#include <functional>

void init(std::vector<std::vector<int>> &adj, int v){
    for(int i = 0; i < v; ++i){
        std::vector<int> n;
        adj.push_back(n);
    }
}

int main(int argc, char* argv[]){
    std::vector<std::vector<int>> adj;

    int v,e;
    std::cin>>v>>e;

    /*Making the adjacency lists*/
    init(adj,v);
    //std::cout<<adj.size()<<'\n'; 
    int a,b;

    for(int i = 0; i < e; ++i){
        std::cin>>a>>b;
        adj[a].push_back(b);
        std::cout<<"insertado\n";
    }

    for(int i = 0; i < v; ++i){
        std::cout<<"a: "<<i<<'\n';
        for(int j = 0; j < adj[i].size(); ++j){
            std::cout<<adj[i][j]<<' ';
        }std::cout<<'\n';
    }
    //return 0;
    /*Minhashing*/
    int arr[v][3];
    int P = 3;//Number of hash functions
    for(int i = 0; i < P; ++i){//For every hash function in P
        std::cout<<"Calculando Minhash\n";
        std::hash<int> ph; //Init of hash function
        /*v-3 = v-k+1, 2-shingle*/
        for(int j = 0; j < v; ++j){//For every adjacency list
        /*abcdef => k = 2 => ab,bc,cd,de,ef => */
        long long int c_j = INT_MAX; // Set inf for the signatures[v_i][i]
            if(adj[j].size() > 2)
                for(int k = 0; k < adj[j].size() - 1; ++k){//Get all the 2-shingles
                    long long int shingle = adj[j][k] + adj[j][k+1];
                    long long int h = ph(shingle);//Hash of the 2-shingle
                    c_j = std::min(h,c_j);//Choose the min between the new hash and the older one
                }
            std::cout<<"llega"<<std::endl;
            arr[j][i] = c_j;
            std::cout<<"signature"<<std::endl;;//Push the final value to the signatures matrix
        }
    }
    std::cout<<"Finalizado\n";
    
    for(int i=0;i<v;i++){
        for(int j=0;j<P;j++){
            std::cout<<arr[i][j]<<"    ";
        }
        std::cout<<std::endl;
    }
    
    return 0;
}