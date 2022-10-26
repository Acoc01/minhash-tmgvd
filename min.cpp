#include <vector>
#include <iostream>
#include <limits.h>
#include <functional>

void init(std::vector<std::vector<int>> *adj, int v){
    for(int i = 0; i < v; ++i){
        std::vector<int> n;
        adj.push_back(n);
    }
}

int main(int argc, char* argv[]){
    std::vector<std::vector<int>> adj;

    int v;
    std::cin>>v;

    /*Making the adjacency lists*/
    init(adj,v);
    
    int a,b;

    while(std::cin>>a>>b){
        adj[a].push_back(b);
    }

    /*Minhashing*/
    std::vector<std::vector<int>> signatures;//Signature Matrix 
    int P = 3;//Number of hash functions
    for(int i = 0; i < P; ++i){//For every hash function in P
        std::hash<int> ph; //Init of hash function
        long long int c_j = INT_MAX; // Set inf for the signatures[v_i][i]
        /*v-3 = v-k+1, 2-shingle*/
        for(int j = 0; j < v; ++j){//For every adjacency list
            for(int k = 0; k < adj[j].size() - 3; ++k){//Get all the 2-shingles
                long long int shingle = adj[j][k] + adj[j][k+1];
                long long int h = ph(shingle);//Hash of the 2-shingle
                c_j = std::min(shingle,c_j);//Choose the min between the new hash and the older one
            }
            signatures[j].push_back(c_j);//Push the final value to the signatures matrix
        }
    }

    return 0;
}