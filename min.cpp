#include <vector>
#include <iostream>
#include <limits.h>
#include <functional>
#include <string>
#include <time.h>

unsigned hash_str(std::string s,int F, int A, int B)
{
   unsigned h = F;
   for(int i = 0; i < s.size(); ++i){
     h = (h * A) ^ (s[i] * B);
   }
   return h; // or return h % C;
}

void init(std::vector<std::vector<std::string>> &adj, int v){
    for(int i = 0; i < v; ++i){
        std::vector<std::string> n;
        adj.push_back(n);
    }
}

int main(int argc, char* argv[]){
    srand(time(NULL));
    std::vector<std::vector<std::string>> adj;

    int v,e;
    std::cin>>v>>e;

    /*Making the adjacency lists*/
    init(adj,v);
    //std::cout<<adj.size()<<'\n'; 
    int a;
    std::string b;

    for(int i = 0; i < e; ++i){
        std::cin>>a;
        getline(std::cin,b);
        adj[a].push_back(b);
        std::cout<<"insertado\n";
    }

    for(int i = 0; i < v; ++i){
        //std::cout<<"a: "<<i<<'\n';
        if(adj[i].size() != 0){
            for(int j = 0; j < adj[i].size(); ++j){
                std::cout<<adj[i][j]<<' ';
            }std::cout<<'\n';
        }
    }
    //return 0;
    /*Minhashing*/
    long long int arr[v][3];
    std::hash<std::string> h1,h2,h3;
    int P = 3;//Number of hash functions
    for(int i = 0; i < P; ++i){//For every hash function in P
        std::cout<<"Calculando Minhash\n";
        /*v-3 = v-k+1, 2-shingle*/
        for(int j = 0; j < v; ++j){//For every adjacency list
        /*abcdef => k = 2 => ab,bc,cd,de,ef => */
        long long int c_j = INT_MAX; // Set inf for the signatures[v_i][i]
            if(adj[j].size() >= 2)
                for(int k = 0; k < adj[j].size() - 1; ++k){//Get all the 2-shingles
                    std::string shingle;
                    shingle.append(adj[j][k]);
                    shingle.append(adj[j][k+1]);
                    std::cout<<shingle<<' ';
                    long long int h;
                    if(i == 0) h = hash_str(shingle,59,7,13);
                    if(i == 1) h = hash_str(shingle,81,59,17);
                    if(i == 2) h = hash_str(shingle,83,17,13);//Hash of the 2-shingle
                    c_j = std::min(h,c_j);//Choose the min between the new hash and the older one
                }std::cout<<'\n';
            std::cout<<"llega"<<std::endl;
            arr[j][i] = c_j;
            std::cout<<"signature"<<std::endl;//Push the final value to the signatures matrix
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