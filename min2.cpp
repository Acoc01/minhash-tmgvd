#include <vector>
#include <iostream>
#include <limits.h>
#include <functional>
#include <string>
#include <fstream>
#include <sstream>


unsigned hash_str(std::string s,int F, int A, int B)
{
   unsigned h = F;
   for(int i = 0; i < s.size(); ++i){
     h = (h * A) ^ (s[i] * B);
   }
   return h; // or return h % C;
}

int main(int argc, char* argv[]){
    std::ifstream myfile; myfile.open("input");

    long long int v,a;
    std::string list;
    long long int l = 0;
    int P = 3;
    myfile>>v;
    std::vector<std::vector<long long int>> arr;
    std::cout<<arr.size()<<'\n';
    while(getline(myfile, list))
    {
        //std::cout<<list<<'\n';
        long long int d,c_j = INT_MAX;
        std::vector<long long int> row;
        std::stringstream lineStream(list);

        while (lineStream >> d)
            row.push_back(d);

        std::vector<std::string>ladj;
        for(int i = 0; i < row.size(); ++i){ 
            std::string aux = std::to_string(row[i]);
            ladj.push_back(aux);
            //std::cout<<aux<<'\n';
        }//std::cout<<'\n';
        std::vector<long long int> hashes;
        for(int i = 0; i < P; ++i){
            if(ladj.size()>=2){
                for(int j = 0; j < ladj.size()-1; ++j){
                    std::string shingle;
                    //std::cout<<ladj[j]<<' '<<ladj[j+1]<<'\n';
                    shingle.append(ladj[j]);
                    shingle.append(ladj[j+1]);
                    long long int h;
                    if(i == 0) h = hash_str(shingle,59,7,13);
                    if(i == 1) h = hash_str(shingle,81,7,17);
                    if(i == 2) h = hash_str(shingle,83,7,27);//Hash of the 2-shingle
                    c_j = std::min(h,c_j);//Choose the min between the new hash and the older one
                }
            }
            hashes.push_back(c_j);
        }
        arr.push_back(hashes);
        l++;
    }
    std::cout<<l<<'\n';
    for(int i=0;i<arr.size();i++){
        for(int j=0;j<P;j++){
            std::cout<<arr[i][j]<<"    ";
        }
        std::cout<<std::endl;
    }
    return 0;
}