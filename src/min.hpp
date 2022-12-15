#ifndef MINHASH_TEST
#define MINHASH_TEST

#include <fstream>
#include <functional>
#include <iostream>
#include <limits.h>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace minhash{
  std::vector<std::vector<long long int> > graph;

  unsigned hash_str(std::string s, int F, int A, int B) {
    unsigned h = F;
    for (int i = 0; i < s.size(); ++i) {
      h = (h * A) ^ (s[i] * B);
    }
    return h; // or return h % C;
  }

  std::vector<std::vector<int>> min(std::string input) {
    const char * c = input.c_str();
    std::ifstream myfile;
    myfile.open(c);

    std::map<long long int, int> mp;
    long long int v, a;
    std::string list;
    long long int l = 0;
    int P = 2;
    //myfile >> v;
    std::vector<std::vector<long long int> > arr;
    // std::cout << arr.size() << '\n';
    while (getline(myfile, list)) {
      // std::cout<<list<<'\n';
      long long int d, c_j = INT_MAX;
      std::vector<long long int> row;
      std::stringstream lineStream(list);

      while (lineStream >> d)
        row.push_back(d);

      std::vector<std::string> ladj;
      for (int i = 0; i < row.size(); ++i) {
        std::string aux = std::to_string(row[i]);
        ladj.push_back(aux);
        // std::cout<<aux<<'\n';
      } // std::cout<<'\n';
      std::vector<long long int> hashes;
      for (int i = 0; i < P; ++i) {
        if (ladj.size() >= 2) {
          for (int j = 0; j < ladj.size() - 1; ++j) {
            std::string shingle;
            // std::cout<<ladj[j]<<' '<<ladj[j+1]<<'\n';
            shingle.append(ladj[j]);
            shingle.append(ladj[j + 1]);
            long long int h;
            if (i == 0)
              h = hash_str(shingle, 59, 7, 13);
            if (i == 1)
              h = hash_str(shingle, 81, 7, 17);
            // if(i == 2) h = hash_str(shingle,83,7,27);//Hash of the 2-shingle
            c_j = std::min(
                h, c_j); // Choose the min between the new hash and the older one
          }
        }
        hashes.push_back(c_j);
        if (c_j != INT_MAX)
          mp[c_j]++;
      }
      arr.push_back(hashes);
      graph.push_back(row);
      l++;
    }
    /*for (int i = 0; i < arr.size(); i++) {
      for (int j = 0; j < P; j++) {
        std::cout << arr[i][j] << "    ";
      }
      std::cout << std::endl;
    }*/

    std::map<long long int, int>::iterator it = mp.begin();
    std::vector<std::vector<int>> clusters;
    for (it; it != mp.end(); ++it) {
      // std::cout<<it->first<<' '<<it->second<<'\n';
      std::vector<int> aux;
      long long int key = it->first;
      int value = it->second;
      if (value > 1 && key != INT_MAX) {
        for (int i = 0; i < arr.size(); ++i) {
          long long int h1 = arr[i][0], h2 = arr[i][1];
          if (h1 == key || h2 == key) {
            aux.push_back(i);
          }
        }
        clusters.push_back(aux);
      }
    }
    /*for (int i = 0; i < clusters.size(); ++i) {
      std::cout<<"Cluster "<<i+1<<": ";
      for (int j = 0; j < clusters[i].size(); ++j) {
          std::cout<<clusters[i][j]<<" ";
      }std::cout<<'\n';
    }*/

    return clusters;
  }
}
#endif