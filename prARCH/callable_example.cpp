#include <iostream>
#include <map>
#include <vector>

//#include <mutex>


using namespace std;

map<string, int> counting_words_worker(const vector<string>& words){
    map<string, int> localm;
    for (size_t  i = 0; i != words.size(); ++i) {
            ++localm[words[i]];
    }
    return localm;
}

map<string, int> reduce_f(map<string, int> m, const map<string, int>& localm) {
    for (auto it = localm.begin(); it != localm.end(); ++it)
        m[it->first] += it->second;

    return m;
}

template <class V, class MF, class RF>
auto func_tmpl(const V& v, MF fn1, RF fn2) ->
            decltype( fn2(fn1(v[0]), fn1(v[0])) ) {
    // V - масив блоків (векторів)
    decltype(func_tmpl(v, fn1, fn2)) res;        // big map<string, int>
    for (size_t i = 0; i < v.size(); ++i) {
        auto x = fn1(v[i]);        //return small   map<string, int>
        res = fn2(move(res), x);
    }
    for (auto i : res){
        cout << i.first << " - " << i.second << endl;
    }
    return res;
}

void callE(const vector<vector<string>> &words){
    func_tmpl(words, counting_words_worker, reduce_f);
}

int main()
{
    vector<string> v1 = {"asdf", "wert", "aaa", "qwerty"};
    vector<string> v2 = {"asdf", "poiuy", "qwewrtdxfb", "qwerty"};
    vector<vector<string>> v3;
    v3.push_back(v1);
    v3.push_back(v2);
    callE(v3);

}

// template <class T>  T reduce_f(T m, const T& localm)
// func_tmpl(words, counting_words_worker, reduce_f<map<string, int>>);