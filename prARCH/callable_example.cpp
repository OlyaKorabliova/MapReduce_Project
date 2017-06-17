#include <iostream>
#include <map>
#include <vector>
#include <math.h>
#include <assert.h>
#include <thread>
#include <mutex>
#include <deque>
#include <condition_variable>


using namespace std;
mutex myMutex;
mutex mx;
condition_variable cv;

deque <map<string, int>> dm;
int num = 0;

void printMap(const map<string, int> &m) {
    for (auto &elem : m) {
        cout << elem.first << " " << elem.second << "\n";
    }
}


map<string, int> reduce_f( map<string, int> m, const map<string, int>& localm ) {

    for (auto it = localm.begin(); it != localm.end(); ++it)
        m[it->first] += it->second;
    return m;
}

deque<map<string, int>> reducer(int num_of_threads){
    unique_lock<mutex> uniqueLock(myMutex);
    while (dm.size() > 1) {
        map<string, int> map1 = dm.front();
        dm.pop_front();
        map<string, int> map2 = dm.front();
        dm.pop_front();
        uniqueLock.unlock();
        map<string, int> map3 = reduce_f(map1, map2);

        uniqueLock.lock();
        dm.push_back(map3);
    }if (dm.size() == 1 && num == num_of_threads) {
        printMap(dm.front());
    } else {
        cv.wait(uniqueLock);
    }
    return dm;
}

map<string, int> counting_words_worker( vector<string>::const_iterator beg, vector<string>::const_iterator  fin ){
    map<string, int> localm;
    for (auto i = beg; i != fin; ++i) {
        ++localm[*i];
    }
    {
        lock_guard<mutex> lg(myMutex);
        dm.push_back(localm);
    }
    cv.notify_one();
    cv.notify_all();
    num += 1;
    return localm;
}



template <class MF, class RF, class I, class N>
auto func_tmpl(I beg, I fin, MF fn1, RF fn2, N num_of_threads) -> decltype( fn2(num_of_threads) ) {
    decltype(func_tmpl(beg, fin, fn1, fn2, num_of_threads)) res;        // big map<string, int>

    size_t delta = (fin - beg) / num_of_threads;
    vector<I> segments;
    for (auto i = beg; i < fin - delta; i += delta) {
        segments.push_back(i);
    }
    segments.push_back(fin);

    thread myThreads [num_of_threads];
    for (auto i = segments.begin(); i < segments.end()-1; ++i) {
        myThreads[i - segments.begin()] = thread(fn1, *i, *(i+1));
        fn2(num_of_threads);
    }

    for (auto& th : myThreads) th.join();
    res = fn2(num_of_threads);
    return res;
}


int main()
{
    vector<string> v = {"aaaa", "bbbbb", "ccccc", "ddddd", "aaaa", "eeeee", "wwwwwwwww", "ccccc", "ccccc", "ccccc"};
    cout<<"Words counter" << endl;
    int num_of_threads = 3;
    func_tmpl(v.begin(), v.end(), counting_words_worker, reducer, num_of_threads);
}
