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
condition_variable cv;

int num = 0;


double func_calculation(int m, double x1, double x2) {
    double sum1 = 0;
    double sum2 = 0;
    double g;
    for (int i = 1; i <= m; ++i)
    {
        sum1 += i * cos((i + 1) * x1 + 1);
        sum2 += i * cos((i + 1) * x2 + 1);
    }

    g = - sum1 * sum2;
    return g;
}

using point_t = pair<double, double>;


class Iter2Ddouble: public iterator<random_access_iterator_tag, point_t >
{
    double x0, x1, y0, y1, m, pr;
    double x, y;
public:
    double get_pr() const {return pr; }
    Iter2Ddouble(const vector<double>& data)
    {
        x0 = data[0];
        x1  = data[1];
        y0 = data[2];
        y1  = data[3];
        m = data[4];
        pr = data[5];

        x = x0;
        y = y0;
    }

    Iter2Ddouble end() const{
        Iter2Ddouble t(*this);
        t.x = x1;
        t.y = y1;
        return t;
    }

    point_t operator*() const
    {
        return make_pair(x + pr / 2.0, y + pr / 2.0);
    }

    Iter2Ddouble& operator++()
    {
        x += pr;
        if(x>x1)
        {
            x = x0;
            y += pr;
        }
        return *this;
    }

    Iter2Ddouble& operator--()
    {
        x -= pr;
        if(x<x0)
        {
            x = x1;
            y -= pr;
        }
        return *this;
    }

    Iter2Ddouble& operator+=(size_t n)
    {
        //! TODO: write without loop!
        for(size_t i = 0; i<n; ++i)
            ++(*this);
        return *this;
    }
    Iter2Ddouble& operator-=(size_t n)
    {
        //! TODO: write without loop!
        for(size_t i = 0; i<n; ++i)
            --(*this);
        return *this;
    }
    ptrdiff_t operator-=(Iter2Ddouble itr)
    {
        //! TODO: write without loop!
        ptrdiff_t diff = 0;
        if( *this < itr) {
            while (*this < itr) {
                --itr;
                --diff;
            }
        }else{
            while ( !(*this < itr) ) {
                ++itr;
                ++diff;
            }

        }

        return diff;
    }

    friend bool operator<(const Iter2Ddouble& itr1, const Iter2Ddouble& itr2);
    friend bool operator==(const Iter2Ddouble& itr1, const Iter2Ddouble& itr2);
};

Iter2Ddouble operator+(Iter2Ddouble itr, size_t n)
{
    return itr += n;
}
Iter2Ddouble operator-(Iter2Ddouble itr, size_t n)
{
    return itr -= n;
}

ptrdiff_t  operator-(Iter2Ddouble itr_l, const Iter2Ddouble& itr_r)
{
    return itr_l -= itr_r;
}

bool operator<(const Iter2Ddouble& itr1, const Iter2Ddouble& itr2)
{
    if( itr1.y == itr2.y ) //! Bad code -- use abs|x - y | < eps
    {
        return itr1.x < itr2.x;
    }
    return itr1.y < itr2.y;
}

bool operator==(const Iter2Ddouble& itr1, const Iter2Ddouble& itr2)
{
    return itr1.y == itr2.y && itr1.x == itr2.x;  //! Bad code -- use abs|x - y | < eps
}
struct func_wrapper
{
    int m;
    deque<double> d;
    func_wrapper(int m_, deque<double> d_): m(m_), d(d_){};
//    double operator()(const point_t& p)
//    {
//        return func_calculation(m, p.first, p.second);
//    }
    double operator()(const Iter2Ddouble& itr1, const Iter2Ddouble& itr2,  deque<double>& d)
    {
        double r = 0;
        for(auto i = itr1; i <itr2; ++i) {
            r += func_calculation(m, (*i).first, (*i).second) * itr1.get_pr() * itr1.get_pr();
            // cout << (*i).first << " " <<  (*i).second << " " << r << endl;
        }
        {
            lock_guard<mutex> lg(myMutex);
            d.push_back(r);
        }
        cv.notify_one();
        cv.notify_all();
        num += 1;
        return r;
    }

};


double thread_integration(double r, const double &sum) {
    r += sum;
    return r;
}


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


template<class RF, class Dd>
Dd reducer_univ(int num_of_threads, deque<Dd> &dm, RF (*fn2)(Dd, const Dd&)) {
    unique_lock<mutex> uniqueLock(myMutex);
    while (dm.size() > 1) {
        Dd d1 = dm.front();
        dm.pop_front();
        Dd d2 = dm.front();
        dm.pop_front();
        uniqueLock.unlock();
        Dd d3 = (*fn2)(d1, d2);

        uniqueLock.lock();
        dm.push_back(d3);
        cout << "1" << endl;
    }
    if (dm.size() != 1 && num != num_of_threads){
        cout << "3" << endl;
        cv.wait(uniqueLock);
    }
    return dm.front();
}

map<string, int> counting_words_worker( vector<string>::const_iterator beg, vector<string>::const_iterator  fin, deque<map<string, int>> &dm ){
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



template <class MF, class RF, class I, class N, class Dd>
Dd func_tmpl(I beg, I fin, MF fn1, deque<Dd> &d, RF fn2, N num_of_threads) {

    Dd res;
    size_t delta = (fin - beg) / num_of_threads;
    vector<I> segments;
    for (auto i = beg; i < fin - delta; i += delta) {
        segments.push_back(i);
    }
    segments.push_back(fin);

    thread myThreads [num_of_threads];
    for (auto i = segments.begin(); i < segments.end()-1; ++i) {
        myThreads[i - segments.begin()] = thread(fn1, *i, *(i+1), ref(d));
        reducer_univ(num_of_threads, d, fn2);
    }

    for (auto& th : myThreads) th.join();
    res = reducer_univ(num_of_threads, d, fn2);
    return res;
}


int main()
{
//    deque <map<string, int>> dm;
//    vector<string> v = {"aaaa", "bbbbb", "ccccc", "ddddd", "aaaa", "eeeee", "wwwwwwwww", "ccccc", "ccccc", "ccccc"};
//    cout<< "Words counter" << endl;
//    int num_of_threads = 3;
//    map<string, int> m = func_tmpl(v.begin(), v.end(), counting_words_worker, dm, reduce_f, num_of_threads);
//    printMap(m);

    vector<double> data = {0, 1, 0, 1, 5, 0.001};
    deque<double> d;
    Iter2Ddouble itr(data);
    cout << "Integral" << endl;
    cout << func_tmpl(itr, itr.end(), func_wrapper(5, d), d, thread_integration, 4) << endl;

}
