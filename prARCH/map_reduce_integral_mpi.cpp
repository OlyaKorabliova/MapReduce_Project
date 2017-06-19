#include <iostream>
#include <map>
#include <vector>
#include <math.h>
#include <assert.h>
#include <thread>
#include <mutex>
#include <deque>
#include <condition_variable>
#include <fstream>

#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <map>
#include <string>
#include <sstream>
#include <unistd.h>
//#include "timing.cpp"


using namespace std;
mutex myMutex;
condition_variable cv;
double num = 0;

double commsize = 0, rankk = 0, len = 0;
char procname[MPI_MAX_PROCESSOR_NAME];

double func_calculation(int m, double x1, double x2) {
    double sum1 = 0;
    double sum2 = 0;
    double g;
    for (int i = 1; i <= m; ++i) {
        sum1 += i * cos((i + 1) * x1 + 1);
        sum2 += i * cos((i + 1) * x2 + 1);
    }
    g = -sum1 * sum2;
    return g;
}

double thread_integration(double &r, const double &sum) {
    r += sum;
    return r;
}


template<class Dd, class RF>
Dd reducer_univ(int num_of_threads, deque<Dd> &dm, RF (*fn2)(Dd &, const Dd &), int num) {
    unique_lock<mutex> uniqueLock(myMutex);
    while (dm.size() > 1) {
        Dd map1 = dm.front();
        dm.pop_front();
        Dd map2 = dm.front();
        dm.pop_front();
        uniqueLock.unlock();
        Dd map3 = (*fn2)(map1, map2);
        uniqueLock.lock();
        dm.push_back(map3);
    }
    if (dm.size() != 1 && num != num_of_threads) {
        cv.wait(uniqueLock);
    }

}


using point_t = pair<double, double>;


class Iter2Ddouble : public iterator<random_access_iterator_tag, point_t> {
    double x0, x1, y0, y1, m, pr;
    double x, y;
public:
    double get_pr() const { return pr; }

    Iter2Ddouble(const vector<double> &data) {
        x0 = data[0];
        x1 = data[1];
        y0 = data[2];
        y1 = data[3];
        m = data[4];
        pr = data[5];

        x = x0;
        y = y0;
    }

    Iter2Ddouble end() const {
        Iter2Ddouble t(*this);
        t.x = x1;
        t.y = y1;
        return t;
    }

    point_t operator*() const {
        return make_pair(x + pr / 2.0, y + pr / 2.0);
    }

    Iter2Ddouble &operator++() {
        x += pr;
        if (x > x1) {
            x = x0;
            y += pr;
        }
        return *this;
    }

    Iter2Ddouble &operator--() {
        x -= pr;
        if (x < x0) {
            x = x1;
            y -= pr;
        }
        return *this;
    }

    Iter2Ddouble &operator+=(size_t n) {
        for (size_t i = 0; i < n; ++i)
            ++(*this);
        return *this;
    }

    Iter2Ddouble &operator-=(size_t n) {
        for (size_t i = 0; i < n; ++i)
            --(*this);
        return *this;
    }

    ptrdiff_t operator-=(Iter2Ddouble itr) {
        ptrdiff_t diff = 0;
        if (*this < itr) {
            while (*this < itr) {
                --itr;
                --diff;
            }
        } else {
            while (!(*this < itr)) {
                ++itr;
                ++diff;
            }

        }

        return diff;
    }

    friend bool operator<(const Iter2Ddouble &itr1, const Iter2Ddouble &itr2);

    friend bool operator==(const Iter2Ddouble &itr1, const Iter2Ddouble &itr2);
};

Iter2Ddouble operator+(Iter2Ddouble itr, size_t n) {
    return itr += n;
}

Iter2Ddouble operator-(Iter2Ddouble itr, size_t n) {
    return itr -= n;
}

ptrdiff_t operator-(Iter2Ddouble itr_l, const Iter2Ddouble &itr_r) {
    return itr_l -= itr_r;
}

bool operator<(const Iter2Ddouble &itr1, const Iter2Ddouble &itr2) {
    if (itr1.y == itr2.y) //! Bad code -- use abs|x - y | < eps
    {
        return itr1.x < itr2.x;
    }
    return itr1.y < itr2.y;
}

bool operator==(const Iter2Ddouble &itr1, const Iter2Ddouble &itr2) {
    return itr1.y == itr2.y && itr1.x == itr2.x;  //! Bad code -- use abs|x - y | < eps
}

struct func_wrapper {
    double n;
    double num;

    func_wrapper(double n_, double num_) : n(n_), num(num_) {};

    double operator()(double x1, double x2, vector<double> v, double n, double &num) {
        double r = 0;
        r += func_calculation(v[0], x1, x2) * v[1] * v[1]; // v[1] = pr, v[0] = m
        {
            lock_guard<mutex> lg(myMutex);
            num += 1;
        }

        cv.notify_one();
        cv.notify_all();
        return r;
    }

};


template<class MF, class RF, class I, class N, class Dd>
Dd func_tmpl(I beg, I fin, vector<Dd> vec, MF fn1, RF fn2, N num_of_threads) {
    cout << "~ " << beg << " " << fin << " " << vec[0] << " " << vec[1] << " " << num_of_threads << endl;
    Dd res;
    int l = (int) num_of_threads;
    size_t delta = (fin - beg) / num_of_threads;
    cout << "del:= " << delta << endl;
    thread myThreads[l];
    double result = 0;
    double from_to[] = {beg, fin, num_of_threads, num};
    if (rankk == 0) {
        double from_to[] = {beg, beg + delta, num_of_threads, num};
        result = fn1(from_to[0], from_to[1], vec, from_to[2], ref(from_to[3]));

        for (int j = 1; j < commsize; ++j) {
            MPI_Send(from_to, 4, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            cout << "from_to[0]= " << from_to[0] << endl;
            cout << "from_to[1]= " << from_to[1] << endl;
            from_to[0] = from_to[1];
            from_to[1] += delta;
        }

        for (int i = 1; i < commsize; ++i) {
            double received;
            MPI_Recv(&received, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            result += received;
        }
        cout << "Result: " << result << endl;
    } else {
        MPI_Recv(from_to, 4, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        result = fn1(from_to[0], from_to[1], vec, from_to[2], ref(from_to[3]));
        MPI_Send(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

    }
//    cout << "Result: " << result << endl;

    MPI_Finalize();
    return 0;

}


int main(int argc, char *argv[]) {
    //auto stage1_start_time = get_current_time_fenced();

    int commsize1, rank1, len1;
    char procname[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize1);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank1);
    MPI_Get_processor_name(procname, &len1);
    rankk = rank1;
    len = len1;
    commsize = commsize1;

    vector<double> data = {0, 4, 5, 0.001};
    double n = 2;       // num of threads
    vector<double> v = {data[2], data[3]};
    cout << "Integral" << endl;
    for (auto i : data) cout << "\tdata: " << " " << i << endl;
    cout << func_tmpl(data[0], data[1], v, func_wrapper(n, num), thread_integration, n) << endl;

}
