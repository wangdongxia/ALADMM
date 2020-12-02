#ifndef SPARSEALLREDUCE_REDUCE_OPERATOR_H
#define SPARSEALLREDUCE_REDUCE_OPERATOR_H

//封装了一些Reduce操作，主要实现了累加，累乘，最大最小操作，同样采用模板和模板特化，最大化程序性能
//由于模板函数不能偏特化，因此需要对各个数值类型分别实现，正在考虑更简洁的实现方式
namespace spar {

struct SumOperator {};

struct MinOperator {};

struct MaxOperator {};

struct ProductOperator {};

template<class O, class T>
void Reduce(T &a, const T &b);

template<>
inline void Reduce<SumOperator>(double &a, const double &b) {
    a += b;
}

template<>
inline void Reduce<MinOperator>(double &a, const double &b) {
    if (b < a) {
        a = b;
    }
}

template<>
inline void Reduce<MaxOperator>(double &a, const double &b) {
    if (b > a) {
        a = b;
    }
}

template<>
inline void Reduce<ProductOperator>(double &a, const double &b) {
    a *= b;
}

template<>
inline void Reduce<SumOperator>(float &a, const float &b) {
    a += b;
}

template<>
inline void Reduce<MinOperator>(float &a, const float &b) {
    if (b < a) {
        a = b;
    }
}

template<>
inline void Reduce<MaxOperator>(float &a, const float &b) {
    if (b > a) {
        a = b;
    }
}

template<>
inline void Reduce<ProductOperator>(float &a, const float &b) {
    a *= b;
}

template<>
inline void Reduce<SumOperator>(int &a, const int &b) {
    a += b;
}

template<>
inline void Reduce<MinOperator>(int &a, const int &b) {
    if (b < a) {
        a = b;
    }
}

template<>
inline void Reduce<MaxOperator>(int &a, const int &b) {
    if (b > a) {
        a = b;
    }
}

template<>
inline void Reduce<ProductOperator>(int &a, const int &b) {
    a *= b;
}

template<>
inline void Reduce<SumOperator>(long &a, const long &b) {
    a += b;
}

template<>
inline void Reduce<MinOperator>(long &a, const long &b) {
    if (b < a) {
        a = b;
    }
}

template<>
inline void Reduce<MaxOperator>(long &a, const long &b) {
    if (b > a) {
        a = b;
    }
}

template<>
inline void Reduce<ProductOperator>(long &a, const long &b) {
    a *= b;
}
}

#endif //SPARSEALLREDUCE_REDUCE_OPERATOR_H
