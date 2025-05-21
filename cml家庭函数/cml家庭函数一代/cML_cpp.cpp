// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // RcppArmadillo 核心头文件，提供了 Armadillo 库的接口
#include <vector>          // C++ 标准库：动态数组 (std::vector)
#include <algorithm>       // C++ 标准库：包含排序等算法 (std::sort)
#include <cmath>           // C++ 标准库：数学函数 (std::pow, std::abs)

struct additem
{
    double value1;
    double value2;
    double calculateSum() const { return value1 + value2; };
};

struct SortItem
{
    double value;      // 要排序的值
    arma::uword index; // 原始索引 (arma::uword 是 Armadillo 使用的无符号整数类型，通常用于索引)

    // 重载小于操作符 (<)，用于自定义排序规则
    // 注意：这里定义为 value > other.value，所以 std::sort 默认会进行降序排序
    bool operator>(const SortItem &other) const
    {
        return value > other.value; // 按 value 降序排序
    }
};
// [[Rcpp::export]]
void cMl_para_test(
    const arma::mat &beta_y_hat)
{
    SortItem test;
    test.value = beta_y_hat(0);
    test.index = 0;
    if (test.index == 0)
    {
        Rcpp::Rcout << "第一个元素是" << test.value << std::endl;
    }
    return;
}