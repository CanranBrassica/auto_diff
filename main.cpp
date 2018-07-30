#include <iostream>
#include <sstream>
#include <functional>
#include <type_traits>
#include <tuple>


namespace kyadet
{
// Argクラスのインスタンスを用いて式を構築
// この式をFuncクラスのコンストラクタに渡し、
// 式の後に使用したArgインスタンスを渡すと、
// その順番で引数をとる関数が生成される。operator()で呼び出し
// また、Funcクラスはdiffメンバメソッドを持ち、
// ここにArgのインスタンスを渡すとそれでを偏微分した
// Funcインスタンスが生成される

// Argクラス
struct Arg
{
    char id_;
    constexpr Arg(char id) : id_(id){};
    constexpr auto func() const
    {
        return [id = this->id_](auto... args) { return id_match_num()(id, args...); };
    }

    friend std::ostream& operator<<(std::ostream& os, const Arg& arg) { return os << arg.id_; }

private:
    template <class... Args>
    constexpr static auto id_match_num() {
        return [](char id, Args &&... args) {
            static_assert(sizeof...(args) != 0, "invalid argument");
            for (auto const &arg : {args...}) {
                if (id == arg.first) {
                    return arg.second;
                }
            }
            //throw std::invalid_argument{"cannot be called"};
            return 0.0;  // cannot be called
        };
    }
};

// 二項演算子クラス
template <class Op, class L, class R>
struct BinOp
{
    L l_;
    R r_;
    constexpr BinOp(L const& l, R const& r) noexcept : l_(l), r_(r) {}

    template <class... IDs>
    constexpr auto func(IDs... ids) const
    {
        return Op::func(l_.func(ids...), r_.func(ids...));
    }

    friend std::ostream& operator<<(std::ostream& os, BinOp<Op, L, R> const& bin)
    {
        return Op::set_os(os, bin.l_, bin.r_);
    }
};

struct Add
{
    template <class L, class R>
    static std::ostream& set_os(std::ostream& os, L const& l, R const& r)
    {
        return os << l << "+" << r;
    }

    template <class LFunc, class RFunc>
    constexpr static auto func(LFunc const& l_func, RFunc const& r_func)
    {
        return [l_func, r_func](auto... args) { return l_func(args...) + r_func(args...); };
    }
};

template <class L, class R>
constexpr auto operator+(L const& l, R const& r) noexcept
{
    return BinOp<Add, L, R>{l, r};
}

template <class T>
struct Func
{
    T e_;
    constexpr Func(T const& expr) : e_(expr) {}

    template <class... Args>
    constexpr auto create(Args...)
    {}
};

}  // namespace kyadet

#include "../demangle/demangle.hpp"

int main()
{
    using namespace kyadet;

    constexpr auto x = Arg{'x'};
    constexpr auto y = Arg{'y'};
    constexpr auto z = x + y + x;

    constexpr auto g = z.func();

    static_assert(g(std::pair<char, double>('x', 1.0), std::pair<char, double>('y', 4.0)) == 6, "");

    /*
    constexpr g = z.func('x', 'y');
    g(1, 4);
    g(1)(4);
    */

    /*
    constexpr auto x = Arg{};
    constexpr auto y = Arg{};
    constexpr auto z = x - y;
    constexpr auto sub = Func{z.func(), x, y};
    static_assert(sub(1, 2) == 1 - 2, "λx.λy.(x - y)");
    constexpr auto sub2 = Func{z.func(), y, x};
    static_assert(sub(1, 2) == 2 - 1, "λy.λx.(x - y)");
     */

    return 0;
}
