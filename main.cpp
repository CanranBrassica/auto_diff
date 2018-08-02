#include <iostream>
#include <cmath>
#include <type_traits>
#include <tuple>

#define USE_STRING 1

namespace kyadet
{

class Var;

template <size_t N>
struct ExprMgr
{
    std::array<Var, N> v_arr;
    constexpr ExprMgr() : v_arr(init_v_arr()){};

    template <class T>
    constexpr auto create_func(T&& expr) const
    {
        return create_func_impl(std::forward<T>(expr), std::is_scalar<T>{});
    }

private:
    constexpr static auto init_v_arr()
    {
        std::array<Var, N> arr{};
        for (size_t index = 0; index < N; ++index) {
            arr.at(index) = Var{index};
        }
        return arr;
    }

    template <class T>
    constexpr static auto create_func_impl(T&& expr, std::true_type)
    {
        return [v = expr](auto... args) noexcept { return v; };
    }

    template <class T>
    constexpr static auto create_func_impl(T&& expr, std::false_type)
    {
        return [f = expr.func()](auto... args)
        /*
                noexcept(noexcept(
                    expr.func()(
                        std::declval<
                            std::array<
                                std::remove_reference_t<
                                    decltype(std::get<0>(std::declval<std::tuple<decltype(args)...>>()))
                                >,
                                N
                            >
                        >()
                    )
                ))*/
        {
            static_assert(sizeof...(args) == N, "A number of args must be N.");
            return f(std::array<
                std::remove_reference_t<
                    decltype(std::get<0>(std::declval<std::tuple<decltype(args)...>>()))>,
                N>{args...});
        };
    }
};  // class ExprMgr

// Varクラス
struct Var
{
    constexpr bool operator==(Var const& other) const { return index_ == other.index_; }

    constexpr auto diff(Var const& v) const
    {
        if (index_ == v.index_) {
            return 1;
        } else {
            return 0;
        }
    }

    constexpr auto func() const
    {
        return [index = index_](auto const& args) noexcept(noexcept(args[std::declval<size_t>()]))
        {                        // noexcept(operator[]) == true
            return args[index];  // argsから該当する部分だけ返す
        };
    }

private:
    template <size_t N>
    friend class ExprMgr;

    constexpr Var(size_t index) : index_(index) {}
    constexpr Var() = default;
    size_t index_ = 0;

#if USE_STRING
public:
    friend std::ostream& operator<<(std::ostream& os, Var const& v)
    {
        return os << "v[" << v.index_ << "]";
    }
#endif
};  // class Var

template <class Op, class X>
struct MonoOp
{
    X x_;
    constexpr MonoOp(X const& x) : x_(x) {}

    constexpr auto func() const
    {
        return Op::func(x_);  // λx.Op(x)   <- 関数オブジェクト
    }

    constexpr auto diff(Var const& v) const
    {
        return Op::diff(x_, v);  // d.Op(x_)/d.v  <- 式木
    }

};  // class MonoOp

template <class Op, class L, class R>
struct BinOp
{
    L l_;
    R r_;
    constexpr BinOp(L const& l, R const& r) : l_(l), r_(r) {}

    constexpr auto func() const
    {
        return Op::func(l_, r_);  // λx.λy.Op(x, y)   <- 関数オブジェクト
    }

    constexpr auto diff(Var const& v) const
    {
        return Op::diff(l_, r_, v);  // d.Op(l_, r_)/d.v  <- 式木
    }
};  // class BinOp

#if USE_STRING
// operator<< for MonoOp
template <class Op, class X>
constexpr std::ostream& operator<<(std::ostream& os, MonoOp<Op, X> const& v)
{
    return set_ostream_impl(os, v);
}

// operator<< for BinOp
template <class Op, class L, class R>
std::ostream& operator<<(std::ostream& os, BinOp<Op, L, R> const& v)
{
    return set_ostream_impl(os, v);
}
#endif

class Exp;

template <class X>
constexpr auto exp(X const& x)
{
    return MonoOp<Exp, X>{x};
}

struct Exp
{
    template <class X>
    static constexpr auto func(X const& x)
    {
        return [x_func = x.func()](auto const& args) noexcept(noexcept(exp(x.func()(args))))
        {
            using std::exp;
            return exp(x_func(args));
        };
    }

    template <class X>
    static constexpr auto diff(X const& x, Var const& v)
    {  // d.exp(x)/dv = exp(x) * dx/dv
        return exp(x) * x.diff(v);
    }
};

#if USE_STRING
template <class X>
std::ostream& set_ostream_impl(std::ostream& os, MonoOp<Exp, X> const& v)
{
    return os << "exp(" << v.x_ << ")";
}
#endif

struct Add
{
private:
    // Expr, Expr
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::false_type, std::false_type)
    {
        return [ l_func = l.func(), r_func = r.func() ](auto const& args) noexcept(noexcept(l.func()(args)) && noexcept(r.func()(args)))
        {
            return l_func(args) + r_func(args);
        };
    }

    // Scalar, Expr
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::true_type, std::false_type)
    {
        return [ l, r_func = r.func() ](auto const& args) noexcept(noexcept(r.func()(args)))
        {
            return l + r_func(args);
        };
    }

    // Expr, Scalar
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::false_type, std::true_type)
    {
        return [ l_func = l.func(), r ](auto const& args) noexcept(noexcept(l.func()(args)))
        {
            return l_func(args) + r;
        };
    }

public:
    template <class L, class R>
    static constexpr auto func(L const& l, R const& r)
    {
        return func_impl(l, r, std::is_scalar<L>{}, std::is_scalar<R>{});
    }

private:
    // Expr, Expr
    template <class L, class R>
    static constexpr auto diff_impl(L const& l, R const& r, Var const& v, std::false_type, std::false_type)
    {
        return l.diff(v) + r.diff(v);  // d(l+r)/dv = dl/dv + dr/dv
    }

    // Scalar, Expr
    template <class L, class R>
    static constexpr auto diff_impl(L const& l, R const& r, Var const& v, std::true_type, std::false_type)
    {
        return r.diff(v);
    }

    // Expr, Scalar
    template <class L, class R>
    static constexpr auto diff_impl(L const& l, R const& r, Var const& v, std::false_type, std::true_type)
    {
        return l.diff(v);
    }

public:
    template <class L, class R>
    static constexpr auto diff(L const& l, R const& r, Var const& v)
    {
        return diff_impl(l, r, v, std::is_scalar<L>{}, std::is_scalar<R>{});
    }
};

template <class L, class R>
constexpr auto operator+(L const& l, R const& r) noexcept
{
    return BinOp<Add, L, R>{l, r};
}

struct Mul
{
private:
    // Expr, Expr
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::false_type, std::false_type)
    {
        return [ l_func = l.func(), r_func = r.func() ](auto const& args) noexcept(noexcept(l.func()(args)) && noexcept(r.func()(args)))
        {
            return l_func(args) * r_func(args);
        };
    }

    // Scalar, Expr
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::true_type, std::false_type)
    {
        return [ l, r_func = r.func() ](auto const& args) noexcept(noexcept(r.func()(args)))
        {
            return l * r_func(args);
        };
    }

    // Expr, Scalar
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::false_type, std::true_type)
    {
        return [ l_func = l.func(), r ](auto const& args) noexcept(noexcept(l.func()(args)))
        {
            return l_func(args) * r;
        };
    }

public:
    template <class L, class R>
    static constexpr auto func(L const& l, R const& r)
    {
        return func_impl(l, r, std::is_scalar<L>{}, std::is_scalar<R>{});
    }

private:
    // Expr, Expr
    template <class L, class R>
    static constexpr auto diff_impl(L const& l, R const& r, Var const& v, std::false_type, std::false_type)
    {
        return l.diff(v) * r + l * r.diff(v);  // d(l*r)/dv = dl/dv * r + l * dr/dv
    }

    // Scalar, Expr
    template <class L, class R>
    static constexpr auto diff_impl(L const& l, R const& r, Var const& v, std::true_type, std::false_type)
    {
        return l * r.diff(v);
    }

    // Expr, Scalar
    template <class L, class R>
    static constexpr auto diff_impl(L const& l, R const& r, Var const& v, std::false_type, std::true_type)
    {
        return l.diff(v) * r;
    }

public:
    template <class L, class R>
    static constexpr auto diff(L const& l, R const& r, Var const& v)
    {
        return diff_impl(l, r, v, std::is_scalar<L>{}, std::is_scalar<R>{});
    }
};

template <class L, class R>
constexpr auto operator*(L const& l, R const& r) noexcept
{
    return BinOp<Mul, L, R>{l, r};
}

#if USE_STRING
template <class L, class R>
std::ostream& set_ostream_impl(std::ostream& os, BinOp<Add, L, R> const& v)
{
    return os << "(" << v.l_ << "+" << v.r_ << ")";
}

template <class L, class R>
std::ostream& set_ostream_impl(std::ostream& os, BinOp<Mul, L, R> const& v)
{
    return os << "(" << v.l_ << "*" << v.r_ << ")";
}
#endif
}  // namespace kyadet

double f2(int a, int b)
{
    return std::exp(a + a + b);
}

int main()
{
    using namespace kyadet;

    constexpr auto mgr = ExprMgr<2>{};

    constexpr auto x = mgr.v_arr[0];
    constexpr auto y = mgr.v_arr[1];
    constexpr auto z = exp(x + x + y);
    constexpr auto dzdx = z.diff(x);
    //std::cout << dzdx << std::endl;

    constexpr auto f = mgr.create_func(z);
    //constexpr auto dfdx = mgr.create_func(dzdx);

    /*
    double d = 0;
    for (int i = 0; i < 10000000; ++i) {
        f(1, 2);
        f2(1, 2);
    }
    */

    std::cout << f(1, 2) << std::endl;

    /*
    std::cout << f(1, 2) << std::endl;
    std::cout << dfdx(1, 2) << std::endl;
     */

    return 0;
}
