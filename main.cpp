#include <iostream>
#include <cmath>
#include <type_traits>
#include <tuple>

#define USE_STRING 1

namespace kyadet
{
class Zero;

class One;

// Varクラス
template <size_t N>
struct Var
{
    template <size_t M>
    constexpr auto diff(Var<M> const& v) const
    {
        return Zero{};
    }

    constexpr auto diff(Var<N> const& v) const
    {
        return One{};
    }

    constexpr auto func() const
    {
        return [](auto const& args) {  // noexcept(operator[]) == true
            return args[N];            // argsから該当する部分だけ返す
        };
    }

private:
    constexpr Var() = default;

};  // class Var

template <size_t Size>
constexpr auto create_v()
{
    return std::tuple_cat(create_v<Size - 1>(), std::tuple<Var<Size - 1>>{Var<Size - 1>{}});
}

template <>
constexpr auto create_v<1>()
{
    return std::tuple<Var<0>>{Var<0>{}};
}


template <size_t N>
struct ExprMgr
{
    decltype(create_v<N>()) v = create_v<N>();

    constexpr ExprMgr() = default;

    template <class Expr>
    constexpr auto create_func(Expr&& expr) const
    {
        return create_func_impl(std::forward<Expr>(expr), std::is_scalar<Expr>{});
    }

private:
    template <class T>
    constexpr static auto create_func_impl(T&& expr, std::true_type)
    {
        return [v = expr](auto... args) noexcept { return v; };
    }

    template <class T>
    constexpr static auto create_func_impl(T&& expr, std::false_type)
    {
        return [f = expr.func()](auto... args)
            // g++ 8.0.0以上ではコンパイルエラーになるのでnoexcept諦める
            noexcept(noexcept(
                expr.func()(
                    std::declval<
                        std::array<
                            std::remove_reference_t<
                                decltype(
                                    std::get<0>(
                                        std::declval<
                                            std::tuple<
                                                decltype(args)...>>()))>,
                            N>>())))
        {
            static_assert(sizeof...(args) == N, "A number of args must be N.");
            return f(
                std::array<
                    std::remove_reference_t<
                        decltype(std::get<0>(
                            std::declval<std::tuple<decltype(args)...>>()))>,
                    N>{args...});
        };
    }
};  // class ExprMgr


struct Zero
{
    constexpr auto func() const noexcept
    {
        return 0;  // λx.Op(x)   <- 関数オブジェクト
    }

    template <size_t N>
    constexpr auto diff(Var<N> const& v) const noexcept
    {
        return Zero{};  // d.zero/d.v  <- 式木
    }
};

struct One
{
    constexpr auto func() const noexcept
    {
        return 1;  // λx.Op(x)   <- 関数オブジェクト
    }

    template <size_t N>
    constexpr auto diff(Var<N> const& v) const noexcept
    {
        return Zero{};  // d.zero/d.v  <- 式木
    }
};

template <class Op, class X>
struct MonoOp
{
    X x_;
    constexpr MonoOp(X const& x) : x_(x) {}

    constexpr auto func() const noexcept
    {
        return Op::func(x_);  // λx.Op(x)   <- 関数オブジェクト
    }

    template <size_t N>
    constexpr auto diff(Var<N> const& v) const noexcept
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

    constexpr auto func() const noexcept
    {
        return Op::func(l_, r_);  // λx.λy.Op(x, y)   <- 関数オブジェクト
    }

    template <size_t N>
    constexpr auto diff(Var<N> const& v) const noexcept
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

class Log;

class Add;

class Sub;

class Mul;

class Div;

template <class X>
constexpr auto exp(X const& x) noexcept
{
    return MonoOp<Exp, X>{x};
}

template <class X>
constexpr auto exp(MonoOp<Log, X> const& logx) noexcept
{
    return logx.x_;
}

template <class X>
constexpr auto log(X const& x) noexcept
{
    return MonoOp<Log, X>{x};
}

template <class X>
constexpr auto log(MonoOp<Exp, X> const& expx) noexcept
{
    return expx.x_;
}

template <class L, class R>
constexpr auto operator+(L const& l, R const& r) noexcept
{
    return BinOp<Add, L, R>{l, r};
}

template <class L, class R>
constexpr auto operator-(L const& l, R const& r) noexcept
{
    return BinOp<Sub, L, R>{l, r};
}

template <class L, class R>
constexpr auto operator*(L const& l, R const& r) noexcept
{
    return BinOp<Mul, L, R>{l, r};
}

template <class R>
constexpr auto operator*(One const&, R const& r) noexcept
{
    return r;
}

template <class L>
constexpr auto operator*(L const& l, One const&)noexcept
{
    return l;
}

template <class R>
constexpr auto operator*(Zero const&, R const&)noexcept
{
    return Zero{};
}

template <class L>
constexpr auto operator*(L const&, Zero const&)noexcept
{
    return Zero{};
}

template <class L, class R>
constexpr auto operator/(L const& l, R const& r) noexcept
{
    return BinOp<Div, L, R>{l, r};
}

// Mono Operator Definition
struct Exp
{
    template <class X>
    static constexpr auto func(X const& x)
    {
        using std::exp;
        return [x_func = x.func()](auto const& args) noexcept(noexcept(exp(x.func()(args))))
        {
            using std::exp;
            return exp(x_func(args));
        };
    }

    template <class X, size_t N>
    static constexpr auto diff(X const& x, Var<N> const& v)
    {  // d.exp(x)/dv = exp(x) * dx/dv
        return exp(x) * x.diff(v);
    }
};  // class Exp

struct Log
{
    template <class X>
    static constexpr auto func(X const& x)
    {
        using std::log;
        return [x_func = x.func()](auto const& args) noexcept(noexcept(log(x.func()(args))))
        {
            using std::log;
            return log(x_func(args));
        };
    }

    template <class X, size_t N>
    static constexpr auto diff(X const& x, Var<N> const& v)
    {  // d.log(x)/dv =  dx/dv / x
        return x.diff(v) / x;
    }
};  // class Exp

// Binary Operator Definition
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
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::false_type, std::false_type)
    {
        return l.diff(v) + r.diff(v);  // d(l+r)/dv = dl/dv + dr/dv
    }

    // Scalar, Expr
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::true_type, std::false_type)
    {
        return r.diff(v);
    }

    // Expr, Scalar
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::false_type, std::true_type)
    {
        return l.diff(v);
    }

public:
    template <class L, class R, size_t N>
    static constexpr auto diff(L const& l, R const& r, Var<N> const& v)
    {
        return diff_impl(l, r, v, std::is_scalar<L>{}, std::is_scalar<R>{});
    }
};  // class Add

struct Sub
{
    template <class L, class R>
    static constexpr auto func(L const& l, R const& r)
    {
        return func_impl(l, r, std::is_scalar<L>{}, std::is_scalar<R>{});
    }

private:
    // Expr, Expr
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::false_type, std::false_type)
    {
        return [ l_func = l.func(), r_func = r.func() ](auto const& args) noexcept(noexcept(l.func()(args)) && noexcept(r.func()(args)))
        {
            return l_func(args) - r_func(args);
        };
    }

    // Scalar, Expr
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::true_type, std::false_type)
    {
        return [ l, r_func = r.func() ](auto const& args) noexcept(noexcept(r.func()(args)))
        {
            return l - r_func(args);
        };
    }

    // Expr, Scalar
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::false_type, std::true_type)
    {
        return [ l_func = l.func(), r ](auto const& args) noexcept(noexcept(l.func()(args)))
        {
            return l_func(args) - r;
        };
    }

public:
    template <class L, class R, size_t N>
    static constexpr auto diff(L const& l, R const& r, Var<N> const& v)
    {
        return diff_impl(l, r, v, std::is_scalar<L>{}, std::is_scalar<R>{});
    }

private:
    // Expr, Expr
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::false_type, std::false_type)
    {
        return l.diff(v) - r.diff(v);  // d(l-r)/dv = dl/dv - dr/dv
    }

    // Scalar, Expr
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::true_type, std::false_type)
    {
        return -r.diff(v);
    }

    // Expr, Scalar
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::false_type, std::true_type)
    {
        return l.diff(v);
    }
};  // class Sub

struct Mul
{
    template <class L, class R>
    static constexpr auto func(L const& l, R const& r)
    {
        return func_impl(l, r, std::is_scalar<L>{}, std::is_scalar<R>{});
    }

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
    template <class L, class R, size_t N>
    static constexpr auto diff(L const& l, R const& r, Var<N> const& v)
    {
        return diff_impl(l, r, v, std::is_scalar<L>{}, std::is_scalar<R>{});
    }

private:
    // Expr, Expr
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::false_type, std::false_type)
    {
        return l.diff(v) * r + l * r.diff(v);  // d(l*r)/dv = dl/dv * r + l * dr/dv
    }

    // Scalar, Expr
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::true_type, std::false_type)
    {
        return l * r.diff(v);
    }

    // Expr, Scalar
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::false_type, std::true_type)
    {
        return l.diff(v) * r;
    }
};  // class Mul

struct Div
{
    template <class L, class R>
    static constexpr auto func(L const& l, R const& r)
    {
        return func_impl(l, r, std::is_scalar<L>{}, std::is_scalar<R>{});
    }

private:
    // Expr, Expr
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::false_type, std::false_type)
    {
        return [ l_func = l.func(), r_func = r.func() ](auto const& args) noexcept(noexcept(l.func()(args)) && noexcept(r.func()(args)))
        {
            return l_func(args) / r_func(args);
        };
    }

    // Scalar, Expr
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::true_type, std::false_type)
    {
        return [ l, r_func = r.func() ](auto const& args) noexcept(noexcept(r.func()(args)))
        {
            return l / r_func(args);
        };
    }

    // Expr, Scalar
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::false_type, std::true_type)
    {
        return [ l_func = l.func(), r ](auto const& args) noexcept(noexcept(l.func()(args)))
        {
            return l_func(args) / r;
        };
    }

public:
    template <class L, class R, size_t N>
    static constexpr auto diff(L const& l, R const& r, Var<N> const& v)
    {
        return diff_impl(l, r, v, std::is_scalar<L>{}, std::is_scalar<R>{});
    }

private:
    // Expr, Expr
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::false_type, std::false_type)
    {
        return l.diff(v) / r - l * r.diff(v) / (r * r);  // d(l/r)/dv = l'/r - lr'/r^2
    }

    // Scalar, Expr
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::true_type, std::false_type)
    {
        return l * r.diff(v) / (r * r);
    }

    // Expr, Scalar
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::false_type, std::true_type)
    {
        return l.diff(v) / r;
    }
};  // class Div

#if USE_STRING
template <size_t N>
std::ostream& operator<<(std::ostream& os, Var<N> const& v)
{
    return os << "v[" << N << "]";
}

std::ostream& operator<<(std::ostream& os, Zero const&)
{
    return os << 0;
}

std::ostream& operator<<(std::ostream& os, One const&)
{
    return os << 1;
}

template <class X>
std::ostream& set_ostream_impl(std::ostream& os, MonoOp<Exp, X> const& v)
{
    return os << "exp(" << v.x_ << ")";
}

template <class L, class R>
std::ostream& set_ostream_impl(std::ostream& os, BinOp<Add, L, R> const& v)
{
    return os << "(" << v.l_ << "+" << v.r_ << ")";
}

template <class L, class R>
std::ostream& set_ostream_impl(std::ostream& os, BinOp<Sub, L, R> const& v)
{
    return os << "(" << v.l_ << "-" << v.r_ << ")";
}

template <class L, class R>
std::ostream& set_ostream_impl(std::ostream& os, BinOp<Mul, L, R> const& v)
{
    return os << "(" << v.l_ << "*" << v.r_ << ")";
}

template <class L, class R>
std::ostream& set_ostream_impl(std::ostream& os, BinOp<Div, L, R> const& v)
{
    return os << "(" << v.l_ << "/" << v.r_ << ")";
}
#endif
}  // namespace kyadet

int main()
{
    using namespace kyadet;

    constexpr auto mgr = ExprMgr<2>{};

    constexpr auto x = std::get<0>(mgr.v);
    constexpr auto y = std::get<1>(mgr.v);
    constexpr auto z = log(exp(x + x / y));
    constexpr auto dzdx = z.diff(x);
    std::cout << z << std::endl;
    std::cout << dzdx << std::endl;

    constexpr auto f = mgr.create_func(z);
    //constexpr auto dfdx = mgr.create_func(dzdx);

    //std::cout << f(1, 2) << std::endl;

    return 0;
}
