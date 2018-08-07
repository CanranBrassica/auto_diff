#pragma once

#include <type_traits>
#include "../var.hpp"
#include "../op_list.hpp"

namespace kyadet::fwd
{

struct Pow
{
    template <class L, class R>
    static constexpr auto func(L const& l, R const& r)
    {
        return func_impl(l, r, std::is_scalar<L>{}, std::is_scalar<R>{});
    }

private:
    // Expr, Expr   l^r = log(exp(l^r)) = log(r * exp(l))
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::false_type, std::false_type)
    {
        return [ l_func = l.func(), r_func = r.func() ](auto const& args) noexcept(noexcept(l.func()(args)) && noexcept(r.func()(args)))
        {
            using std::exp;
            using std::log;
            return log(r_func(args) * exp(l_func(args)));
        };
    }

    // Scalar, Expr
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::true_type, std::false_type)
    {
        return [ l, r_func = r.func() ](auto const& args) noexcept(noexcept(r.func()(args)))
        {
            using std::pow;
            return pow(l ,r_func(args));
        };
    }

    // Expr, Scalar
    template <class L, class R>
    static constexpr auto func_impl(L const& l, R const& r, std::false_type, std::true_type)
    {
        return [ l_func = l.func(), r ](auto const& args) noexcept(noexcept(l.func()(args)))
        {
            using std::pow;
            return pow(l_func(args) , r);
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
        return log(r * exp(l)).diff(v);   // l^r = log(exp(l^r)) = log(r * exp(l))
    }

    // Scalar, Expr
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::true_type, std::false_type)
    {
        return pow(l, r) * log(l) * r.diff(v);
    }

    // Expr, Scalar
    template <class L, class R, size_t N>
    static constexpr auto diff_impl(L const& l, R const& r, Var<N> const& v, std::false_type, std::true_type)
    {
        return r * pow(l, r - One{}) * l.diff(v);
    }

};  // class Pow
}  // namespace kyadet::fwd
