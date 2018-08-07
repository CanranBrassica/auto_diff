#pragma once

#include "../var.hpp"
#include "../op_list.hpp"

namespace kyadet::fwd
{
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

    template <class X, size_t N>
    static constexpr auto diff(X const& x, Var<N> const& v)
    {  // d.exp(x)/dv = exp(x) * dx/dv
        return exp(x) * x.diff(v);
    }
};  // class Exp
}  // namespace kyadet::fwd
