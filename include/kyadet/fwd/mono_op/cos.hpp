#pragma once

#include "../var.hpp"
#include "../op_list.hpp"

namespace kyadet::fwd
{
struct Cos
{
    template <class X>
    static constexpr auto func(X const& x)
    {
        return [x_func = x.func()](auto const& args) noexcept(noexcept(sin(x.func()(args))))
        {
            using std::cos;
            return cos(x_func(args));
        };
    }

    template <class X, size_t N>
    static constexpr auto diff(X const& x, Var<N> const& v)
    {
        // d.cos(x)/dv = -sin(x) * dx/dv
        using std::sin;
        return -sin(x) * x.diff(v);
    }
};  // class Sin
}  // namespace kyadet::fwd
