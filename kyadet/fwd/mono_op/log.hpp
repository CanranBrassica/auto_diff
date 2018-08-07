#pragma once

#include "../var.hpp"
#include "../op_list.hpp"

namespace kyadet::fwd
{
struct Log
{
    template <class X>
    static constexpr auto func(X const& x)
    {
        return [x_func = x.func()](auto const& args) noexcept(noexcept(log(x.func()(args))))
        {
            using std::log;
            return log(x_func(args));
        };
    }

    template <class X, size_t N>
    static constexpr auto diff(X const& x, Var<N> const& v)
    {   // d.log(x)/dv =  dx/dv / x
        return x.diff(v) / x;
    }
};  // class Log
}  // namespace kyadet::fwd
