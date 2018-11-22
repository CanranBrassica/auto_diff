#pragma once

#include "var.hpp"

namespace kyadet::fwd
{
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
}  // namespace kyadet::fwd

