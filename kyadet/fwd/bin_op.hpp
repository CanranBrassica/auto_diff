#pragma once

#include "var.hpp"

namespace kyadet::fwd
{
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
}  // namespace kyadet::fwd
