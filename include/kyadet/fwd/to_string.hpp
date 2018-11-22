#pragma once

#include <ostream>

#include "var.hpp"
#include "mono_op.hpp"
#include "bin_op.hpp"
#include "op_list.hpp"

namespace kyadet::fwd
{
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

template <size_t N>
std::ostream& operator<<(std::ostream& os, [[maybe_unused]] Var<N> const& v)
{
    return os << "v_" << N;
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

template <class X>
std::ostream& set_ostream_impl(std::ostream& os, MonoOp<Log, X> const& v)
{
    return os << "log(" << v.x_ << ")";
}

template <class X>
std::ostream& set_ostream_impl(std::ostream& os, MonoOp<Sin, X> const& v)
{
    return os << "sin(" << v.x_ << ")";
}

template <class X>
std::ostream& set_ostream_impl(std::ostream& os, MonoOp<Cos, X> const& v)
{
    return os << "cos(" << v.x_ << ")";
}

template <class L, class R>
std::ostream& set_ostream_impl(std::ostream& os, BinOp<Add, L, R> const& v)
{
    return os << "(" << v.l_ << " + " << v.r_ << ")";
}

template <class L, class R>
std::ostream& set_ostream_impl(std::ostream& os, BinOp<Sub, L, R> const& v)
{
    return os << "(" << v.l_ << " - " << v.r_ << ")";
}

template <class L, class R>
std::ostream& set_ostream_impl(std::ostream& os, BinOp<Mul, L, R> const& v)
{
    return os << "(" << v.l_ << " * " << v.r_ << ")";
}

template <class L, class R>
std::ostream& set_ostream_impl(std::ostream& os, BinOp<Div, L, R> const& v)
{
    return os << "(" << v.l_ << " / " << v.r_ << ")";
}

template <class L, class R>
std::ostream& set_ostream_impl(std::ostream& os, BinOp<Pow, L, R> const& v)
{
    return os << "(" << v.l_ << "^" << v.r_ << ")";
}
}  // namespace kyadet::fwd
