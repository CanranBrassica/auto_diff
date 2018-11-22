#pragma once

#include <cmath>
#include <exception>

#include "constant.hpp"
#include "mono_op.hpp"
#include "bin_op.hpp"

namespace kyadet::fwd
{

// mono op
class Exp;

class Log;

class Sin;

class Cos;

// bin op
class Add;

class Sub;

class Mul;

class Div;

class Pow;

// exp
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

// log
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

template <>
constexpr auto log(One const&) noexcept
{
    return Zero{};
}

constexpr auto log(Zero const&)
{
    //throw std::exception{"division by zero"};
}

template <class X>
constexpr auto sin(X const& x) noexcept
{
    return MonoOp<Sin, X>{x};
}

template <class X>
constexpr auto cos(X const& x) noexcept
{
    return MonoOp<Cos, X>{x};
}

// operator+
template <class L, class R>
constexpr auto operator+(L const& l, R const& r) noexcept
{
    return BinOp<Add, L, R>{l, r};
}

template <class L>
constexpr auto operator+(L const& l, Zero const&) noexcept
{
    return l;
}

template <class R>
constexpr auto operator+(Zero const&, R const& r) noexcept
{
    return r;
}

// operator-
template <class L, class R>
constexpr auto operator-(L const& l, R const& r) noexcept
{
    return BinOp<Sub, L, R>{l, r};
}

template <class L>
constexpr auto operator-(L const& l, Zero const&) noexcept
{
    return l;
}

template <class R>
constexpr auto operator-(Zero const&, R const& r) noexcept
{
    return r;
}

// operator*
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

// operator/
template <class L, class R>
constexpr auto operator/(L const& l, R const& r) noexcept
{
    return BinOp<Div, L, R>{l, r};
}

template <class R>
constexpr auto operator/(Zero const&, R const&) noexcept
{
    return Zero{};
}

template <class L>
constexpr auto operator/(L const&, Zero const&)
{
    throw std::exception{"division by zero"};
}

// pow
template <class L, class R>
constexpr auto pow(L const& l, R const& r) noexcept
{
    return BinOp<Pow, L, R>{l, r};
}

template <class R>
constexpr auto pow(Zero const&, R const&) noexcept
{
    return Zero{};
}

template <class R>
constexpr auto pow(One const&, R const&) noexcept
{
    return One{};
}

template <class L>
constexpr auto pow(L const&, Zero const&) noexcept
{
    return One{};
}

template <class L>
constexpr auto pow(L const& l, One const&) noexcept
{
    return l;
}

}  // namespace kyadet::fwd

#include "mono_op/exp.hpp"
#include "mono_op/log.hpp"
#include "mono_op/sin.hpp"
#include "mono_op/cos.hpp"

#include "bin_op/add.hpp"
#include "bin_op/sub.hpp"
#include "bin_op/mul.hpp"
#include "bin_op/div.hpp"
#include "bin_op/pow.hpp"
