#pragma once

#include <cstddef>

#include "var.hpp"

namespace kyadet::fwd {

struct Zero {
    constexpr auto func() const noexcept {
        return Zero{};  // λx.Op(x)   <- 関数オブジェクト
    }

    template<size_t M>
    constexpr auto diff(Var <M> const &v) const noexcept {
        return Zero{};  // d.zero/d.v  <- 式木
    }

    template<class T>
    constexpr operator T() const noexcept { return T{0}; }
};

struct One {
    constexpr auto func() const noexcept {
        return One{};  // λx.Op(x)   <- 関数オブジェクト
    }

    template<size_t M>
    constexpr auto diff(Var <M> const &v) const noexcept {
        return Zero{};  // d.zero/d.v  <- 式木
    }

    template<class T>
    constexpr operator T() const noexcept { return T{1}; }
};

template<class T>
struct Constant {
    T v_;

    constexpr Constant(T v) : v_(v) {}

    constexpr auto func() const noexcept {
        return Constant{v_};  // λx.Op(x)   <- 関数オブジェクト
    }

    template<size_t M>
    constexpr auto diff(Var <M> const &v) const noexcept {
        return Zero{};  // d.zero/d.v  <- 式木
    }

    constexpr operator T() const noexcept { return v_; }    // cast operator
};


}  // namespace kyadet::fwd
