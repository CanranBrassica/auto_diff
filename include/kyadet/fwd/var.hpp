#pragma once

#include <cstddef>

namespace kyadet::fwd
{
class Zero;
class One;

template <size_t N>
struct Var
{
    template <size_t M>
    constexpr auto diff([[maybe_unused]] Var<M> const& v) const
    {
        return Zero{};
    }

    constexpr auto diff([[maybe_unused]] Var<N> const& v) const
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

}  // namespace kyadet::fwd
