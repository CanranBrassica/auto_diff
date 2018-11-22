#pragma once

#include <tuple>
#include <type_traits>
#include "var.hpp"

namespace kyadet::fwd
{

template <size_t Size>
constexpr auto create_v()
{
    return std::tuple_cat(create_v<Size - 1>(), std::tuple<Var<Size - 1>>{Var<Size - 1>{}});
}

template <>
constexpr auto create_v<1>()
{
    return std::tuple<Var<0>>{Var<0>{}};
}


template <size_t N>
struct ExprMgr
{
    decltype(create_v<N>()) v = create_v<N>();

    constexpr ExprMgr() = default;

    template <class Expr>
    constexpr auto create_func(Expr&& expr) const
    {
        return create_func_impl(std::forward<Expr>(expr), std::is_scalar<Expr>{});
    }

private:
    template <class T>
    constexpr static auto create_func_impl(T&& expr, std::true_type)
    {
        return [v = expr](auto... args) noexcept { return v; };
    }

    template <class T>
    constexpr static auto create_func_impl(T&& expr, std::false_type)
    {
        return [f = expr.func()](auto... args)
            // g++ 8.0.0以上ではコンパイルエラーになるのでnoexcept諦める
            noexcept(noexcept(
                expr.func()(
                    std::declval<
                        std::array<
                            std::remove_reference_t<
                                decltype(
                                    std::get<0>(
                                        std::declval<
                                            std::tuple<
                                                decltype(args)...>>()))>,
                            N>>())))
        {
            static_assert(sizeof...(args) == N, "A number of args must be N.");
            return f(
                std::array<
                    std::remove_reference_t<
                        decltype(std::get<0>(
                            std::declval<std::tuple<decltype(args)...>>()))>,
                    N>{args...});
        };
    }
};  // class ExprMgr
}  // namespace kyadet::fwd
