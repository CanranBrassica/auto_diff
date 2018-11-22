#define USE_STRING 1

#include "kyadet/fwd.hpp"

#include <iostream>

int main()
{
    using namespace kyadet::fwd;

    constexpr auto mgr = ExprMgr<2>{};

    constexpr auto x = std::get<0>(mgr.v);
    constexpr auto y = std::get<1>(mgr.v);
    constexpr auto z = sin(log(pow(x, y)));
    constexpr auto dzdx = z.diff(x);
    std::cout << z << std::endl;
    std::cout << dzdx << std::endl;

    constexpr auto f = mgr.create_func(dzdx);

    std::cout << f(1.0, 2.0) << std::endl;

    return 0;
}
