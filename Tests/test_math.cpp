#include <gtest/gtest.h>
#include <OpenScofo.hpp>

namespace OpenScofo {

// Kappa
TEST(OnlineForwardKappaTest, A2Properties) {
    OnlineForward f(44100, 1024, 512);

    double a1 = f.A2(1.0);
    double a2 = f.A2(2.0);

    EXPECT_GT(a2, a1);
}

} // namespace OpenScofo
