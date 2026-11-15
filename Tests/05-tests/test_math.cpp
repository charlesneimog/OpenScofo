#include <gtest/gtest.h>
#include <OpenScofo.hpp>

namespace OpenScofo {

TEST(ParseFloatingPoint, AcceptsCompleteDecimalValues) {
    double doubleValue = 0.0;
    float floatValue = 0.0f;

    EXPECT_TRUE(ParseDouble("-12.5e2", doubleValue));
    EXPECT_DOUBLE_EQ(doubleValue, -1250.0);
    EXPECT_TRUE(ParseFloat("0.125", floatValue));
    EXPECT_FLOAT_EQ(floatValue, 0.125f);
}

TEST(ParseFloatingPoint, RejectsTextOutsideTheNumber) {
    double value = 42.0;

    EXPECT_FALSE(ParseDouble("", value));
    EXPECT_FALSE(ParseDouble(" 1.5", value));
    EXPECT_FALSE(ParseDouble("1.5 trailing", value));
}

} // namespace OpenScofo
