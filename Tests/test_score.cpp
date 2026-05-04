// #include <gtest/gtest.h>
// #include <OpenScofo.hpp>
//
// std::filesystem::path assets = TEST_DATA_DIR;
//
// // Test fixture for OpenScofo
// class Score : public ::testing::Test {
//   protected:
//     OpenScofo::Score *score;
//
//     void SetUp() override {
//         score = new OpenScofo::Score();
//     }
//
//     void TearDown() override {
//         delete score;
//     }
// };
//
// // 1. Test basic initialization
// TEST_F(Score, Initialization) {
//     ASSERT_NE(score, nullptr);
// }
//
// // WRong trill
// TEST_F(Score, WrongTRILL) {
//     (score->Parse(assets / "wrongtrill.txt"));
// }
//
// // 1. Score Parse
// TEST_F(Score, ScoreParse) {
//     // OpenScofo::States states = score->Parse(assets / "bwv-1013.txt");
// }
//
// TEST_F(Score, WrongFFTSize) {
//     EXPECT_THROW(score->Parse(assets / "wrongfftsize.txt"), std::runtime_error);
// }
//
// TEST_F(Score, WrongHopSize) {
//     EXPECT_THROW(score->Parse(assets / "wronghopsize.txt"), std::runtime_error);
// }
