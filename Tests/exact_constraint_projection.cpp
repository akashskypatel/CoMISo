#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <CoMISo/Utils/ExactConstraintProjection.hh>

TEST(ExactConstraint, basic) {
    COMISO::ExactConstraintProjection proj;
    const auto A = (Eigen::MatrixXd(2, 3) <<
        1, 1, 0,
        0, 1, -1).finished();
    const auto b = (Eigen::VectorXd(2) << 3, 0).finished();
    const auto x = (Eigen::VectorXd(3) << 1.3, 1.5, 1.6).finished();
    proj.initialize(A, b);
    Eigen::VectorXd xp = x;
    bool success = proj.project(xp);
    ASSERT_EQ(success, true);
    EXPECT_EQ(xp[1], xp[2]);
    EXPECT_NE(A*x, b);
    EXPECT_EQ(A*xp, b);
    EXPECT_DOUBLE_EQ(xp[0], 1.4);
    EXPECT_DOUBLE_EQ(xp[1], 1.6);
    EXPECT_DOUBLE_EQ(xp[1], 1.6);
}
