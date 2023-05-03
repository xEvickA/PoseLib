#ifndef POSELIB_RELPOSE_6PT_PLANAR_H_
#define POSELIB_RELPOSE_6PT_PLANAR_H_

#include "PoseLib/camera_pose.h"

#include <Eigen/Dense>
#include <vector>

namespace poselib {

// Computes the planar essential matrix from 6 point correspondences. Returning up to 9 solutions.
int relpose_6pt_planar(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2,
                       std::vector<Eigen::Matrix3d> *output);
    
int relpose_6pt_planar(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2,
                       Eigen::Matrix<std::complex<double>, 4, 30> &sols);

}; // namespace poselib

#endif