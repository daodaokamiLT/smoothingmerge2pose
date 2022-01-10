// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// @author: Yi Liu
// @date:   2021-01-15
//

#pragma once
#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>

#include "numeric_types.h"

// namespace means smoothing vio vps fusion
namespace svv_fusion {
/*
 * For Eigen
 */

template <typename T, int rows, int cols, int options = Eigen::RowMajor>
using Matrix = Eigen::Matrix<T, rows, cols, options>;

template <typename T>
using Matrix2 = Matrix<T, 2, 2>;

template <typename T>
using Matrix3 = Matrix<T, 3, 3>;

template <typename T>
using Matrix4 = Matrix<T, 4, 4>;

template <typename T>
using Matrix3x2 = Matrix<T, 3, 2>;

template <typename T>
using Matrix3x4 = Matrix<T, 3, 4>;

template <typename T, int rows>
using Vector = Matrix<T, rows, 1, Eigen::ColMajor>;

template <typename T>
using Vector2 = Vector<T, 2>;

template <typename T>
using Vector3 = Vector<T, 3>;

template <typename T>
using Vector4 = Vector<T, 4>;

template <typename T>
using Vector6 = Vector<T, 6>;

/*==================== Quaternion types. ====================*/
template <typename T>
using Quaternion = Eigen::Quaternion<T>;

using Quaternionf = Quaternion<float>;
using Quaterniond = Quaternion<double>;

/*==================== AngleAxis types ====================*/
template <typename T>
using AngleAxis = Eigen::AngleAxis<T>;

using AngleAxisf = AngleAxis<float>;
using AngleAxisd = AngleAxis<double>;

/*==================== se3 pose types. ====================*/
template <typename T>
using PoseSE3 = Matrix3x4<T>;

using PoseSE3f = PoseSE3<float>;
using PoseSE3d = PoseSE3<double>;

using Mat2f = Matrix2<float>;
using Mat2d = Matrix2<double>;
using Mat3f = Matrix3<float>;
using Mat3d = Matrix3<double>;
using Mat4f = Matrix4<float>;
using Mat4d = Matrix4<double>;
using Mat34d = PoseSE3d;

using Mat32d = Matrix<double, 3, 2>;
using Mat36d = Matrix<double, 3, 6>;
using Mat6d = Matrix<double, 6, 6>;
using Mat63d = Matrix<double, 6, 3>;
using Mat69d = Matrix<double, 6, 9>;
using Mat9d = Matrix<double, 9, 9>;
using Mat96d = Matrix<double, 9, 6>;
using Mat126d = Matrix<double, 12, 6>;
using Mat612d = Matrix<double, 6, 12>;
using Mat312d = Matrix<double, 3, 12>;

using Mat12d = Matrix<double, 12, 12>;
using Mat15d = Matrix<double, 15, 15>;
using Mat152d = Matrix<double, 15, 2>;
using Mat153d = Matrix<double, 15, 3>;
using Mat156d = Matrix<double, 15, 6>;
using Mat1512d = Matrix<double, 15, 12>;

using MatXf = Matrix<float, Eigen::Dynamic, Eigen::Dynamic>;
using MatXd = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

using Vec2f = Vector2<float>;
using Vec2d = Vector2<double>;
using Vec3f = Vector3<float>;
using Vec3d = Vector3<double>;
using Vec4f = Vector4<float>;
using Vec4d = Vector4<double>;
using Vec6d = Vector6<double>;
using Vec9d = Eigen::Matrix<double, 9, 1, Eigen::ColMajor>;
using Vec15d = Eigen::Matrix<double, 15, 1, Eigen::ColMajor>;

/*==================== redefinetion types ====================*/
using Vec2_t = Vector2<NumericType>;
using Vec3_t = Vector3<NumericType>;
using Vec6_t = Vector6<NumericType>;

using Mat2_t = Matrix2<NumericType>;
using Mat3_t = Matrix3<NumericType>;
using Mat4_t = Matrix4<NumericType>;
using Mat36_t = Matrix<NumericType, 3, 6>;
using Mat6_t = Matrix<NumericType, 6, 6>;
using Mat312_t = Matrix<NumericType, 3, 12>;
using Mat126_t = Matrix<NumericType, 12, 6>;
using Mat612_t = Matrix<NumericType, 6, 12>;
using MatX_t = Matrix<NumericType, Eigen::Dynamic, Eigen::Dynamic>;
using PoseSE3_t = PoseSE3<NumericType>;
using Quaternion_t = Quaternion<NumericType>;

struct Posed_t{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Quaterniond q_wc;
  Vec3d t_wc;
  double timestamp;
  double q_var;
  double t_var;

  Posed_t(){
    q_wc = Quaterniond::Identity();
    t_wc = Vec3d::Constant(std::nan(""));
    q_var = std::numeric_limits<double>::max();
    t_var = std::numeric_limits<double>::max();
  }

  Posed_t(Quaterniond& q, Vec3d t){
    q_wc = q;
    t_wc = t;
    q_var = std::numeric_limits<double>::max();
    t_var = std::numeric_limits<double>::max();
  }

  Posed_t(Quaterniond& q, Vec3d t, double qv, double tv){
    q_wc = q;
    t_wc = t;
    q_var = qv;
    t_var = tv;
  }
};

void inline deltaPosed(const Posed_t& p1, const Posed_t& p0, Posed_t& p10){
  p10.q_wc = p1.q_wc.toRotationMatrix().transpose() * p0.q_wc.toRotationMatrix();
  p10.t_wc = p1.q_wc.toRotationMatrix().transpose() * (p0.t_wc - p1.t_wc);
  p10.q_var = p1.q_var;
  p10.t_var = p1.t_var;
}


struct Posef_t{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Quaternionf q_wc;
  Vec3f t_wc;
  double timestamp;
  double q_var;
  double t_var;
};

// px, py, pz, qw, qx, qy, qz
template <typename Scalar,
          typename = enable_if_t<std::is_floating_point<Scalar>::value>>
inline Matrix3x4<Scalar> TransArrayToPose(const std::array<Scalar, 7>& pose) {
  Quaternion<Scalar> quat;
  quat.w() = pose[3];
  quat.x() = pose[4];
  quat.y() = pose[5];
  quat.z() = pose[6];

  Matrix3x4<Scalar> pose_eig = Matrix3x4<Scalar>::Identity();
  pose_eig.template block<3, 1>(0, 3)[0] = pose[0];
  pose_eig.template block<3, 1>(0, 3)[1] = pose[1];
  pose_eig.template block<3, 1>(0, 3)[2] = pose[2];
  pose_eig.template block<3, 3>(0, 0) = quat.toRotationMatrix();
  return pose_eig;
}

/*==================== data invalid values ====================*/
const Vec3f kInvalid3dVecf = Vec3f::Constant(std::nan(""));
const Vec3d kInvalid3dVecd = Vec3d::Constant(std::nan(""));
const Mat3f kInvalidMat3f = Mat3f::Constant(std::nan(""));
const Mat3d kInvalidMat3d = Mat3d::Constant(std::nan(""));
const PoseSE3f kInvalidPose3dVecf = PoseSE3f::Constant(std::nan(""));
const PoseSE3d kInvalidPose3dVecd = PoseSE3d::Constant(std::nan(""));
const Vec3d kInvalidVec3_t = Vec3d::Constant(std::nan(""));
const Mat3d kInvalidMat3_t = Mat3d::Constant(std::nan(""));
const PoseSE3d kInvalidPoseSE3_t = PoseSE3d::Constant(std::nan(""));

/*==================== eigen alignment ====================*/

/*
 *  traits
 */

template <class Derived>
struct IsEigenDense
    : public std::is_base_of<Eigen::DenseBase<Derived>, Derived> {};

template <typename Mat>
struct Check<Mat, enable_if_t<IsEigenDense<Mat>::value>> {
  static_assert(std::is_floating_point<typename Mat::Scalar>::value,
                "value must be a float!");
  static bool IsValid(const Mat &value) { return !value.hasNaN(); }
};

/// If the Vector type is of fixed size, then IsFixedSizeVector::value will be
/// true.
template <typename Vector, int NumDimensions,
          typename = typename std::enable_if<
              Vector::RowsAtCompileTime == NumDimensions &&
              Vector::ColsAtCompileTime == 1>::type>
struct IsFixedSizeVector : std::true_type {};

template <typename Matrix, int RowDimensions, int ColDimensions,
          typename = typename std::enable_if<
              Matrix::RowsAtCompileTime == RowDimensions &&
              Matrix::ColsAtCompileTime == ColDimensions>::type>
struct IsFixedSizeMatrix : std::true_type {};

template <typename T>
using AlignedVector = std::vector<T, Eigen::aligned_allocator<T>>;

template <typename T>
using AlignedDeque = std::deque<T, Eigen::aligned_allocator<T>>;

template <typename K, typename V>
using AlignedMap = std::map<K, V, std::less<K>,
                            Eigen::aligned_allocator<std::pair<K const, V>>>;

template <typename K, typename V>
using AlignedUnorderedMap =
    std::unordered_map<K, V, std::hash<K>, std::equal_to<K>,
                       Eigen::aligned_allocator<std::pair<K const, V>>>;

/*==================== vector eigen ====================*/
using VVec2f = AlignedVector<Vec2f>;
using VVec2d = AlignedVector<Vec2d>;
using VVec3f = AlignedVector<Vec3f>;
using VVec3d = AlignedVector<Vec3d>;
using VVec6d = AlignedVector<Vec6d>;
using VMat2f = AlignedVector<Mat2f>;
using VMat2d = AlignedVector<Mat2d>;
using VMat3d = AlignedVector<Mat3d>;
using VMat6d = AlignedVector<Mat6d>;
using VMat34d = AlignedVector<Mat34d>;
using VMat63d = AlignedVector<Mat63d>;
using VVec3_t = AlignedVector<Vec3_t>;
using VMat2_t = AlignedVector<Mat2_t>;

/*==================== map eigen ====================*/
using VMap3f = AlignedMap<PointID, Vec3f>;
using VMap3d = AlignedMap<PointID, Vec3d>;

template <class MATRIX>
bool EqualWithAbsTol(const Eigen::DenseBase<MATRIX>& A,
                     const Eigen::DenseBase<MATRIX>& B, double tol = 1e-9) {
  const size_t n1 = A.cols(), m1 = A.rows();
  const size_t n2 = B.cols(), m2 = B.rows();

  if (m1 != m2 || n1 != n2) return false;

  for (size_t i = 0; i < m1; i++)
    for (size_t j = 0; j < n1; j++) {
      if (!fpEqual(A(i, j), B(i, j), tol, false)) {
        return false;
      }
    }
  return true;
}

template <class V>
bool AssertEqual(const V& expected, const V& actual, double tol = 1e-9) {
  if (EqualWithAbsTol(expected, actual, tol)) return true;
  std::cout << "Not equal:\n";
  std::cout << "Expected:\n" << expected << "\n";
  std::cout << "Actual:\n" << actual << "\n";
  return false;
}
}  // namespace crispy