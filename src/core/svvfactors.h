/*******************************************************
 * Copyright (C) 2019, Aerial Robotics Group, Hong Kong University of Science and Technology
 * 
 * This file is part of VINS.
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Qin Tong (qintonguav@gmail.com)
 *******************************************************/

#pragma once
#include <ceres/ceres.h>
#include <ceres/rotation.h>


namespace svv_fusion{
template <typename T> inline
void QuaternionInverse(const T q[4], T q_inverse[4])
{
	q_inverse[0] = q[0];
	q_inverse[1] = -q[1];
	q_inverse[2] = -q[2];
	q_inverse[3] = -q[3];
};


struct TError
{
	TError(double t_x, double t_y, double t_z, double var)
				  :t_x(t_x), t_y(t_y), t_z(t_z), var(var){}

	template <typename T>
	bool operator()(const T* tj, T* residuals) const
	{
		residuals[0] = (tj[0] - T(t_x)) / T(var);
		residuals[1] = (tj[1] - T(t_y)) / T(var);
		residuals[2] = (tj[2] - T(t_z)) / T(var);

		return true;
	}

	static ceres::CostFunction* Create(const double t_x, const double t_y, const double t_z, const double var) 
	{
	  return (new ceres::AutoDiffCostFunction<
	          TError, 3, 3>(
	          	new TError(t_x, t_y, t_z, var)));
	}

	double t_x, t_y, t_z, var;

};

struct RelativeRTError
{
	RelativeRTError(double t_x, double t_y, double t_z, 
					double q_w, double q_x, double q_y, double q_z,
					double t_var, double q_var)
				  :t_x(t_x), t_y(t_y), t_z(t_z), 
				   q_w(q_w), q_x(q_x), q_y(q_y), q_z(q_z),
				   t_var(t_var), q_var(q_var){}

	template <typename T>
	bool operator()(const T* const w_q_i, const T* ti, const T* w_q_j, const T* tj, T* residuals) const
	{
		T t_w_ij[3];
		t_w_ij[0] = tj[0] - ti[0];
		t_w_ij[1] = tj[1] - ti[1];
		t_w_ij[2] = tj[2] - ti[2];

		T i_q_w[4];
		QuaternionInverse(w_q_i, i_q_w);

		T t_i_ij[3];
		ceres::QuaternionRotatePoint(i_q_w, t_w_ij, t_i_ij);

		residuals[0] = (t_i_ij[0] - T(t_x)) / T(t_var);
		residuals[1] = (t_i_ij[1] - T(t_y)) / T(t_var);
		residuals[2] = (t_i_ij[2] - T(t_z)) / T(t_var);

		T relative_q[4];
		relative_q[0] = T(q_w);
		relative_q[1] = T(q_x);
		relative_q[2] = T(q_y);
		relative_q[3] = T(q_z);

		T q_i_j[4];
		ceres::QuaternionProduct(i_q_w, w_q_j, q_i_j);

		T relative_q_inv[4];
		QuaternionInverse(relative_q, relative_q_inv);

		T error_q[4];
		ceres::QuaternionProduct(relative_q_inv, q_i_j, error_q); 

		residuals[3] = T(2) * error_q[1] / T(q_var);
		residuals[4] = T(2) * error_q[2] / T(q_var);
		residuals[5] = T(2) * error_q[3] / T(q_var);

		return true;
	}

	static ceres::CostFunction* Create(const double t_x, const double t_y, const double t_z,
									   const double q_w, const double q_x, const double q_y, const double q_z,
									   const double t_var, const double q_var) 
	{
	  return (new ceres::AutoDiffCostFunction<
	          RelativeRTError, 6, 4, 3, 4, 3>(
	          	new RelativeRTError(t_x, t_y, t_z, q_w, q_x, q_y, q_z, t_var, q_var)));
	}

	double t_x, t_y, t_z, t_norm;
	double q_w, q_x, q_y, q_z;
	double t_var, q_var;
};

struct RelativeCRTError{
	RelativeCRTError(double t_x, double t_y, double t_z,
					double q_w, double q_x, double q_y, double q_z,
					double t_var, double q_var)
		: t_x(t_x), t_y(t_y), t_z(t_z),
		  q_w(q_w), q_x(q_x), q_y(q_y), q_z(q_z),
		  t_var(t_var), q_var(q_var) {}

	template <typename T>
	bool operator()(const T* const w_qwi, const T* w_twi, const T* q_wi, const T* t_wi, T* residuals) const
	{

		T wqwi[4];
		T wtwi[3];

		ceres::QuaternionRotatePoint(w_qwi, t_wi, wtwi);
		wtwi[0] += w_twi[0];
		wtwi[1] += w_twi[1];
		wtwi[2] += w_twi[2];

		ceres::QuaternionProduct(w_qwi, q_wi, wqwi);

		residuals[0] = (wtwi[0] - T(t_x)) / T(t_var);
		residuals[1] = (wtwi[1] - T(t_y)) / T(t_var);
		residuals[2] = (wtwi[2] - T(t_z)) / T(t_var);

		T relative_q[4];
		relative_q[0] = T(q_w);
		relative_q[1] = T(q_x);
		relative_q[2] = T(q_y);
		relative_q[3] = T(q_z);

		T relative_q_inv[4];
		QuaternionInverse(relative_q, relative_q_inv);

		T error_q[4];
		ceres::QuaternionProduct(relative_q_inv, wqwi, error_q); 

		residuals[3] = T(2) * error_q[1] / T(q_var);
		residuals[4] = T(2) * error_q[2] / T(q_var);
		residuals[5] = T(2) * error_q[3] / T(q_var);

		return true;
	}

	static ceres::CostFunction* Create(const double t_x, const double t_y, const double t_z,
									   const double q_w, const double q_x, const double q_y, const double q_z,
									   const double t_var, const double q_var)
	{
		return (new ceres::AutoDiffCostFunction<RelativeCRTError, 6, 4, 3, 4, 3>(
			new RelativeCRTError(t_x, t_y, t_z, q_w, q_x, q_y, q_z, t_var, q_var)
		));
	}

	double t_x, t_y, t_z, t_norm;
	double q_w, q_x, q_y, q_z;
	double t_var, q_var;
};


struct RelativeCRTError2{
	RelativeCRTError2(double t_x0, double t_y0, double t_z0,
					double q_w0, double q_x0, double q_y0, double q_z0,
					double t_var0, double q_var0,
					double t_x1, double t_y1, double t_z1,
					double q_w1, double q_x1, double q_y1, double q_z1,
					double t_var1, double q_var1)
		: t_x0(t_x0), t_y0(t_y0), t_z0(t_z0),
		  q_w0(q_w0), q_x0(q_x0), q_y0(q_y0), q_z0(q_z0),
		  t_var0(t_var0), q_var0(q_var0),
		  t_x1(t_x1), t_y1(t_y1), t_z1(t_z1),
		  q_w1(q_w1), q_x1(q_x1), q_y1(q_y1), q_z1(q_z1),
		  t_var1(t_var1), q_var1(q_var1) {}

	template <typename T>
	bool operator()(const T* const w_qwi, const T* w_twi, T* residuals) const
	{
		
		T q_wi[4];
		q_wi[0] = T(q_w1);
		q_wi[1] = T(q_x1);
		q_wi[2] = T(q_y1);
		q_wi[3] = T(q_z1);
		T t_wi[3];
		t_wi[0] = T(t_x1);
		t_wi[0] = T(t_y1);
		t_wi[0] = T(t_z1);
		T wqwi[4];
		T wtwi[3];

		ceres::QuaternionRotatePoint(w_qwi, t_wi, wtwi);
		wtwi[0] += w_twi[0];
		wtwi[1] += w_twi[1];
		wtwi[2] += w_twi[2];

		ceres::QuaternionProduct(w_qwi, q_wi, wqwi);

		residuals[0] = (wtwi[0] - T(t_x0)) / T(t_var0);
		residuals[1] = (wtwi[1] - T(t_y0)) / T(t_var0);
		residuals[2] = (wtwi[2] - T(t_z0)) / T(t_var0);

		T relative_q[4];
		relative_q[0] = T(q_w0);
		relative_q[1] = T(q_x0);
		relative_q[2] = T(q_y0);
		relative_q[3] = T(q_z0);

		T relative_q_inv[4];
		QuaternionInverse(relative_q, relative_q_inv);

		T error_q[4];
		ceres::QuaternionProduct(relative_q_inv, wqwi, error_q); 

		residuals[3] = T(2) * error_q[1] / T(q_var0);
		residuals[4] = T(2) * error_q[2] / T(q_var0);
		residuals[5] = T(2) * error_q[3] / T(q_var0);

		return true;
	}

	static ceres::CostFunction* Create(double t_x0, double t_y0, double t_z0,
					double q_w0, double q_x0, double q_y0, double q_z0,
					double t_var0, double q_var0,
					double t_x1, double t_y1, double t_z1,
					double q_w1, double q_x1, double q_y1, double q_z1,
					double t_var1, double q_var1)
	{
		return (new ceres::AutoDiffCostFunction<RelativeCRTError2, 6, 4, 3>(
			new RelativeCRTError2(t_x0, t_y0, t_z0, q_w0, q_x0, q_y0, q_z0, t_var0, q_var0, 
								  t_x1, t_y1, t_z1, q_w1, q_x1, q_y1, q_z1, t_var1, q_var1)
		));
	}
	
	double t_x0, t_y0, t_z0, t_norm0;
	double q_w0, q_x0, q_y0, q_z0;
	double t_var0, q_var0;
	double t_x1, t_y1, t_z1, t_norm1;
	double q_w1, q_x1, q_y1, q_z1;
	double t_var1, q_var1;
};
}