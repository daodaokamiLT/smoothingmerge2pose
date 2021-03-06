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

struct RTError
{
	RTError(double t_x, double t_y, double t_z, 
			double q_w, double q_x, double q_y, double q_z, double t_var, double q_var)
				  :t_x(t_x), t_y(t_y), t_z(t_z), q_w(q_w), q_x(q_z), q_y(q_y), q_z(q_z), t_var(t_var), q_var(q_var){}

	template <typename T>
	bool operator()(const T* w_q_j, const T* tj, T* residuals) const
	{
		residuals[0] = (tj[0] - T(t_x)) / T(t_var);
		residuals[1] = (tj[1] - T(t_y)) / T(t_var);
		residuals[2] = (tj[2] - T(t_z)) / T(t_var);
		
		T q[4];
		q[0] = T(q_w);
		q[1] = T(q_x);
		q[2] = T(q_y);
		q[3] = T(q_z);
		
		T q_inv[4];
		QuaternionInverse(q, q_inv);

		T error_q[4];
		ceres::QuaternionProduct(q_inv, w_q_j, error_q);

		residuals[3] = T(2) * error_q[1] / T(q_var);
		residuals[4] = T(2) * error_q[2] / T(q_var);
		residuals[5] = T(2) * error_q[3] / T(q_var);

		return true;
	}

	static ceres::CostFunction* Create(const double t_x, const double t_y, const double t_z,
									   const double q_w, const double q_x, const double q_y, const double q_z,
									   const double t_var, const double q_var) 
	{
	  return (new ceres::AutoDiffCostFunction<RTError, 6, 4, 3>(
	          	new RTError(t_x, t_y, t_z, q_w, q_x, q_y, q_z, t_var, q_var)));
	}

	double t_x, t_y, t_z, q_w, q_x, q_y, q_z, t_var, q_var;

};

struct ChangeCoordinateError{
	ChangeCoordinateError(const double t0_x, const double t0_y,const double t0_z,
	const double q0_w,const  double q0_x,const double q0_y,const double q0_z,
	const double t0_var,const double q0_var,
	const double t1_x, const double t1_y,const double t1_z,
	const double q1_w,const  double q1_x,const double q1_y,const double q1_z,
	const double t1_var,const double q1_var):t0_x(t0_x), t0_y(t0_y), t0_z(t0_z), 
			q0_w(q0_w), q0_x(q0_x), q0_y(q0_y), q0_z(q0_z), 
			t0_var(t0_var), q0_var(q0_var),
			t1_x(t1_x), t1_y(t1_y), t1_z(t1_z), 
			q1_w(q1_w), q1_x(q1_x), q1_y(q1_y), q1_z(q1_z), 
			t1_var(t1_var), q1_var(q1_var){}

	template <typename T>
	bool operator()(const T* const q_i_j, const T* t_i_j, T* residuals) const
	{
		T qwi_0[4]; 
		qwi_0[0] = T(q0_w);qwi_0[1] = T(q0_x);qwi_0[2] = T(q0_y);qwi_0[3] = T(q0_z);
		T twi_0[3]; 
		twi_0[0] = T(t0_x);twi_0[1] = T(t0_y);twi_0[2] = T(t0_z);

		T qwj_1[4]; 
		qwj_1[0] = T(q1_w);qwj_1[1] = T(q1_x);qwj_1[2] = T(q1_y);qwj_1[3] = T(q1_z);
		T twj_1[3]; 
		twj_1[0] = T(t1_x);twj_1[1] = T(t1_y);twj_1[2] = T(t1_z);

		T qwj_0[4];
		ceres::QuaternionProduct(qwi_0, q_i_j, qwj_0);

		T twj_0[3];
		ceres::QuaternionRotatePoint(qwi_0, t_i_j, twj_0);
		twj_0[0] += twi_0[0];
		twj_0[1] += twi_0[1];
		twj_0[2] += twi_0[2];

		T tmp[3];
		tmp[0] = (twj_1[0] - twj_0[0]);
		tmp[1] = (twj_1[1] - twj_0[1]);
		tmp[2] = (twj_1[2] - twj_0[2]);

		T delta_t_01[3];
		ceres::QuaternionRotatePoint(qwj_0, tmp, delta_t_01);
		residuals[0] = delta_t_01[0];
		residuals[1] = delta_t_01[1];
		residuals[2] = delta_t_01[2];

		T qwj_0_inv[4];
		QuaternionInverse(qwj_0, qwj_0_inv);

		T err_q01[4];
		ceres::QuaternionProduct(qwj_0_inv, qwj_1, err_q01);

		residuals[3] = T(2) * err_q01[1];
		residuals[4] = T(2) * err_q01[2];
		residuals[5] = T(2) * err_q01[3];
		
		return true;
	}

	static ceres::CostFunction* Create(const double t0_x, const double t0_y,const double t0_z,
	const double q0_w,const  double q0_x,const double q0_y,const double q0_z,
	const double t0_var,const double q0_var,
	const double t1_x, const double t1_y,const double t1_z,
	const double q1_w,const  double q1_x,const double q1_y,const double q1_z,
	const double t1_var,const double q1_var) 
	{
	  return (new ceres::AutoDiffCostFunction<
	          ChangeCoordinateError, 6, 4, 3>(
	          	new ChangeCoordinateError(t0_x, t0_y, t0_z, q0_w, q0_x, q0_y, q0_z, t0_var, q0_var, 
				  			t1_x, t1_y, t1_z, q1_w, q1_x, q1_y, q1_z, t1_var, q1_var)));
	}

	double t0_x, t0_y, t0_z, q0_w, q0_x, q0_y, q0_z, t0_var, q0_var;
	double t1_x, t1_y, t1_z, q1_w, q1_x, q1_y, q1_z, t1_var, q1_var;
};

struct ChangeCoordinateTError{
	ChangeCoordinateTError(
		const double t0_x, const double t0_y,const double t0_z,
		const double t0_var,
		const double t1_x, const double t1_y,const double t1_z,
		const double t1_var)
		   :t0_x(t0_x), t0_y(t0_y), t0_z(t0_z), 
			t0_var(t0_var),
			t1_x(t1_x), t1_y(t1_y), t1_z(t1_z), 
			t1_var(t1_var){}

	// qvec is  x, y, z, w ???
	// but real w, x, y, z ???
	template <typename T>
	bool operator()(const T* const q_i_j, const T* t_i_j, T* residuals) const
	{
		T twi_0[3]; 
		twi_0[0] = T(t0_x);twi_0[1] = T(t0_y);twi_0[2] = T(t0_z);

		T twj_1[3]; 
		twj_1[0] = T(t1_x);twj_1[1] = T(t1_y);twj_1[2] = T(t1_z);

		T twj_0[3];
		ceres::QuaternionRotatePoint(q_i_j, twi_0, twj_0);
		twj_0[0] += t_i_j[0];
		twj_0[1] += t_i_j[1];
		twj_0[2] += t_i_j[2];

		residuals[0] = (twj_1[0] - twj_0[0]);
		residuals[1] = (twj_1[1] - twj_0[1]);
		residuals[2] = (twj_1[2] - twj_0[2]);

		return true;
	}

	static ceres::CostFunction* Create(
		const double t0_x, const double t0_y,const double t0_z,
		const double t0_var,
		const double t1_x, const double t1_y,const double t1_z,
		const double t1_var) 
	{
	  return new ceres::AutoDiffCostFunction<ChangeCoordinateTError, 3, 4, 3>(
	          		new ChangeCoordinateTError(t0_x, t0_y, t0_z, t0_var, t1_x, t1_y, t1_z, t1_var)
				);
	}

	double t0_x, t0_y, t0_z, t0_var;
	double t1_x, t1_y, t1_z, t1_var;
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

struct ChangeCoordinateRTError{
	ChangeCoordinateRTError(const double t_x, const double t_y, const double t_z, 
				const double q_w, const double q_x, const double q_y, const double q_z, 
				const double t_var, const double q_var):
				t_x(t_x), t_y(t_y), t_z(t_z), t_var(t_var), 
				q_w(q_w), q_x(q_x), q_y(q_y), q_z(q_z), q_var(q_var){}

	template <typename T>
	bool operator()(const T* const q_wi_wj, const T* t_wi_wj, 
					const T* q_wj, const T* t_wj, T* residuals) const
	{
		T q_wi[4];
		q_wi[0] = T(q_w);
		q_wi[1] = T(q_x);
		q_wi[2] = T(q_y);
		q_wi[3] = T(q_z);

		T t_wi[3];
		t_wi[0] = T(t_x);
		t_wi[1] = T(t_y);
		t_wi[2] = T(t_z);

		T q_wij[4];
		ceres::QuaternionProduct(q_wi_wj, q_wj, q_wij);
		T t_wij[3];
		ceres::QuaternionRotatePoint(q_wi_wj, t_wj, t_wij);
		t_wij[0] += t_wi_wj[0];
		t_wij[1] += t_wi_wj[1] ;
		t_wij[2] += t_wi_wj[2] ;
		residuals[0] = (t_wij[0] - t_wi[0])/t_var;
		residuals[1] = (t_wij[1] - t_wi[1])/t_var;
		residuals[2] = (t_wij[2] - t_wi[2])/t_var;
		T q_wij_inv[4], delta_ij[4];
		QuaternionInverse(q_wij, q_wij_inv);
		ceres::QuaternionProduct(q_wij_inv, q_wi, delta_ij);
		residuals[3] = T(2) * delta_ij[1] / q_var;
		residuals[4] = T(2) * delta_ij[2] / q_var;
		residuals[5] = T(2) * delta_ij[3] / q_var;

		return true;
	}


	static ceres::CostFunction* Create(const double t_x, const double t_y,const double t_z,
			const double q_w,const  double q_x,const double q_y,const double q_z,
			const double t_var,const double q_var) 
	{
	  return (new ceres::AutoDiffCostFunction<
	          ChangeCoordinateRTError, 6, 4, 3, 4, 3>(
	          	new ChangeCoordinateRTError(t_x, t_y, t_z, q_w, q_x, q_y, q_z, t_var, q_var)));
	}
	
	double t_x, t_y, t_z, t_var;
	double q_w, q_x, q_y, q_z, q_var;
};
}