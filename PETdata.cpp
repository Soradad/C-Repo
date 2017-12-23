#include <iostream>
#include<fstream>
#include<math.h>
#include"Eigen/Eigen"
#include"PETdata.h"

#define VOLTAGE_THRESHOLD_1 40
#define VOLTAGE_THRESHOLD_2 110
#define VOLTAGE_THRESHOLD_3 180
#define VOLTAGE_THRESHOLD_4 270


#define LN_VOLTAGE_THRESHOLD_1 log(VOLTAGE_THRESHOLD_1)
#define LN_VOLTAGE_THRESHOLD_2 log(VOLTAGE_THRESHOLD_2)
#define LN_VOLTAGE_THRESHOLD_3 log(VOLTAGE_THRESHOLD_3)
#define LN_VOLTAGE_THRESHOLD_4 log(VOLTAGE_THRESHOLD_4)

#define TC_PRECISION 1

using namespace std;
using namespace Eigen;

namespace pet_sd{


	Pet_data::Pet_data(raw_data raw)
	{
		using namespace std;

		//copy raw data from structure raw_data

		for (int i = 0; i < 8; i++)
		{
			Pet_data::m_dthreshold_abs_time[i] = raw.m_dthreshold_abs_time[i];
		}
		Pet_data::m_sBDMID = raw.m_sBDMID;
		Pet_data::m_sCrystalID = raw.m_sCrystalID;

	}

	void Pet_data::show_abs_time()
	{
		using namespace std;

		for (int i = 0; i < 8; i++)
		{
			std::cout << m_dthreshold_abs_time[i] << endl;
		}
	}


	Pet_singles::Pet_singles(raw_data raw)
	{
		
		using namespace std;
		using namespace Eigen;

		//copy raw data from structure raw_data

		for (int i = 0; i < 8; i++)
		{
			Pet_singles::m_dthreshold_abs_time[i] = raw.m_dthreshold_abs_time[i];
		}
		Pet_singles::m_sBDMID = raw.m_sBDMID;
		Pet_singles::m_sCrystalID = raw.m_sCrystalID;



		//double quadratic_const_H;//quadratic falling model parameter,voltage=H-k*sqrt(t-start_time) which are not used
		//double quadratic_const_k;


		double exp_const_t0;//exponent falling model parameter, voltage=exp(t0-tn)


		double lin_const_k;
		double lin_const_h;//linear rising model parameter,voltage=k*(tn-ts)=k*tn+h

		double lin_const_abs_ts;
		double lin_const_rela_ts;
		double lin_exp_const_rela_tc;


		//double lg_vol_thresh_1 = log(VOLTAGE_THRESHOLD_1);
		//double lg_vol_thresh_2 = log(VOLTAGE_THRESHOLD_2);
		//double lg_vol_thresh_3 = log(VOLTAGE_THRESHOLD_3);
		//double lg_vol_thresh_4 = log(VOLTAGE_THRESHOLD_4);

		MatrixXd _Mat_linear_rising_model_A(4, 2);//matrix of [1,t1;1,t2;1,t3;1,t4]

		MatrixXd _Mat_linear_rising_model_h_k(2, 1);//matrix of [h;k]

		MatrixXd _Mat_exponent_falling_model_A(4, 2);//matrix if[1,ln(t1);......]

		MatrixXd _Mat_exponent_falling_model_tn(4, 1);//matrix if[1,ln(t1);......]

		MatrixXd _Mat_exponent_falling_model_t0_reciprocal_k(2, 1);


		MatrixXd _Mat_of_THRESHOLD(4, 1);//matrix of 4 thresholds

		MatrixXd _Mat_of_ln_THRESHOLD(4, 1);//matrix of 4 thresholds

		_Mat_of_THRESHOLD << VOLTAGE_THRESHOLD_1, VOLTAGE_THRESHOLD_2, VOLTAGE_THRESHOLD_3, VOLTAGE_THRESHOLD_4;

		

		double first_thresh_time = Pet_singles::m_dthreshold_abs_time[0];//make the first threshold time the benchmark


		for (int j = 0; j < 8; j++)
		{//make time relative to the first threshold time

			Pet_singles::m_dthreshold_relative_time[j] = Pet_singles::m_dthreshold_abs_time[j] - first_thresh_time;

		}






		_Mat_linear_rising_model_A << 1, Pet_singles::m_dthreshold_relative_time[0],
			1, Pet_singles::m_dthreshold_relative_time[1],
			1, Pet_singles::m_dthreshold_relative_time[2],
			1, Pet_singles::m_dthreshold_relative_time[3];


		_Mat_linear_rising_model_h_k = ((_Mat_linear_rising_model_A.transpose() * _Mat_linear_rising_model_A).inverse())*(_Mat_linear_rising_model_A.transpose()*_Mat_of_THRESHOLD);//calculating the linear rising model parameter, x=(ATA)'ATb,x=h_k,A=t,b=THRESHOLD


		lin_const_h = _Mat_linear_rising_model_h_k(0, 0);


		Pet_singles::m_linear_ringsing_model_para_k		=	lin_const_k		=	_Mat_linear_rising_model_h_k(1, 0);



		lin_const_rela_ts = -(lin_const_h / lin_const_k);//calculate the relative_ts which relative to first_thresh_time

		Pet_singles::m_dstart_time	=	lin_const_abs_ts	=	first_thresh_time + lin_const_rela_ts;


		for (int j = 0; j < 8; j++)
		{//make abs_ts the benchmark ,which we define it the arriving time of this pulse

			Pet_singles::m_dthreshold_relative_time[j] = Pet_singles::m_dthreshold_abs_time[j] - Pet_singles::m_dstart_time;

		}






		_Mat_exponent_falling_model_tn << Pet_singles::m_dthreshold_relative_time[4], Pet_singles::m_dthreshold_relative_time[5], Pet_singles::m_dthreshold_relative_time[6], Pet_singles::m_dthreshold_relative_time[7];

		_Mat_exponent_falling_model_A << 1, -LN_VOLTAGE_THRESHOLD_4,
			1, -LN_VOLTAGE_THRESHOLD_3,
			1, -LN_VOLTAGE_THRESHOLD_2,
			1, -LN_VOLTAGE_THRESHOLD_1;

		_Mat_exponent_falling_model_t0_reciprocal_k=((_Mat_exponent_falling_model_A.transpose() * _Mat_exponent_falling_model_A).inverse())*(_Mat_exponent_falling_model_A.transpose()*_Mat_exponent_falling_model_tn);

		Pet_singles::m_exponent_falling_model_para_t0 = _Mat_exponent_falling_model_t0_reciprocal_k(0, 0);

		Pet_singles::m_exponent_falling_model_para_k = 1/_Mat_exponent_falling_model_t0_reciprocal_k(1, 0);

		double lin_exp_const_rela_tc_lower_bound = Pet_singles::m_dthreshold_relative_time[3];
		double lin_exp_const_rela_tc_upper_bound = Pet_singles::m_dthreshold_relative_time[4];
		double error_crossing_y = INFINITY;

		for (; error_crossing_y > TC_PRECISION;)
		{
			lin_exp_const_rela_tc = (lin_exp_const_rela_tc_lower_bound + lin_exp_const_rela_tc_upper_bound) / 2;
			error_crossing_y = exp(m_exponent_falling_model_para_k*(m_exponent_falling_model_para_t0- lin_exp_const_rela_tc))-m_linear_ringsing_model_para_k*lin_exp_const_rela_tc;
			if (error_crossing_y > 0)//given that the function is monotone decreasing
			{
				lin_exp_const_rela_tc_lower_bound = lin_exp_const_rela_tc;
			}
			else
			{
				lin_exp_const_rela_tc_upper_bound = lin_exp_const_rela_tc;
			}
		}
		Pet_singles::m_model_peak_relative_time_tc = lin_exp_const_rela_tc;



		


		





		Pet_singles::m_energy = exp(m_exponent_falling_model_para_k *( m_exponent_falling_model_para_t0- m_model_peak_relative_time_tc) ) / m_exponent_falling_model_para_k
			+0.5*(m_exponent_falling_model_para_k * m_model_peak_relative_time_tc * m_model_peak_relative_time_tc);
		
		//it's a math regression of minimum average integration error with the exponent falling model parameter
		  


	}



	

	void Pet_singles::show_relative_time()
	{
		using namespace std;

		std::cout << "start_time = " << m_dstart_time << std::endl;
		for (int i = 0; i < 8; i++)
		{
			std::cout << m_dthreshold_relative_time[i]<< endl;
		}

	}

	void Pet_singles::show_parameter()
	{
		using namespace std;
		std::cout << "m_linear_ringsing_model_para_k = "<< m_linear_ringsing_model_para_k << endl;
		std::cout << "m_exponent_falling_model_para_t0 = " << m_exponent_falling_model_para_t0 << endl;
		std::cout << "m_exponent_falling_model_para_k = " << m_exponent_falling_model_para_k << endl;
		std::cout << "m_energy = " << m_energy << endl;

	}

};




