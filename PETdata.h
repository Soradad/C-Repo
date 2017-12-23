#pragma once

#include<iostream>
namespace pet_sd {


	typedef struct raw_data
	{
		unsigned short m_sBDMID;

		unsigned short m_sCrystalID;

		double m_dthreshold_abs_time[8];
	}raw_data;



	class Pet_data
	{
	public:

		unsigned short m_sBDMID;

		unsigned short m_sCrystalID;

		double m_dthreshold_abs_time[8];

		Pet_data() {};

		Pet_data(raw_data a);

		void show_abs_time();

	};

	class Pet_singles:public Pet_data
	{public:
		double m_dstart_time;

		double m_dthreshold_relative_time[8];


		double m_linear_ringsing_model_para_k;


		double m_exponent_falling_model_para_t0;

		double m_exponent_falling_model_para_k;

		double m_model_peak_relative_time_tc;
		
		double m_energy;


		void show_relative_time();

		void show_parameter();

		Pet_singles(raw_data a);

	};


}
