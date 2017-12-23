
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<string>
#include"PETdata.h"

using namespace pet_sd;
using namespace std;
union trans{
	raw_data m_transfered_data;
	char m_untransfered_data[2 + 2 + 8 * 8];
};

int main()
{
	
	raw_data a;
	ifstream inf;
	ofstream outf;

	inf.open("./BDM136.samples", ios_base::binary);
	outf.open("./BDM136.singles", ios_base::binary|ios_base::trunc);

	

	for (int i=0; i<5000;i++)
	{


		
		inf.read((char*)(&a.m_sBDMID), sizeof(a.m_sBDMID));
		inf.read((char*)(&a.m_sCrystalID), sizeof(a.m_sCrystalID));
		inf.read((char*)(&a.m_dthreshold_abs_time), sizeof(double) * 8);

		Pet_singles b(a);
		

		outf.write((char*)(&b.m_energy), sizeof(double));
		
		{
			cout << i <<" singles has completed"<< endl;
		}

		//b.show_relative_time();
		//b.show_parameter();
	}

	inf.close();
	outf.close();

	system("PAUSE");
	return 0;
}