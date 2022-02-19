#include <iostream>
#include "fraction.hpp"

int main()
{
	using namespace Fraction;
	//levelset for (i, j), (i+1, j), (i + 1, j + 1), (i, j + 1)
	std::array<fType, 4> phi2d{ -0.5, 0.5, -0.5, 0.5 };
	std::cout << get_ms_area(phi2d) << std::endl;
	std::cout << get_ms_len(phi2d) << std::endl;
	
	//levelset for (i,j,k), (i+1,j,k), (i+1,j+1,k), (i,j+1,k), 
	//             (i,j,k+1), (i+1,j,k+1), (i+1,j+1,k+1), (i,j+1,k+1)
	std::array<fType, 8> phi3d{ -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5 };
	std::cout << get_mc_vol(phi3d) << std::endl;
	std::cout << get_mc_area(phi3d) << std::endl;
}
