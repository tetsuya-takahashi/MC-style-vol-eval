# MC-style-vol-eval

A single header-only C++17 library (depending on Eigen: https://eigen.tuxfamily.org/index.php?title=Main_Page) to evaluate volumes and surfaces of 3D grid-based level-set data in the Marching-Cubes-style. Corresponding Marching-Square-style area and perimeter evaluation is also included.

"Fast Marching-Cubes-Style Volume Evaluation for Level Set Surfaces", Journal of Computer Graphics Techniques (https://www.jcgt.org/published/0011/02/02/)

Tetsuya Takahashi (Adobe),
Christopher Batty (University of Waterloo)

# Usage
- Please include Eigen to use the library.

```
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
  
```
