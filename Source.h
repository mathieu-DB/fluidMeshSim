/*
* Author : Mathieu David-Babin
*/
#include <Eigen/Dense>
using namespace Eigen;

class Source {
	
	public:
	Vector3d location = Vector3d().Zero();
	double amount = 0;
	int vertex = 0;

	Source(int v, float a) {
		vertex = v;
		amount = a;
	}
};