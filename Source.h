#include <Eigen/Dense>
using namespace Eigen;

class Source {
	
	

	public:
	Vector3d location = Vector3d().Zero();
	double amount = 0;

	Source(Vector3d p, float a) {
		location = p;
		amount = a;
	}
};