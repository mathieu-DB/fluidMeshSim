/**
*	Author : Mathieu David-Babin
*   License ...
* 
*   Set of methods for manipulating a range of objects
*/

#include <Eigen/Dense>
using namespace Eigen;

namespace tools {

	void RotateVector(Vector3d& v, Vector3d axis, Vector3d n1, Vector3d n2) {
		double c = (1 - n1.dot(n2));
		double s = (n1.cross(n2)).norm();
		Matrix3d R;
		R.setZero();
		double x = axis.x();
		double y = axis.y();
		double z = axis.z();

		R(0, 0) = c * x * x;
		R(0, 1) = c * x * y - s * z;
		R(0, 2) = c * x * z - s * y;
		R(1, 0) = R.row(0)[1];
		R(1, 1) = c * y * y;
		R(1, 2) = c * y * z - s * x;
		R(2, 0) = R.row(0)[2];
		R(2, 1) = R.row(1)[2];
		R(2, 2) = c * z * z;

		v = R * v;
	}

	void GetNormalizedTriangleNormal(Vector3d v0, Vector3d v1, Vector3d v2, Vector3d& n) {
		n = ((v2 - v0).cross(v1 - v0)).normalized();
	}


}