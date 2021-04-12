/**
*	Author : Mathieu David-Babin
*   License ...
* 
*   Set of methods for manipulating a range of objects
*/

#include <Eigen/Dense>
#include <map>
#include <list>

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

	void SetVectorToZero(std::vector<double> v) {
		int s = v.size();
		for (int i = 0; i < s; i++) {
			v[i] = 0;
		}
	}

	/*
	* Takes a mesh with duplicate vertices along the seam and creates a new one by fusing both sides of the seam
	* NOT OPTIMAL !!!
	* 
	* @param V_split		The original vertices (with duplicates)
	* @param V_fused		he empty matrix to receive the fused vertices
	* @param F_split		The original faces (using duplicate vertices)
	* @param F_fused		The empty matrix to receive the newly defined faces
	* @param V_mapSF		Map that will hold the relation between the old vertices and the new ones [old] --> [new]
	* @param dupe_map	Map that will hold the relation between the duplicate vertices
	* @param bnd			List of the old vertices placed along the seam
	*/
	void CreateFusedMesh(MatrixXd V_split, MatrixXd& V_fused, MatrixXi F_split, MatrixXi& F_fused, std::map<int, int>& V_mapSF, std::map<int,int> V_mapFS, std::map<int,int>& dupe_map, VectorXi bnd) {
		// 1. Find the dupes
		std::vector<int> l = {};
		std::vector<pair<int, int>> dupes;
		int c = 0;
		std::vector<int> nbnd;		//indices of vectors not placed on the seam

		//Fill nbnd
		for (int i = 0; i < V_split.rows(); i++) {
			for (int j = 0; j < bnd.size(); j++) {
				if (i == bnd[j]) {
					nbnd.emplace_back(i);
					j = bnd.size();
				}
			}
		}

		//Fill l with sequence from 0 to |bnd|
		for (int i = 0; i < bnd.size(); i++) {
			l.emplace_back(i);
			
		}
		//Use l to identify the indices in bnd that have not been identified as dupes
		while (!l.empty()) {
			Vector3d v = V_split.row(bnd[l[0]]);
			int s = l.size();
			for (int i = 1; i < s; i++) {
				if (v == V_split.row(bnd[l[i]])) {
					dupes.emplace_back(bnd[l[0]], bnd[l[i]]);
					dupe_map[bnd[l[0]]] = bnd[l[i]];
					dupe_map[bnd[l[i]]] = bnd[l[0]];
					l.erase(l.begin() + i);
					l.erase(l.begin());
					i = s;
				}
			}
			//No duplicate of the first element found
			if (s == l.size()) {
				dupes.emplace_back(bnd[l[0]], -1);
				l.erase(l.begin());
			}
		}
		int size = dupes.size() + nbnd.size();
		V_fused.resize(size, Eigen::NoChange);
		for (int i : nbnd) {
			V_mapSF[i] = c;
			V_mapFS[c] = i;
			V_fused.row(c) = V_split.row(i);
			c++;
		}

		for (pair<int,int> p : dupes) {
			V_mapSF[p.first] = c;
			V_mapSF[p.second] = c;
			V_mapFS[c] = p.first;
			V_fused.row(c) = V_split.row(p.first);
			c++;
		}

		//Getting the faces with the update indexes
		for (int i = 0; i < F_split.rows(); i++) {
			F_fused.row(i) = Eigen::VectorXi(V_mapSF[F_split.row(i)[0]], V_mapSF[F_split.row(i)[1]], V_mapSF[F_split.row(i)[2]]);
		}


	}
}