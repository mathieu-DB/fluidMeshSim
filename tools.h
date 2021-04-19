/**
*	Author : Mathieu David-Babin
*   License ...
* 
*   Set of methods for manipulating a range of objects
*/
#ifndef __tools_h__
#define __tools_h__

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <map>
#include <list>
#include <vector>
#include <utility>

using namespace Eigen;

namespace tools {
	typedef Eigen::Triplet<double> T;

	inline void RotateVector(Vector3d& v, Vector3d axis, Vector3d n1, Vector3d n2) {
		/*double c = (1 - n1.dot(n2));
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

		v = R * v;*/
	}

	inline void GetNormalizedTriangleNormal(Vector3d v0, Vector3d v1, Vector3d v2, Vector3d& n) {
		n = ((v2 - v0).cross(v1 - v0)).normalized();
	}

	//void SetVectorToZero(std::vector<double> v) {
	//	int s = v.size();
	//	for (int i = 0; i < s; i++) {
	//		v[i] = 0;
	//	}
	//}

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
	inline void CreateFusedMesh(MatrixXd V_split, MatrixXd& V_fused, MatrixXi F_split, MatrixXi& F_fused, std::map<int, int>& V_mapSF, std::map<int,int>& V_mapFS, std::map<int,int>& dupe_map, VectorXi bnd) {
		// 1. Find the dupes
		std::vector<int> l = {};
		std::vector<std::pair<int, int>> dupes;
		int c = 0;
		std::vector<int> nbnd;		//indices of vectors not placed on the seam

		//Fill nbnd
		bool found = false;
		for (int i = 0; i < V_split.rows(); i++) {
			for (int j = 0; j < bnd.size(); j++) {
				if (i == bnd[j]) {
					found = true;
				}
			}
			if (!found) {
				nbnd.emplace_back(i);
			}
			found = false;
		}

		//Fill l with sequence from 0 to |bnd|
		for (int i = 0; i < bnd.size(); i++) {
			l.emplace_back(i);
			
		}
		//Use l to identify the indices in bnd that have not been identified as dupes
		while (!l.empty()) {
			Vector3d v;
			Vector3d v2;
			v= V_split.row(bnd[l[0]]);
			int s = l.size();
			for (int i = 1; i < s; i++) {
				v2 = V_split.row(bnd[l[i]]);
				if (v.isApprox(v2)) {
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
		V_fused.resize(size,V_split.cols());
		for (int i : nbnd) {
			V_mapSF[i] = c;
			V_mapFS[c] = i;
			V_fused.row(c) = V_split.row(i);
			c++;
		}

		for (std::pair<int,int> p : dupes) {
			V_mapSF[p.first] = c;
			if (p.second >= 0) {
				V_mapSF[p.second] = c;
			}
			V_mapFS[c] = p.first;
			V_fused.row(c) = V_split.row(p.first);
			c++;
		}

		//Getting the faces with the update indexes
		F_fused.resize(F_split.rows(), F_split.cols());
		for (int i = 0; i < F_split.rows(); i++) {
			F_fused.row(i) = Eigen::Vector3i(V_mapSF[F_split.row(i)[0]], V_mapSF[F_split.row(i)[1]], V_mapSF[F_split.row(i)[2]]);
		}


	}

	inline void CreateSparseWeightMatrix(MatrixXd V, MatrixXi F, Eigen::SparseMatrix<double>& S) {
		int a, b, c;
		double dab, dac, dbc;
		int numF = F.rows();
		int numV = V.rows();
		S.resize(numV, numV);
		std::vector<T> triplets;
		for (int i = 0; i < numF; i++) {
			a = F.row(i)[0];
			b = F.row(i)[1];
			c = F.row(i)[2];

			dab = (V.row(a) - V.row(b)).norm();
			dac = (V.row(a) - V.row(c)).norm();
			dbc = (V.row(b) - V.row(c)).norm();

			triplets.emplace_back(a, a, dab + dac);
			triplets.emplace_back(b, b, dab + dbc);
			triplets.emplace_back(c, c, dac + dbc);

			triplets.emplace_back(a, c, dac);
			triplets.emplace_back(c, a, dac);
			triplets.emplace_back(a, b, dab);
			triplets.emplace_back(b, a, dab);
			triplets.emplace_back(b, c, dbc);
			triplets.emplace_back(c, b, dbc);
		}

		S.setFromTriplets(triplets.begin(), triplets.end());
	}

	inline double ComputeMeshSurfaceArea(MatrixXd V, MatrixXi F) {
		double r = 0;

		int a, b, c;
		if (V.cols() == 3) {
			Vector3d ab;
			Vector3d ac;
			for (int i = 0; i < F.rows(); i++) {
				a = F.row(i)[0];
				b = F.row(i)[1];
				c = F.row(i)[2];
				ab = V.row(a) - V.row(b);
				ac = V.row(a) - V.row(c);
				r += (ab.cross(ac)).norm() / 2;
			}
		}
		else {
			VectorXd ab;
			VectorXd ac;
			for (int i = 0; i < F.rows(); i++) {
				a = F.row(i)[0];
				b = F.row(i)[1];
				c = F.row(i)[2];

				ab = V.row(a) - V.row(b);
				ac = V.row(a) - V.row(c);

				r += std::abs((ab[0] * ac[1] - ac[0] * ab[1]) / 2);
			}
		}
		return r;
	}
}

#endif