/**
* Author: Mathieu David-Babin
* License: ...
* 
* Based on the Java program given for McGill's COMP559 Assignment 3 
* Stable fluids which is itself based on [STAM03] Real-Time Fluid Dynamics for Games
* 
* The original program is designed for a 2D square grid while this version
* tries to extend it to work on a 3D triangle mesh.
* 
*/
#ifndef fluid_h
#define fluid_h

#include <list>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Source.h"
#include "trimesh.h"

using namespace std;
	class fluid {
	private:
		const static int DIM = 3;
		const static int iterations = 10;

		int faces;
		
		trimesh::trimesh_t mesh;



		/**
		* Velocity Field non-staggered, packed. Each velocity vector is set 
		* at a vertex of the mesh
		*/
		vector<double> U0[DIM];

		/** Temporary velocity field*/
		vector<double> U1[DIM];

		/** Temperature (packed)*/
		vector<double> temperature0[DIM];

		/** Temporary Temperature (packed)*/
		vector<double> temperature1[DIM];

		/** Array of vertices for the mesh*/
		Eigen::Matrix3d V;

		/** Array of triangles and their associated vertex indexes for the mesh*/
		Eigen::Matrix3i F;

		Eigen::SparseMatrix<double> L;

		/** Time elapsed in the fluid simulation*/
		double elapsed;

		/** Sources of heat and cold*/
		std::list<Source> sources;
	public:
		fluid(trimesh::trimesh_t mesh, Eigen::Matrix3Xd v, Eigen::Matrix3Xi f);

		void setup();

		int IX(int i, int j);

		void setBoundary(int b, double x[]);

		void getVelocity(Vector3d x, int t, Vector3d vel);

		void getVelocity(Vector3d x, int t, vector<double> U[], Vector3d vel);

		double interpolate(Vector3d x, int t, vector<double> s);

		int traceParticle(Vector3d x0, int t, double h, Vector3d x1);

		int traceParticle(Vector3d x0, int t, vector<double> U[], double h, Vector3d x1);

		void diffuse(vector<double> S1, vector<double> S0, int b, double diff, double dt);

		void transport(vector<double> s1, vector<double> s0, vector<double> U[], double dt);

		void project(vector<double> U[]);

		void addForce(vector<double> U[], double dt, Vector3d x, Vector3d f, int v);

		void addSource(vector<double> S, double dt, Vector3d x, double amount, int v);

		double getReferenceTemperature();

		void addTemperatureForce(vector<double> U[], double dt);

		void velocityStep(double dt);

		void scalarStep(double dt);

		void step();

		bool isInTriangle(Vector3d x, int t);

	};

#endif // !fluid_h
