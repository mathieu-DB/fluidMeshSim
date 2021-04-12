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
#include <map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Source.h"
#include "trimesh.h"

using namespace std;
	class Fluid {
	private:
		const static int DIM = 3;
		int N = 16;
		double dx = 1;

		int faces;
		
		trimesh::trimesh_t mesh;
		vector<double> div;
		vector<double> p;

		
		/**
		* Velocity Field non-staggered, packed. Each velocity vector is set 
		* at a vertex of the mesh
		*/
		vector<double> U0[DIM];

		/** Temporary velocity field*/
		vector<double> U1[DIM];

		/** Temperature (packed)*/
		vector<double> temperature0;

		/** Temporary Temperature (packed)*/
		vector<double> temperature1;

		/** Array of vertices for the mesh*/
		Eigen::MatrixXd V;

		/** Array of triangles and their associated vertex indexes for the mesh*/
		Eigen::MatrixXi F;

		/* Altered version of the V matrix that had the duped vertices along the seam fused*/
		Eigen::MatrixXd V_fused;

		/* Altered version of the F matrix with vertex indexes corresponding to the ones in V_fused*/
		Eigen::MatrixXi F_fused;

		/*Map containing the relation from the vertices in V_fused [V] --> [V_fused]*/
		map<int, int> V_mapSF;

		/*Map containing the relation from the vertices in V_fused [V_fused] --> [V]*/
		map<int, int> V_mapFS;

		/*Map associating dupe vertices in V with one and the other*/
		map<int, int> dupe_map;
			

		/** Array of the vertices in 2D parametrized form*/
		Eigen::MatrixXd uv;

		Eigen::SparseMatrix<double> L;

		Eigen::SparseMatrix<double> L_split;

		/** Time elapsed in the fluid simulation*/
		double elapsed;

		/** Sources of heat and cold*/
		std::list<Source> sources;
	public:

		double doubleVariable = 0;
		float viscosity = 1e-6;
		float diffusion = 1e-6;
		float bouyancy = 0.1;
		int iterations = 30;
		float timeStep = 0.1;
		float gravity = 1;


		bool velocityDiffuse = true;
		bool velocityProject = true;
		bool velocityAdvect = true;
		bool scalarDiffuse = true;
		bool scalarAdvect = true;


		Fluid(trimesh::trimesh_t mesh, Eigen::MatrixXd v, Eigen::MatrixXi f, Eigen::MatrixXd uv, Eigen::MatrixXd v_fused, Eigen::MatrixXi f_fused, map<int,int> v_mapsf, map<int,int> v_mapfs, map<int,int> dupe_map);

		void setup();

		int IX(int i, int j);

		void setBoundary(int b, vector<double> x);

		void getVelocity(Vector3d x,  Vector3d vel);

		void getVelocity(Vector3d x,  vector<double> U[], Vector3d vel);

		double interpolate(Vector3d x,  vector<double> s);

		void traceParticle(Vector3d x0, double h, Vector3d &x1);

		void traceParticleFromVertex(int v, vector<double> U[], double h, Vector3d& x1);

		void traceParticle(Vector3d x0, vector<double> U[], double h, Vector3d &x1);

		void diffuse(vector<double> &S1, vector<double> &S0, int b, double diff, double dt);

		void transport(vector<double> &s1, vector<double> &s0, vector<double> U[], double dt);

		void project(vector<double> U[]);

		void addForce(vector<double> U[], double dt, Vector3d x, Vector3d f, int v);

		void addSource(vector<double> &S, double dt, Vector3d x, double amount, int v);

		double getReferenceTemperature();

		void addTemperatureForce(vector<double> U[], double dt);

		void velocityStep(double dt);

		void scalarStep(double dt);

		void step();

		bool isInTriangle(Vector3d x, int t);

		void createSource(int v, double amount);

		double interpolateTempForVertex(Vector2d x);

		double getMaxTemp();

		double getMinTemp();

		bool isDupe(int v);

		int indentifyTriangle(int startV, Vector3d x, int& t);
	};

#endif // !fluid_h
