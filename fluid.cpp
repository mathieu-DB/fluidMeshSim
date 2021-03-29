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
#include "fluid.h"
#include "tools.h"
#include <igl/cotmatrix.h>
using namespace trimesh;

/**
* Constructor
* 
* @param mesh Halfedge mesh structure
* @param v vertices
* @param f faces
*/
fluid::fluid(trimesh::trimesh_t mesh, Eigen::Matrix3Xd v, Eigen::Matrix3Xi f)
{
	mesh = mesh;
	V = v;
	F = f;
	faces = mesh.getNumFaces();
	igl::cotmatrix(V, F, L);
}

/**
* Setups the varius memory structures of the program
* (i.e. U0, U1, temperature0, temperature1)
*/
void fluid::setup()
{
	elapsed = 0;
	for (int i = 0; i < V.rows(); i++) {
		U0[0].emplace_back(0);
		U0[1].emplace_back(0);
		U0[2].emplace_back(0);

		U1[0].emplace_back(0);
		U1[1].emplace_back(0);
		U1[2].emplace_back(0);

		temperature0[0].emplace_back(0);
		temperature0[1].emplace_back(0);
		temperature0[2].emplace_back(0);
		temperature1[0].emplace_back(0);
		temperature1[1].emplace_back(0);
		temperature1[2].emplace_back(0);

	}
}

int fluid::IX(int i, int j)
{
	return 0;
}

void fluid::setBoundary(int b, double x[])
{
}

/**
* Gets the velocity at the given point using interpolation
*
* @param x
* @param vel
*/
void fluid::getVelocity(Vector3d x, int t, Vector3d vel)
{
	getVelocity(x, t, U0, vel);
}

/**
* Gets the velocity in the provided velocity field at the given point using
* interpolation
*
* @param x
* @param U
* @param vel
*/
void fluid::getVelocity(Vector3d x, int t, vector<double> U[], Vector3d vel)
{
	vel[0] = interpolate(x, t, U[0]);
	vel[1] = interpolate(x, t, U[1]);
	vel[2] = interpolate(x, t, U[2]);
}

/**
* Interpolates with the vertices values using barycentric coordinates
* to point x 
* @param x Point needing an interpolation
* @param t Triangle index
* @param s Scalar values to interpolate
* @return interpolated value
*/
double fluid::interpolate(Vector3d x, int t, vector<double> s)
{
	Vector3i vertices = F.row(t);
	Vector3d v0 = V.row(vertices[0]);
	Vector3d v1 = V.row(vertices[1]);
	Vector3d v2 = V.row(vertices[2]);
	double A = ((v1 - v0).cross(v2 - v0)).norm() / 2.0;
	double A0 = ((x - v2).cross(v1 - x)).norm() / 2.0;
	double A1 = ((x - v0).cross(v2 - x)).norm() / 2.0;
	double A2 = ((x - v0).cross(v1 - x)).norm() / 2.0;
	double r = (s[vertices[0]] * A0 + s[vertices[1]] * A1 + s[vertices[2]] * A2) / A;
	return r;
}

/**
* Performs a simple Forward Euler particle trace using the current velocity
* field.
*
* @param x0 Current particle location
* @param h  Time step
* @param x1 Final particle location
*/
int fluid::traceParticle(Vector3d x0, int t, double h, Vector3d x1)
{
	return traceParticle(x0, t, U0, h, x1);
}

/**
* Performs a simple particle trace using Forward Euler. Up to the caller to
* provide a positive or negative time step depending if they want to trace
* forward (i.e., filaments) or backwards (advection term). Note that this
* method should assure that the resulting point be on the mesh.
*
* @param x0 Starting point
* @param t  
* @param U  Velocity field
* @param h  Time step
* @param x1 Resulting point
*/
int fluid::traceParticle(Vector3d x0, int t, vector<double> U[], double h, Vector3d x1)
{
	//TODO Find way to move on the mesh following a given direction
	//TODO Extend so the foward euler can travel across more than 2 triangles
	x1.setZero();
	Vector3d vec(0, 0, 0);
	getVelocity(x0, t, vec);
	vec *= h;
	x1 = x0 + vec;
	if (isInTriangle(x1, t)) {
		return t;
	}
	else {
		trimesh_t::halfedge_t he = mesh.face_halfedge(t);
		while (V.row(he.to_vertex) != x0) {
			he = mesh.next_halfedge(he);
		}
		he = mesh.next_halfedge(mesh.next_halfedge(he));
		Vector3d b0 = V.row(mesh.from_vertex(he));
		Vector3d b = (V.row(he.to_vertex) - b0);
		
		int tt = mesh.get_halfedge(he.opposite_he).face;
		double a0x = x0.x();
		double a0y = x0.y();
		double b0x = b0.x();
		double b0y = b0.y();
		double k = (1 / (b.y() - (b.x() * vec.y() / vec.x()))) * (a0y - b0y + (b0x * vec.y() / vec.x()) + (a0x * vec.y() / vec.x()));
		Vector3d trans = b * k;
		double rest = (vec - (trans - x0)).norm();
		Vector3d w = b.normalized();
		Vector3d n1;
		Vector3d n2;
		tools::GetNormalizedTriangleNormal(V.row(F.row(t)[0]), V.row(F.row(t)[1]), V.row(F.row(t)[2]), n1);
		tools::GetNormalizedTriangleNormal(V.row(F.row(tt)[0]), V.row(F.row(tt)[1]), V.row(F.row(tt)[2]), n2);
		tools::RotateVector(vec, w, n1, n2);
		vec = vec.normalized() * rest;
		x1 = trans + vec;

		return tt;
	}
	return 0;
}

/**
* Diffuse the given scalar field by the given amount. 
*
* @param S1   diffused quantities
* @param S0   initial quantities
* @param diff diffusion coefficient
* @param dt   time step
*/
void fluid::diffuse(vector<double> S1, vector<double> S0, int b, double diff, double dt)
{
	//TODO Remove parameter b
	double diffRate = dt * diff;
	int i, j, k;
	for (k = 0; k < iterations; k++) {
		for (i = 0; i < V.rows(); i++) {
			vector<index_t> vertices = mesh.vertex_vertex_neighbors(i);
			int size = vertices.size();
			double t = 0;
			for (j = 0; j < size; j++) {
				t += S1[vertices[j]];
			}
			t *= diffRate;
			t += S0[i];
			S1[i] = t / (1 + diffRate * size);
		}
	}
}

/**
* Advects / transports scalar quantities by the given velocities.
*
* @param s1 Final advected quantities
* @param s0 Current grid of quantities
* @param U  Velocity field
* @param dt Time step
*/
void fluid::transport(vector<double> s1, vector<double> s0, vector<double> U[], double dt)
{
	Vector3d p1(0, 0, 0);
	Vector3d p2(0, 0, 0);
	int t;
	for (int i = 0; i < V.rows(); i++) {
		p2.setZero();
		p1 = V.row(i);
		t = mesh.vertex_face(i);
		t = traceParticle(p1, t, dt, p2);
		s1[i] = interpolate(p2, t, s0);
		
	}
}

/**
* Does the Poisson solve to make sure that velocities U respect incompressible
* flow
*
* @param U
*/
void fluid::project(vector<double> U[])
{
	/*vector<float> div;
	vector<float> p;
	int v = V.rows();

	for (int i = 0; i < v; i++) {

	}*/
}

/**
* Adds a force at a given point in the provided velocity field.
*
* @param U
* @param dt
* @param x
* @param f
*/
void fluid::addForce(vector<double> U[], double dt, Vector3d x, Vector3d f, int v)
{
	addSource(U[0], dt, x, f.x(), v);
	addSource(U[1], dt, x, f.y(), v);
	addSource(U[2], dt, x, f.z(), v);
}

/**
* Adds some time step scaled amount to the provided scalar field. 
*
* @param S      quantity field to modify
* @param dt     time step
* @param x      position
* @param amount amount
*/
void fluid::addSource(vector<double> S, double dt, Vector3d x, double amount, int v)
{
	//Using the Cotangent Laplacian Matrix to get the vertices' weights
	//Iterate over non-zero coefficients of the sparse laplacian
	//in the column corresponding to vertex v
	//The non-zero coefficients are the vertices adjacent to v
	double div = std::abs(L.coeff(v,v));
	for (SparseMatrix<double>::InnerIterator it(L, v); it; ++it)
	{
		if (it.row() != v) {
			S[it.row()] += amount * dt * it.value() / div;
		}
		it.value();
		it.row();   // row index
		it.col();   // col index (here it is equal to k)
		it.index(); // inner index, here it is equal to it.row()
	}


}

/**
* Gets the average temperature of the continuum. (use this in computing
* buoyancy forces)
*
* @return average temperature
*/
double fluid::getReferenceTemperature()
{
	double t = 0;
	for (int i = 0; i < V.rows(); i++) {
		t += temperature0->at(i);
	}
	return t/V.rows();
}

/**
* Applies buoyancy forces to the velocity field due to temperature. Use Foster
* and Metaxis [1997] Equation 2: F_{bv} = beta g_v (T_0 - T_k)
*
* @param U
* @param dt
*/
void fluid::addTemperatureForce(vector<double> U[], double dt)
{
	double refTemp = getReferenceTemperature();
	double beta = 0;
	
	
}

/**
* Performs the velocity step
* 
* @param dt
*/
void fluid::velocityStep(double dt)
{
	double visc = 1;
	vector<double> temp [DIM];
	//Diffuse step
	diffuse(U1[0], U0[0], 0, visc, dt);
	diffuse(U1[1], U0[1], 0, visc, dt);
	diffuse(U1[2], U0[2], 0, visc, dt);
	std::copy(std::begin(U1), std::end(U1), std::begin(temp));
	std::copy(std::begin(U0), std::end(U0), std::begin(U1));
	std::copy(std::begin(temp), std::end(temp), std::begin(U0));

	//Should technecaly project ...

	transport(U1[0], U0[0], U1, dt);
	transport(U1[1], U0[1], U1, dt);
	transport(U1[2], U0[2], U1, dt);

	std::copy(std::begin(U1), std::end(U1), std::begin(temp));
	std::copy(std::begin(U0), std::end(U0), std::begin(U1));
	std::copy(std::begin(temp), std::end(temp), std::begin(U0));

	//Should project again ....

}

/**
* performs the scalar step
*
* @param dt
*/
void fluid::scalarStep(double dt)
{

}
void fluid::step()
{
}

/**
* Checks that a point in the triangle's tangent space is inside
* the triangles edges 
* 
* @param x The point
* @param t The triangle's index
* @return wether or not the point is in the triangle
*/
bool fluid::isInTriangle(Vector3d x, int t) {
	Vector3d a = V.row(F.row(t)[0]);
	Vector3d b = V.row(F.row(t)[1]);
	Vector3d c = V.row(F.row(t)[2]);
	Vector3d p = x;
	Vector3d ab = b - a;
	Vector3d ac = c - a;
	Vector3d cb = c - b;

	Vector3d pa = p - a;
	Vector3d pb = p - b;
	Vector3d pc = p - c;

	Vector3d abpa = ab.cross(pa);
	Vector3d acpc = ac.cross(pc);
	Vector3d cbpb = cb.cross(pb);

	if (abpa.dot(acpc) >= 0 && abpa.dot(cbpb) >= 0 && acpc.dot(cbpb)>=0) {
		return true;
	}
	return false;
}



