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


//TODO Reproduce 2D grid behavior from assignement 3
//TODO Apply grid simulation to parametrized mesh ignoring 3rd component
//TODO Add 3rd dimension (possibly only applies to temperature and boyency)



/**
* Constructor
* 
* @param mesh Halfedge mesh structure
* @param v vertices
* @param f faces
*/
Fluid::Fluid(trimesh::trimesh_t mesh, Eigen::MatrixXd v, Eigen::MatrixXi f, Eigen::MatrixXd uv)
{
	V = v;
	F = f;
	Fluid::uv = uv;
	mesh = mesh;
	igl::cotmatrix(V, F, L);
	faces = F.rows();

}

/**
* Setups the varius memory structures of the program
* (i.e. U0, U1, temperature0, temperature1)
*/
void Fluid::setup()
{
	elapsed = 0;
	dx = 1.0 / N;
	int np2s = (N + 2) * (N + 2);
	for (int i = 0; i < V.rows(); i++) {
		div.emplace_back(0);
		p.emplace_back(0);
		U0[0].emplace_back(0);
		U0[1].emplace_back(0);
		U0[2].emplace_back(0);

		U1[0].emplace_back(0);
		U1[1].emplace_back(0);
		U1[2].emplace_back(0);

		temperature0.emplace_back(0);
		/*temperature0[1].emplace_back(0);
		temperature0[2].emplace_back(0);*/
		temperature1.emplace_back(0);
		//temperature1[1].emplace_back(0);
		//temperature1[2].emplace_back(0);

	}
}

int Fluid::IX(int i, int j)
{
	return i*(N+2)+j;
}

void Fluid::setBoundary(int b, vector<double> x)
{
	int i;
	for (i = 1; i <= N; i++) {
		x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
		x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
		x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
	}
	x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
	x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
	x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

/**
* Gets the velocity at the given point using interpolation
*
* @param x
* @param vel
*/
void Fluid::getVelocity(Vector3d x, Vector3d vel)
{
	getVelocity(x, U0, vel);
}

/**
* Gets the velocity in the provided velocity field at the given point using
* interpolation
*
* @param x
* @param U
* @param vel
*/
void Fluid::getVelocity(Vector3d x, vector<double> U[], Vector3d vel)
{
	vel[0] = interpolate(x,U[0]);
	vel[1] = interpolate(x,U[1]);
}

/**
* Interpolates with the vertices values using barycentric coordinates
* to point x 
* @param x Point needing an interpolation
* @param t Triangle index
* @param s Scalar values to interpolate
* @return interpolated value
*/
double Fluid::interpolate(Vector3d x, vector<double> s)
{
	double ir = ((x[0] * N) + 0.5);
	double jr = ((x[1] * N) + 0.5);
	ir = std::min(std::max(ir, 0.5), N + 0.5);
	jr = std::min(std::max(jr, 0.5), N + 0.5);
	int i = (int)ir;
	int j = (int)jr;

	double a1 = ir - i;
	double a0 = 1 - a1;

	double b1 = jr - j;
	double b0 = 1 - b1;

	return a0 * (b0 * s[IX(i, j)] + b1 * s[IX(i, j + 1)]) + a1 * (b0 * s[IX(i + 1, j)] + b1 * s[IX(i + 1, j + 1)]);
}

/**
* Performs a simple Forward Euler particle trace using the current velocity
* field.
*
* @param x0 Current particle location
* @param h  Time step
* @param x1 Final particle location
*/
void Fluid::traceParticle(Vector3d x0, double h, Vector3d &x1)
{
	traceParticle(x0, U0, h, x1);
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
void Fluid::traceParticle(Vector3d x0, vector<double> U[], double h, Vector3d &x1)
{
	
	Eigen::Vector3d vec;
	vec.setZero();
	getVelocity(x0, U, vec);
	vec *= h;
	x1 = x0 + vec;

}

/**
* Diffuse the given scalar field by the given amount. 
*
* @param S1   diffused quantities
* @param S0   initial quantities
* @param diff diffusion coefficient
* @param dt   time step
*/
void Fluid::diffuse(vector<double> &S1, vector<double> &S0, int b, double diff, double dt)
{
	int i, j, k;
	// (1/dx)^2 == N*N
	double diffRate = dt * diff * N * N;
	for (k = 0; k < iterations; k++) {

		for (int i = 0; i < V.rows(); i++) {
			double c = S0[i];
			int n = 0;
			for (SparseMatrix<double>::InnerIterator it(L, i); it; ++it)
			{
				if (it.row() != i) {
					c += diffRate * S1[it.row()];
					n++;
				}
				S1[i] += c / (1 + n * diffRate);
			}

		}


		//for (i = 1; i <= N; i++) {
		//	for (j = 1; j <= N; j++) {
		//		S1[IX(i, j)] = (S0[IX(i, j)]
		//			+ diffRate * (S1[IX(i - 1, j)] + S1[IX(i + 1, j)] + S1[IX(i, j - 1)] + S1[IX(i, j + 1)]))
		//			/ (1 + 4 * diffRate);
		//	}
		//}
		//setBoundary(b, S1);

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
void Fluid::transport(vector<double> &s1, vector<double> &s0, vector<double> U[], double dt)
{
	int i, j;
	double x, y;
	Vector3d p1;
	p1.setZero();
	Vector3d p2;
	
	for (i = 0; i < V.rows(); i++) {
		p2.setZero();
		p1 = uv.row(i);

		traceParticle(p1, U, dt, p2);
		//TODO find wich triangle the new point belongs to:
		//1. Find edge attached to original vertex that is closest to new point
		//2. Verify if point is in one of the 2 triangles attached to said edge
		//	if it is -> done
		//	else ->  REPEAT 1.  with the edge's 2nd vertex as a source
		// Repeat until found.


	}
	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			p2.setZero();

			x = i * dx - dx * 0.5f;
			y = j * dx - dx * 0.5f;
			p1[0] = x;
			p1[1] = y;
			p1[2] = 0;

			traceParticle(p1, U, dt, p2);

			s1[IX(i, j)] = interpolate(p2, s0);


		}
	}
}

/**
* Does the Poisson solve to make sure that velocities U respect incompressible
* flow
*
* @param U
*/
void Fluid::project(vector<double> U[])
{
	
	for (int i = 0; i < V.rows(); i++) {
		/*for (SparseMatrix<double>::InnerIterator it(L, i); it; ++it)
		{
			if (it.row() != i) {
				div[i] +
			}
		}*/
		p[i] = 0;
	}
	for (int k = 0; k < iterations; k++) {
		for (int i = 0; i < V.rows(); i++) {
			p[i] = L.coeff(i, i);
			int n = 0;
			for (SparseMatrix<double>::InnerIterator it(L, i); it; ++it)
			{
				if (it.row() != i) {
					p[i] += p[it.row()];
					n++;
				}
			}
			p[i] /= n;
		}
	}

	for (int i = 0; i < V.rows(); i++) {
		int n = 0;
		for (SparseMatrix<double>::InnerIterator it(L, i); it; ++it)
		{
			if (it.row() != i) {
				U[0][i] += p[it.row()];
				U[1][i] += p[it.row()];
				n++;
			}
			U[0][i] /= n;
			U[1][i] /= n;
		}
	}
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			div[IX(i, j)] = -0.5f * dx
				* (U[0][IX(i + 1, j)] - U[0][IX(i - 1, j)] + U[1][IX(i, j + 1)] - U[1][IX(i, j - 1)]);
			p[IX(i, j)] = 0;
		}
	}
	setBoundary(0, div);
	setBoundary(0, p);
	for (int k = 0; k < iterations; k++) {
		for (int i = 1; i <= N; i++) {
			for (int j = 1; j <= N; j++) {
				p[IX(i, j)] = (div[IX(i, j)] + p[IX(i + 1, j)] + p[IX(i - 1, j)] + p[IX(i, j + 1)]
					+ p[IX(i, j - 1)]) / 4.0;
			}
		}
		setBoundary(0, p);
	}

	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			U[0][IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
			U[1][IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
		}
	}
	setBoundary(1, U[0]);
	setBoundary(2, U[1]);
}

/**
* Adds a force at a given point in the provided velocity field.
*
* @param U
* @param dt
* @param x
* @param f
*/
void Fluid::addForce(vector<double> U[], double dt, Vector3d x, Vector3d f, int v)
{
	addSource(U[0], dt, x, f.x(), v);
	addSource(U[1], dt, x, f.y(), v);
}

/**
* Adds some time step scaled amount to the provided scalar field. 
*
* @param S      quantity field to modify
* @param dt     time step
* @param x      position
* @param amount amount
*/
void Fluid::addSource(vector<double> &S, double dt, Vector3d x, double amount, int v)
{
	//Use cotan matrix to distribute in a weighted fashion
	if (v >= 0) {
		double det =std::abs( L.coeff(v, v));
		for (SparseMatrix<double>::InnerIterator it(L, v); it; ++it)
		{
			if (it.row() != v) {
				S[it.row()] += it.value() / det * amount * dt;
			}
		}
	}
	else {
		double ir = ((x.x() * N) + 0.5);
		double jr = ((x.y() * N) + 0.5);
		int i = (int)ir;
		int j = (int)jr;

		double alpha1 = (1 - (ir - i)) * (1 - (jr - j));
		double alpha2 = (ir - i) * (1 - (jr - j));
		double alpha3 = (ir - i) * (jr - j);
		double alpha4 = (1 - (ir - i)) * (jr - j);
		S[IX(i, j)] += alpha1 * amount * dt;
		S[IX(i, j + 1)] += alpha4 * amount * dt;
		S[IX(i + 1, j)] += alpha2 * amount * dt;
		S[IX(i + 1, j + 1)] += alpha3 * amount * dt;
	}
	

}

/**
* Gets the average temperature of the continuum. (use this in computing
* buoyancy forces)
*
* @return average temperature
*/
double Fluid::getReferenceTemperature()
{
	int count = 0;
	double referenceTemperature = 0;
	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= N; j++) {
			referenceTemperature += temperature0[IX(i, j)];
			count++;
		}
	}
	referenceTemperature /= count;
	return referenceTemperature;
}

/**
* Applies buoyancy forces to the velocity field due to temperature. Use Foster
* and Metaxis [1997] Equation 2: F_{bv} = beta g_v (T_0 - T_k)
*
* @param U
* @param dt
*/
void Fluid::addTemperatureForce(vector<double> U[], double dt)
{

	double referenceTemperature = getReferenceTemperature();
	double beta = bouyancy;
	
	//Using the z component of the original vertices to add temperature force to the 
	//velocity vector. We find the lowest vertex attached to the current one (if lower
	//than it) and we have the force go that way
	for (int i = 0; i < V.rows(); i++) {
			int lowest = i;
			Eigen::VectorXd p = V.row(i);
			double lowestZ = p[2];
			for (SparseMatrix<double>::InnerIterator it(L, i); it; ++it)
			{
				
				Eigen::VectorXd p2= V.row(it.row());
				if(p2[2] < lowestZ){
					lowest = it.row();
					lowestZ = p2[2];
				}
			}
			if (lowest != i) {
				Eigen::VectorXd vec = uv.row(i) - uv.row(lowest);
				vec.normalize();
				vec *= beta * dt * (referenceTemperature - 0.5f * (temperature0[i] + temperature0[lowest]));
				U[0][i] += vec[0];
				U[1][i] += vec[1];
			}
	}

	// TODO: Objective 7: change velocities based on the temperature. Don't forget
	// to set Boundaries after modifying velocities!
	//for (int i = 1; i <= N; i++) {
	//	for (int j = 1; j <= N; j++) {
	//		//buoyancy times step size times temperature delta
	//		U[1][IX(i, j)] += beta * dt * (referenceTemperature - 0.5f * (temperature0[IX(i, j + 1)] + temperature0[IX(i, j)]));
	//	}
	//}
	//setBoundary(2, U[1]);
}

/**
* Performs the velocity step
* 
* @param dt
*/
void Fluid::velocityStep(double dt)
{
	double visc = viscosity;
	vector<double> temp[DIM];
	if (velocityDiffuse) {
		diffuse(U1[0], U0[0], 1, visc, dt);
		diffuse(U1[1], U0[1], 2, visc, dt);
		std::copy(std::begin(U1), std::end(U1), std::begin(temp));
		std::copy(std::begin(U0), std::end(U0), std::begin(U1));
		std::copy(std::begin(temp), std::end(temp), std::begin(U0));
	}
	if (velocityProject) {
		project(U0);
	}
	if (velocityAdvect) {
		transport(U1[0], U0[0], U1, dt);
		transport(U1[1], U0[1], U1, dt);
		std::copy(std::begin(U1), std::end(U1), std::begin(temp));
		std::copy(std::begin(U0), std::end(U0), std::begin(U1));
		std::copy(std::begin(temp), std::end(temp), std::begin(U0));
		setBoundary(1, U0[0]);
		setBoundary(2, U0[1]);
	}
	if (velocityProject) {
		project(U0);
	}

}

/**
* performs the scalar step
*
* @param dt
*/
void Fluid::scalarStep(double dt)
{
	vector<double> temp;
	if (scalarDiffuse) {
		double diff = diffusion;
		diffuse(temperature1, temperature0, 0, diff, dt);
		temp = temperature1;
		temperature1 = temperature0;
		temperature0 = temp;
	}

	if (scalarAdvect) {
		transport(temperature1, temperature0, U0, dt);
		temp = temperature1;
		temperature1 = temperature0;
		temperature0 = temp;
		setBoundary(0, temperature0);
	}

}
void Fluid::step()
{
	double dt = timeStep;
	for (Source s : sources) {
		addSource(temperature0, dt, s.location, s.amount, s.vertex);
	}
	addTemperatureForce(U0, dt);
	velocityStep(dt);
	scalarStep(dt);
	elapsed += dt;
}

/**
* Checks that a point in the triangle's tangent space is inside
* the triangles edges 
* 
* @param x The point
* @param t The triangle's index
* @return wether or not the point is in the triangle
*/
bool Fluid::isInTriangle(Vector3d x, int t) {
	return false;
}

void Fluid::createSource(int v, double amount)
{	
	int v = 0;
	double d2 = 0;
	double t;
	
	Source s(v, amount);
	s.vertex = v;
	sources.push_back(s);
}

double Fluid::interpolateTempForVertex(Vector2d x)
{
	double ir = ((x[0] * 0.05) + 0.5);
	double jr = ((x[1] * 0.05) + 0.5);
	ir = std::min(std::max(ir, 0.5), N + 0.5);
	jr = std::min(std::max(jr, 0.5), N + 0.5);
	int i = (int)ir;
	int j = (int)jr;

	double a1 = ir - i;
	double a0 = 1 - a1;

	double b1 = jr - j;
	double b0 = 1 - b1;

	return a0 * (b0 * temperature0[IX(i, j)] + b1 * temperature0[IX(i, j + 1)]) + a1 * (b0 * temperature0[IX(i + 1, j)] + b1 * temperature0[IX(i + 1, j + 1)]);
}

double Fluid::getMaxTemp() {
	double m = -DBL_MAX;
	for (double t : temperature0) {
		if (t > m) m = t;
	}
	return m;
}

double Fluid::getMinTemp() {
	double m = DBL_MAX;
	for (double t : temperature0) {
		if (t < m) m = t;
	}
	return m;
}


