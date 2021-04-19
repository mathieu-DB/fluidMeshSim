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
#include <igl/cotmatrix.h>
#include <algorithm> 
#include <limits>
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
Fluid::Fluid(trimesh::trimesh_t mesh,
	Eigen::MatrixXd v,
	Eigen::MatrixXi f,
	Eigen::MatrixXd uv,
	Eigen::MatrixXd v_fused,
	Eigen::MatrixXi f_fused,
	map<int, int> v_mapsf,
	map<int, int> v_mapfs,
	map<int, int> dupe_map)
{
	V = v;
	F = f;
	V_fused = v_fused;
	F_fused = f_fused;
	V_mapFS = v_mapfs;
	V_mapSF = v_mapsf;
	Fluid::dupe_map = dupe_map;
	Fluid::uv = uv;
	Fluid::mesh = mesh;
	tools::CreateSparseWeightMatrix(V_fused, F_fused, L);
	tools::CreateSparseWeightMatrix(V, F, L_split);
	//igl::cotmatrix(V_fused, F_fused, L);
	//igl::cotmatrix(V, F, L_split);
	faces = F.rows();
	V_split_rows = V.rows();
	V_fused_rows = V_fused.rows();
	Area = tools::ComputeMeshSurfaceArea(V_fused, F_fused);

}

/**
* Setups the varius memory structures of the program
* (i.e. U0, U1, temperature0, temperature1)
*/
void Fluid::setup()
{
	elapsed = 0;
	U0[0].resize(V_fused.rows());
	U0[1].resize(V_fused.rows());
	U0[2].resize(V_fused.rows());
	std::fill(U0[0].begin(), U0[0].end(), 0);
	std::fill(U0[1].begin(), U0[1].end(), 0);
	std::fill(U0[2].begin(), U0[2].end(), 0);
	U1[0].resize(V_fused.rows());
	U1[1].resize(V_fused.rows());
	U1[2].resize(V_fused.rows());
	std::fill(U1[0].begin(), U1[0].end(), 0);
	std::fill(U1[1].begin(), U1[1].end(), 0);
	std::fill(U1[2].begin(), U1[2].end(), 0);

	temperature0.resize(V_fused.rows());
	temperature1.resize(V_fused.rows());
	std::fill(temperature0.begin(), temperature0.end(), 0);
	std::fill(temperature1.begin(), temperature1.end(), 0);

	div.resize(V_fused.rows());
	p.resize(V_fused.rows());



	//for (int i = 0; i < V_fused.rows(); i++) {
	//	div.emplace_back(0);
	//	p.emplace_back(0);
	//	U0[0].emplace_back(0);
	//	U0[1].emplace_back(0);
	//	U0[2].emplace_back(0);

	//	U1[0].emplace_back(0);
	//	U1[1].emplace_back(0);
	//	U1[2].emplace_back(0);

	//	temperature0.emplace_back(0);
	//	/*temperature0[1].emplace_back(0);
	//	temperature0[2].emplace_back(0);*/
	//	temperature1.emplace_back(0);
	//	//temperature1[1].emplace_back(0);
	//	//temperature1[2].emplace_back(0);

	//}
}

int Fluid::IX(int i, int j)
{
	return 0;
}

void Fluid::setBoundary(int b, vector<double> x)
{
	int i;
	//for (i = 1; i <= N; i++) {
	//	x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
	//	x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
	//	x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
	//	x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
	//}
	//x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
	//x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
	//x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
	//x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

/**
* Gets the velocity at the given point using interpolation
*
* @param x
* @param vel
*/
void Fluid::getVelocity(Vector3d x, Vector3d vel)
{
	//getVelocity(x, U0, vel);
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
	/*vel[0] = interpolate(x,U[0]);
	vel[1] = interpolate(x,U[1]);*/
}

/**
* Interpolates with the vertices values using barycentric coordinates
* to point x 
* @param x Point needing an interpolation  (Assuming 2D vector)
* @param t Triangle index
* @param s Scalar values to interpolate
* @return interpolated value
*/
double Fluid::interpolate(VectorXd x, int t, vector<double> s)
{
	VectorXd a = uv.row(F.row(t)[0]);
	VectorXd b = uv.row(F.row(t)[1]);
	VectorXd c = uv.row(F.row(t)[2]);

	VectorXd ca = c - a;
	VectorXd ba = b - a;
	VectorXd bc = b - c;
	VectorXd xa = x - a;
	VectorXd xb = x - b;
	VectorXd xc = x - c;

	double T = std::abs( (ca[0] * ba[1] - ba[0] * ca[1])) / 2;
	double alpha = std::abs((bc[0] * xc[1] - xc[0] * bc[1]) )/ 2;
	double beta = std::abs((ca[0] * xa[1] - xa[0] * ca[1]) )/ 2;
	double charli = std::abs(ba[0] * xb[1] - xb[0] * ba[1]) / 2;

	return (alpha * s[F_fused.row(t)[0]] + beta * s[F_fused.row(t)[1]] + charli * s[F_fused.row(t)[2]]) / T;
	
}

/**
* Performs a simple Forward Euler particle trace using the current velocity
* field.
*
* @param x0 Current particle location
* @param h  Time step
* @param x1 Final particle location
*/
void Fluid::traceParticle(Vector3d x0, double h, Vector3d& x1)
{
	//traceParticle(x0, U0, h, x1);
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
void Fluid::traceParticle(Vector3d x0, vector<double> U[], double h, Vector3d& x1)
{
	
	/*Eigen::Vector3d vec;
	vec.setZero();
	getVelocity(x0, U, vec);
	vec *= h;
	x1 = x0 + vec;*/

}

void Fluid::traceParticleFromVertex(int v, vector<double> U[], double h, Vector2d& x1) {
	/*Eigen::Vector2d vec;
	vec.setZero();
	vec[0] = U[0][v];
	vec[1] = U[1][v];
	vec *= h;
	x1 = uv.row(V_mapFS[v]) + vec;*/
}

/**
* Diffuse the given scalar field by the given amount. 
*
* @param S1   diffused quantities
* @param S0   initial quantities
* @param diff diffusion coefficient
* @param dt   time step
*/
void Fluid::diffuse(vector<double>& S1, vector<double>& S0, int b, double diff, double dt)
{
	int k;
	// (1/dx)^2 == N*N
	double diffRate = dt * diff*V_fused_rows;
	for (k = 0; k < iterations; k++) {


		for (int i = 0; i < V_fused_rows; i++) {
			double c = S0[i];
			int n = 0;
			double t = L.coeff(i,i);
			for (SparseMatrix<double>::InnerIterator it(L, i); it; ++it)
			{
				if (it.row() != i) {
					c += diffRate * it.value()/t * S1[it.row()];
					n++;
				}
				else {
					t = it.value();
				}
				
			}S1[i] = c / ((1 + n * diffRate));

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
void Fluid::transport(vector<double>& s1, vector<double>& s0, vector<double> U[], double dt)
{
	int i;
	double x, y;
	Vector2d p1;
	p1.setZero();
	Vector2d p2;

	for (i = 0; i < V_fused_rows; i++) {
		p2.setZero();

		p1 = uv.row(V_mapFS[i]);
		
		Eigen::Vector2d vec;
		vec.setZero();
		vec[0] = U[0][i];
		vec[1] = U[1][i];
		vec *= dt;

		int j = V_mapFS[i];
				
		p2 = p1 + vec;
		int t = -1;
		if (p2.isApprox(p1)) {
			if(scalarAdvect){
				temperature1[i] = temperature0[i];
			}
			if (velocityAdvect) {
				U1[0][i] = U0[0][i];
				U1[1][i] = U0[1][i];
			}
			//s1[i] = s0[i];
		}
		else {
			int v1 = indentifyTriangle(j, p2, t, -1);
			if (t < 0) {
				if (isDupe(v1)) {
					int v2 = dupe_map[v1];
					p2 = uv.row(v2);
					p2+=vec;
					int v = indentifyTriangle(v2, p2, t, -1);
					if (t < 0) {
						//vec = -vec;
						//p2 = uv.row(v2);
						//p2 += vec;
						//v = indentifyTriangle(v2, p2, t, -1);
						//if (t < 0) {
						//	// << "Triangle lookup failed on both vertice ....." << endl;
						//	//exit(1);
						//}
						
					}
				}
			}
			if (t >= 0) {
				if (scalarAdvect) {
					temperature1[i] = interpolate(p2, t, temperature0);
				}
				if (velocityAdvect) {
					U1[0][i] = interpolate(p2, t, U0[0]);
					U1[1][i] = interpolate(p2, t, U0[1]);
				}
				//s1[i] = interpolate(p2, t, s0);
			}
			else {
				if (scalarAdvect) {
					temperature1[i] = temperature0[i];
				}
				if (velocityAdvect) {
					U1[0][i] = U0[0][i];
					U1[1][i] = U0[1][i];
				}
				//s1[i] = s0[i];
			}
			
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
	Vector2d x(1, 0);
	Vector2d y(0, 1);
	std::fill(div.begin(), div.end(), 0);
	for (int i = 0; i < V_split_rows; i++) {
		for (SparseMatrix<double>::InnerIterator it(L_split, i); it; ++it)
		{
			if (it.row() != i) {
				VectorXd s = (uv.row(it.row()) - uv.row(i));
				double d = s.norm();
				s /= d;
				div[V_mapSF[i]] += d*(s.dot(x) * U0[0][V_mapSF[it.row()]] + s.dot(y) * U0[1][V_mapSF[it.row()]]);
			}
			
		}
		div[V_mapSF[i]] *= -0.5;
		p[V_mapSF[i]] = 0;
	}

	for (int k = 0; k < iterations; k++) {
		for (int i = 0; i < V_split_rows; i++) {
			p[V_mapSF[i]] = div[V_mapSF[i]];
			int n = 0;
			for (SparseMatrix<double>::InnerIterator it(L_split, i); it; ++it)
			{
				if (it.row() != i) {
					p[V_mapSF[i]] += p[V_mapSF[it.row()]];
					n++;
				}
			}
			p[V_mapSF[i]] /= n;
		}
	}

	for (int i = 0; i < V_split_rows; i++) {
		int n = 0;
		for (SparseMatrix<double>::InnerIterator it(L_split, i); it; ++it)
		{
			if (it.row() != i) {
				U0[0][V_mapSF[i]] -= 0.5* p[V_mapSF[it.row()]];
				U0[1][V_mapSF[i]] -= 0.5* p[V_mapSF[it.row()]];
				n++;
			}
			
		}
		U0[0][V_mapSF[i]] /= n;
		U0[1][V_mapSF[i]] /= n;
	}
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
void Fluid::addSource(vector<double>& S, double dt, Vector3d x, double amount, int v)
{
	//Use cotan matrix to distribute in a weighted fashion
	if (v >= 0) {
		/*double det = L.coeff(v, v);
		for (SparseMatrix<double>::InnerIterator it(L, v); it; ++it)
		{
			if (it.row() != v) {
				S[it.row()] += it.value() / det * amount * dt;
			}
		}*/
		S[v] += amount * dt;

	}
	
	else {
		/*double ir = ((x.x() * N) + 0.5);
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
		S[IX(i + 1, j + 1)] += alpha3 * amount * dt;*/
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
	for (int i = 0; i < V_fused_rows; i++) {
		referenceTemperature += temperature0[i];
	}
	referenceTemperature /= V_fused_rows;
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
	for (int i = 0; i < V_fused_rows; i++) {
			int lowest = i;
			Eigen::VectorXd p = V_fused.row(i);
			double lowestZ = p[2];
			for (SparseMatrix<double>::InnerIterator it(L, i); it; ++it)
			{
				
				Eigen::VectorXd p2= V_fused.row(it.row());
				if(p2[2] > lowestZ){
					lowest = it.row();
					lowestZ = p2[2];
				}
			}
			if (lowest != i) {
				int j = V_mapFS[i];
				int uv_low = V_mapFS[lowest];
				//if currently analysed vertex is a dupe from the seam
				//find to which dupe the "lowest" vertex is neighbour in uv
				if (isDupe(j)) {
					if (isDupe(uv_low) ){
						int low2 = dupe_map[uv_low];
						double d1 = (uv.row(j) - uv.row(uv_low)).squaredNorm();
						double d2 = (uv.row(j) - uv.row(low2)).squaredNorm();
						uv_low = (d1 < d2) ? uv_low : low2;
					}
					else {
						int j2 = dupe_map[j];
						double d1 = (uv.row(j) - uv.row(uv_low)).squaredNorm();
						double d2 = (uv.row(j2) - uv.row(uv_low)).squaredNorm();
						lowest = uv_low;
						j = (d1 > d2) ? j2 : j;
					}
					
					
				}
				Eigen::VectorXd vec =   uv.row(uv_low)-uv.row(j);
				vec.normalize();
				vec *= beta * dt * (referenceTemperature - 0.5f * (temperature0[i] + temperature0[lowest]));
				U[0][i] += vec[0];
				U[1][i] += vec[1];
			}
	}

}

/**
* Performs the velocity step
* 
* @param dt
*/
void Fluid::velocityStep(double dt)
{
	double visc = viscosity;
	vector<double> temp[3];
	if (velocityDiffuse) {
		diffuse(U1[0], U0[0], 1, visc, dt);
		diffuse(U1[1], U0[1], 2, visc, dt);
		U1[0].swap(U0[0]);
		U1[1].swap(U0[1]);
	/*	std::copy(std::begin(U1), std::end(U1), std::begin(temp));
		std::copy(std::begin(U0), std::end(U0), std::begin(U1));
		std::copy(std::begin(temp), std::end(temp), std::begin(U0));*/
	}
	if (velocityProject) {
		project(U0);
	}
	if (velocityAdvect) {
		transport(U1[0], U0[0], U1, dt);
		//transport(U1[1], U0[1], U1, dt);
		U1[0].swap(U0[0]);
		U1[1].swap(U0[1]);
	/*	std::copy(std::begin(U1), std::end(U1), std::begin(temp));
		std::copy(std::begin(U0), std::end(U0), std::begin(U1));
		std::copy(std::begin(temp), std::end(temp), std::begin(U0));*/
		if (scalarAdvect) {
			temperature1.swap(temperature0);
		}

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
		temperature1.swap(temperature0);
		/*temp = temperature1;
		temperature1 = temperature0;
		temperature0 = temp;*/
	}

	if (scalarAdvect) {
		//transport(temperature1, temperature0, U0, dt);
		//temperature1.swap(temperature0);
		//temp = temperature1;
		//temperature1 = temperature0;
		//temperature0 = temp;
		//setBoundary(0, temperature0);
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
bool Fluid::isInTriangle(VectorXd x, int t) {
	
	VectorXd a = uv.row(F.row(t)[0]);
	VectorXd b = uv.row(F.row(t)[1]);
	VectorXd c = uv.row(F.row(t)[2]);

	VectorXd ca = c - a;
	VectorXd ab = a - b;
	VectorXd bc = b - c;
	VectorXd xa = x - a;
	VectorXd xb = x - b;
	VectorXd xc = x - c;
	
	double alpha = (bc[0] * xc[1]) - (xc[0] * bc[1]);
	double beta = (ca[0] * xa[1]) - (xa[0] * ca[1]);
	double charli = (ab[0] * xb[1]) - (xb[0] * ab[1]);

	if (alpha * beta >= 0 && beta * charli >= 0 && alpha * charli >= 0) return true;

	return false;
}

void Fluid::createSource(int v, double amount)
{	
	
	double d2 = 0;
	double t;
	
	Source s(v, amount);
	sources.push_back(s);
}


double Fluid::interpolateTempForVertex(Vector2d x)
{
	/*double ir = ((x[0] * 0.05) + 0.5);
	double jr = ((x[1] * 0.05) + 0.5);
	ir = std::min(std::max(ir, 0.5), N + 0.5);
	jr = std::min(std::max(jr, 0.5), N + 0.5);
	int i = (int)ir;
	int j = (int)jr;

	double a1 = ir - i;
	double a0 = 1 - a1;

	double b1 = jr - j;
	double b0 = 1 - b1;

	return a0 * (b0 * temperature0[IX(i, j)] + b1 * temperature0[IX(i, j + 1)]) + a1 * (b0 * temperature0[IX(i + 1, j)] + b1 * temperature0[IX(i + 1, j + 1)]);*/
	return 0;
}

double Fluid::getTempAtVertex(int v) {
	return temperature0[V_mapSF[v]];
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

bool Fluid::isDupe(int v) {
	return (dupe_map.find(v) == dupe_map.end()) ? false : true;
}

/*
* Finds in which triangle the point x lies by searching from point startV  (all in uv coord)
* 
* @param startV Index in the V_split matrix of the starting vertex
* @param x		Point with unknown triangle
* @param t		Returning value for the triangle index in F_split
* @return		Returns -1 if triangle is found otherwise returns the index of the boundary vertex
				closest to the point (the point is outside the boundary and needs to be teleported to the boundary's dupe)
*/
int Fluid::indentifyTriangle(int startV, VectorXd x, int& t, int pred) {

	
	
	//Check all neighbors (using cot matrix)
	Vector2d pointVector = x;
	pointVector-=uv.row(startV);
	double closeDist = std::numeric_limits<double>::max();
	int closest = startV;
	int n = 0;
	//Find closest edge
	for (SparseMatrix<double>::InnerIterator it(L_split, startV); it; ++it)
	{
		//For neighbour vertex, build edge and check if x is in that direction
		if (startV != it.row()) {
			Vector2d edge = uv.row(it.row());
			edge-=uv.row(startV);
			double dot = pointVector.dot(edge);
			//ignore vertices not in the same general direction
			if (dot >= 0) {
				//using vector rejection to get the distance of the the point to the edge
				double d = (pointVector - (dot / edge.dot(edge)) * edge).squaredNorm();
				if (d < closeDist) {
					closest = it.row();
					closeDist = d;
				}			
			}
			n++;
		}
		
	}
	if (closest == pred) {
		t = -1;
		return closest;
	}
	if (closest != startV) {
		index_t hei = mesh.find_halfedge_of_neib_vertex(startV, closest, n+1);
		//Check both triangles around the edge to see if point is in them
		int f = mesh.get_he_face(hei);
		/*if (f == 23 && startV == 14) {
			cout << x << endl;
			cout << uv.row(startV) << endl;
			cout << uv.row(closest) << endl;
		}*/
		int f1 = f;
		if (f != -1) { //Triangle 1 exists
			if (isInTriangle(x, f)) {
				t = f;
				return -1;
			}
			else {
				f = mesh.get_he_face(mesh.get_opposite(hei));
				if (f != -1) {
					if (isInTriangle(x, f)) {
						t = f;
						return -1;
					}
					else { //Both triangles exist but point is further ahead
						//Call recursivly on the following vertex
						return indentifyTriangle(closest, x, t, startV);
					}
				}
				else { //triangle 2 doesn't exist
					double d_sv = (uv.row(closest) - uv.row(startV)).squaredNorm();
					double d_sx = pointVector.squaredNorm();

					//If point x is further from start than the following vertex is recurse
					if (d_sv < d_sx) {
						return indentifyTriangle(closest, x, t, startV);
					}
					else { //Check on which side of the edge the point lies.
						VectorXi vertices = F.row(f1);
						//Get 3rd vertex from triangle
						int v3 = (vertices[0] != startV && vertices[0] != closest) ? vertices[0] : (vertices[1] != startV && vertices[1] != closest) ? vertices[1] : vertices[2];
						
						//Check if x is on the same side of the edge as v3
						Vector2d cs = uv.row(closest) - uv.row(startV);
						Vector2d v3s = uv.row(v3) - uv.row(startV);
						double closest_v3 = v3s[0] * pointVector[1] - pointVector[0] * v3s[1];
						double closest_x = cs[0] * pointVector[1] - pointVector[0] * cs[1];

						if (closest_v3 * closest_x < 0) { //Same side
							//Recurse
							return indentifyTriangle(closest, x, t, startV);
						}
						else {
							// x is outside the boundary, caller needs to recompute from dupe vector
							t = -1;
							return closest;
						}
					}

				}
			}
		}
		else { //Triangle 1 doesn't exist
			f = mesh.get_he_face(mesh.get_opposite(hei));
			if (f != -1) {
				if (isInTriangle(x, f)) {
					t = f;
					return -1;
				}
				else { //Not in triangle 2
					double d_sv = (uv.row(closest) - uv.row(startV)).squaredNorm();
					double d_sx = pointVector.squaredNorm();

					//If point x is further from start than the following vertex is recurse
					if (d_sv < d_sx) {	
						return indentifyTriangle(closest, x, t, startV);
					}
					else { //Check on which side of the edge the point lies.
						VectorXi vertices = F.row(f);
						//Get 3rd vertex from triangle
						int v3 = (vertices[0] != startV && vertices[0] != closest) ? vertices[0] : (vertices[1] != startV && vertices[1] != closest) ? vertices[1] : vertices[2];

						//Check if x is on the same side of the edge as v3
						Vector2d cs = uv.row(closest) - uv.row(startV);
						Vector2d v3s = uv.row(v3) - uv.row(startV);
						double closest_v3 = v3s[0] * pointVector[1] - pointVector[0] * v3s[1];
						double closest_x = cs[0] * pointVector[1] - pointVector[0] * cs[1];

						if (closest_v3 * closest_x < 0) { //Same side
							//Recurse
							return indentifyTriangle(closest, x, t, startV);
						}
						else {
							// x is outside the boundary, caller needs to recompute from dupe vector
							t = -1;
							return closest;
						}
					}

				}
			}
			else {
				cout << "Error in halfedge structure: Edge not associated with any triangles !!!" << endl;
				cout << "See Fluid::indentifyTriangle in fluid.cpp" << endl;
			}
		}
	}

	
	return -2;
}


