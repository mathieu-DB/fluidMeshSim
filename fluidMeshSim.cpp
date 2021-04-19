#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <igl/read_triangle_mesh.h>
#include <igl/triangle/scaf.h>
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/readOBJ.h>
#include <igl/Timer.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/MappingEnergyType.h>
#include <igl/doublearea.h>
#include <igl/PI.h>
#include <igl/flipped_triangles.h>
#include <igl/topological_hole_fill.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <map>


#include "fluid.h"


Eigen::MatrixXd V_uv, uv_bnd, uv_edges1, uv_edges2;
Eigen::VectorXi bnd;
igl::Timer timer;
igl::triangle::SCAFData scaf_data;
Eigen::MatrixXd V, V_split , V_fused;
Eigen::MatrixXi F, F_split , F_fused;
std::map<int, int> V_mapSF{};
std::map<int, int> V_mapFS{};
std::map<int, int> dupe_map{};

bool show_uv = false;
float uv_scale = 0.2f;
bool open = true;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
    if (key == '1')
        show_uv = false;
    else if (key == '2')
        show_uv = true;

    if (key == ' ')
    {
        timer.start();
        igl::triangle::scaf_solve(scaf_data, 1);
        
        std::cout << "time = " << timer.getElapsedTime() << std::endl;
    }

    const auto& V_uv = uv_scale * scaf_data.w_uv.topRows(V.rows());
    if (show_uv)
    {
        
        viewer.data().clear();
        viewer.data().clear_points();
        for (int i = 0; i < uv_bnd.rows(); i++) {
            int j = bnd[i];

            uv_bnd.row(i) = Eigen::Vector3d(V_uv.row(j)[0], V_uv.row(j)[1],0);
            if (i < uv_bnd.rows() - 1) {
                uv_edges1.row(i) = uv_bnd.row(i);
            }
            if (i > 0) {
                uv_edges2.row(i - 1) = uv_bnd.row(i);
            }
        }
        viewer.data().set_mesh(V_uv, F);
        viewer.data().set_uv(V_uv);
        viewer.data().add_points(uv_bnd, Eigen::RowVector3d(1, 0, 0));

        viewer.core().align_camera_center(V_uv, F);

    }
    else
    {
        viewer.data().clear();
        for (int i = 0; i < uv_bnd.rows(); i++) {
            int j = bnd[i];

            uv_bnd.row(i) = Eigen::Vector3d(V.row(j)[0], V.row(j)[1], V.row(j)[2]);
            if (i < uv_bnd.rows() - 1) {
                uv_edges1.row(i) = uv_bnd.row(i);
            }
            if (i > 0) {
                uv_edges2.row(i - 1) = uv_bnd.row(i);
            }
        }
        viewer.data().set_mesh(V, F);
        viewer.data().set_uv(V_uv);
        viewer.core().align_camera_center(V, F);
        viewer.data().add_points(uv_bnd, Eigen::RowVector3d(1, 0, 0));
        viewer.data().add_edges(uv_edges1, uv_edges2, Eigen::RowVector3d(1, 0, 0));
        viewer.data().point_size = 5;
    }

    
    viewer.data().compute_normals();

    return false;
}

int main(int argc, char* argv[])
{
    using namespace std;
    // Load a mesh in OFF format
    igl::readOBJ("C:/Users/m_dav/Documents/UdeM/COMP559/FreshStart/libigl/tutorial/data/camel_b.obj", V, F);

    Eigen::MatrixXd bnd_uv, uv_init;

    Eigen::VectorXd M;
    igl::doublearea(V, F, M);
    std::vector<std::vector<int>> all_bnds;
    igl::boundary_loop(F, all_bnds);

    // Heuristic primary boundary choice: longest
    auto primary_bnd = std::max_element(all_bnds.begin(), all_bnds.end(), [](const std::vector<int>& a, const std::vector<int>& b) { return a.size() < b.size(); });

    bnd = Eigen::Map<Eigen::VectorXi>(primary_bnd->data(), primary_bnd->size());

    
    //Original mesh is build with a seam and a mesh flattenning in mind
    //vertices along the seam are pre duplicated to allow a cut
    //Need to create a fused mesh to use cotmatrix and to walk along the vertices

    tools::CreateFusedMesh(V, V_fused, F, F_fused, V_mapSF, V_mapFS, dupe_map, bnd);


    igl::map_vertices_to_circle(V, bnd, bnd_uv);

    trimesh::trimesh_t mesh;
    trimesh::trimesh_t uv_mesh;
    std::vector< trimesh::triangle_t > triangles;

    int kNumVertices = V.rows();
    int kNumFaces = F.rows();
    triangles.resize(kNumFaces);
    for (int i = 0; i < kNumFaces; ++i) {
        triangles[i].v[0] = F(i, 0);
        triangles[i].v[1] = F(i, 1);
        triangles[i].v[2] = F(i, 2);
    }

    std::vector< trimesh::edge_t > edges;
    trimesh::unordered_edges_from_triangles(triangles.size(), &triangles[0], edges);
    uv_mesh.build(kNumVertices, triangles.size(), &triangles[0], edges.size(), &edges[0]);


    bnd_uv *= sqrt(M.sum() / (2 * igl::PI));
    if (all_bnds.size() == 1)
    {
        if (bnd.rows() == V.rows()) // case: all vertex on boundary
        {
            uv_init.resize(V.rows(), 2);
            for (int i = 0; i < bnd.rows(); i++)
                uv_init.row(bnd(i)) = bnd_uv.row(i);
        }
        else
        {
            igl::harmonic(V, F, bnd, bnd_uv, 1, uv_init);
            if (igl::flipped_triangles(uv_init, F).size() != 0)
                igl::harmonic(F, bnd, bnd_uv, 1, uv_init); // fallback uniform laplacian
        }
    }
    else
    {
        // if there is a hole, fill it and erase additional vertices.
        all_bnds.erase(primary_bnd);
        Eigen::MatrixXi F_filled;
        igl::topological_hole_fill(F, bnd, all_bnds, F_filled);
        igl::harmonic(F_filled, bnd, bnd_uv, 1, uv_init);
        uv_init.conservativeResize(V.rows(), 2);
    }

    Eigen::VectorXi b; Eigen::MatrixXd bc;
    igl::triangle::scaf_precompute(V, F, uv_init, scaf_data, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 0);

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;




    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    
    viewer.data().set_mesh(V, F);
    V_uv = uv_scale * scaf_data.w_uv.topRows(V.rows());
    uv_bnd.resize(bnd_uv.rows(), 3);
    uv_edges1.resize(bnd_uv.rows() - 1, 3);
    uv_edges2.resize(bnd_uv.rows() - 1, 3);
    for (int i = 0; i < uv_bnd.rows(); i++) {
        int j = bnd[i];

        uv_bnd.row(i) = Eigen::Vector3d(V.row(j)[0], V.row(j)[1], V.row(j)[2]);
        if (i < uv_bnd.rows() - 1) {
            uv_edges1.row(i) = uv_bnd.row(i);
        }
        if (i > 0) {
            uv_edges2.row(i-1) = uv_bnd.row(i);
        }
    }
    viewer.core().is_animating = true;
    igl::triangle::scaf_solve(scaf_data, 7);
    const auto& V_uv = uv_scale * scaf_data.w_uv.topRows(V.rows());


    
    Fluid fluid(uv_mesh, V, F, V_uv, V_fused, F_fused, V_mapSF, V_mapFS, dupe_map);
    fluid.setup();

    int N = 16;
    Eigen::MatrixX3d v_temps;
    Eigen::MatrixX3d ones;
    v_temps.resize(V.rows(), Eigen::NoChange);
    ones.resize(V.rows(), Eigen::NoChange);
    v_temps.setZero();
    ones.setOnes();


    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        //menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Controls", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // Expose variable directly ...
            ImGui::Checkbox("Velocity Diffuse", &fluid.velocityDiffuse);
            ImGui::Checkbox("Velocity Project", &fluid.velocityProject);
            ImGui::Checkbox("Velocity Advect", &fluid.velocityAdvect);
            ImGui::Checkbox("Scalar Diffuse", &fluid.scalarDiffuse);
            ImGui::Checkbox("Scalar Advect", &fluid.scalarAdvect);


            ImGui::SliderFloat("Viscosity", &fluid.viscosity, 1e-8, 1);
            ImGui::SliderFloat("Diffusion", &fluid.diffusion, 1e-8, 1);
            ImGui::SliderFloat("Bouyancy", &fluid.bouyancy, -1, 1);
            ImGui::SliderInt("Iterations", &fluid.iterations, 0, 1000);
            ImGui::SliderFloat("Time step", &fluid.timeStep, 0.001, 1);
            ImGui::SliderFloat("Gravity", &fluid.gravity, -10, 15);


        }
    };

   // viewer.data().set_colors(v_temps);
    viewer.callback_key_down = &key_down;

    // Enable wireframe
    viewer.data().show_lines = true;

    // Draw checkerboard texture
    viewer.data().show_texture = true;
    Eigen::Vector3d s(0.5, 0.5, 0);
    fluid.createSource(1236 , 20);
    fluid.createSource(50, -20);

    std::cerr << "Press space for running an iteration." << std::endl;
    std::cerr << "Press 1 for Mesh 2 for UV" << std::endl;
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool
    {
        // Create orbiting animation
        fluid.step();
        double max = std::max(fluid.getMaxTemp(),1.0);
        double min = std::abs(std::min(fluid.getMinTemp(),-1.0));
        for (int i = 0; i < V.rows(); i++) {
            int v = V_mapSF[i];
            double t = fluid.getTempAtVertex(v);
            s = Eigen::Vector3d(t > 0 ? t : 0, 0.125, t < 0 ? -t : 0);
            /*if (t == 0) {
                s[0] = max;
                s[1] = max;
                s[2] = max;
            }
            else if(t>0){
                s[0] = max;
                s[1] =  max-t;
                s[2] = max-t;
            }
            else {
                s[0] = min - std::abs(t);
                s[1] = min - std::abs(t);
                s[2] = min;
            }*/
            
            v_temps.row(i) = s;
        }
        viewer.data().set_colors(v_temps);
        return false;
    };
    cout << "1 : " << V_uv.row(1) << endl;
    cout << "2 : " << V_uv.row(2) << endl;
    cout << "521 : " << V_uv.row(521) << endl;
    cout << "1509 : " << V_uv.row(1509) << endl;
    // Launch the viewer
    viewer.launch();
    


}