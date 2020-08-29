#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "geometrycentral/surface/direction_fields.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/utilities/utilities.h"

#include "animated_quantity.h"
#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "flow.h"
#include "streamlines.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

VertexData<Vector2> velocityField;

// Some algorithm parameters
float dt  = 1. / 24.;
int steps = 8;

void drawEigenvectorInterpolation() {
    geometry->requireCotanLaplacian();
    geometry->requireVertexLumpedMassMatrix();

    std::vector<Vector<double>> evecs = smallestKEigenvectorsPositiveDefinite(
        geometry->cotanLaplacian, geometry->vertexLumpedMassMatrix, 3);

    polyscope::addAnimatedVertexScalarQuantity<Vector<double>>(
        *psMesh, "twoEigs", {evecs[1], evecs[2]});

    psMesh->addVertexScalarQuantity("v1", evecs[1]);
    psMesh->addVertexScalarQuantity("v2", evecs[2]);
}

void computeFlow(bool semiLagrangian) {
    VertexData<double> randomField(*mesh);
    for (Vertex v : mesh->vertices()) {
        randomField[v] = unitRand();
    }

    geometry->requireMeshLengthScale();
    double scale = geometry->meshLengthScale;
    std::vector<VertexData<double>> frames{randomField};
    for (size_t iF = 1; iF < (size_t)steps; ++iF) {
        if (semiLagrangian) {
            frames.push_back(advectSemiLagrangian(
                *mesh, *geometry, velocityField, frames[frames.size() - 1],
                dt / scale / 1000.));
        } else {
            frames.push_back(advect(*mesh, *geometry, velocityField,
                                    frames[frames.size() - 1], dt));
        }
    }


    polyscope::addAnimatedVertexScalarQuantity<VertexData<double>>(
        *psMesh, "flow", frames);
}
// float theta = atan(@P.x, @P.z);
// float ymax  = 0.3;
// float ymin  = -0.4;
// float m     = 20;
// i @inband   = 0;

// if (@P.y >= ymin && @P.y <= ymax) {
//     @psi1re = cos(m * theta);
//     @psi1im = sin(m * theta);
//     @inband = 1;
// } else {
//     @psi1re = 1;
//     @psi1im = 0;
// }
// @psi2re = 0.01;
// @psi2im = 0;

void computeSchroedingerFlow() {
    VertexData<double> randomField(*mesh);
    for (Vertex v : mesh->vertices()) {
        randomField[v] = unitRand();
    }

    SchroedingerSolver ernst(*mesh, *geometry, dt, 0.05);

    Wavefunction psi{VertexData<std::complex<double>>(*mesh),
                     VertexData<std::complex<double>>(*mesh)};
    std::complex<double> i(0, 1);
    VertexData<double> inBand(*mesh);
    for (Vertex v : mesh->vertices()) {
        Vector3 pos  = geometry->inputVertexPositions[v];
        double theta = atan2(pos.x, pos.z);
        double ymax  = 0.3;
        double ymin  = -0.4;
        double m     = 20;
        if (pos.y >= ymin && pos.y <= ymax) {
            psi[0][v] = exp(i * m * theta);
            inBand[v] = 1;
        } else {
            psi[0][v] = 1;
            inBand[v] = 0;
        }
        psi[1][v] = 0.01;
    }
    ernst.normalize(psi);
    psi = ernst.pressureProject(psi);
    ernst.normalize(psi);
    psMesh->addVertexScalarQuantity("inBand", inBand);

    std::vector<Wavefunction> frames = ernst.schroedingerFlow(psi, steps);

    std::vector<FaceData<Vector2>> faceVectors;
    std::vector<VertexData<Vector3>> vertexPhaseColors;
    for (size_t iF = 0; iF < frames.size(); ++iF) {
        faceVectors.push_back(
            ernst.reconstructFaceVectors(ernst.getVelocityField(frames[iF])));
        vertexPhaseColors.push_back(ernst.computeSpin(frames[iF]));
        // psMesh->addFaceIntrinsicVectorQuantity(
        //     "Frame" + std::to_string(iF),
        //     ernst.reconstructFaceVectors(frames[iF]));
    }

    polyscope::addAnimatedFaceIntrinsicVectorQuantity<FaceData<Vector2>>(
        *psMesh, "ISF", faceVectors);
    polyscope::addAnimatedVertexColorQuantity(*psMesh, "ISF phase",
                                              vertexPhaseColors);
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

    if (ImGui::Button("Compute Flow (eulerian)")) {
        computeFlow(false);
    }
    if (ImGui::Button("Compute Flow (semilagrangian)")) {
        computeFlow(true);
    }

    ImGui::SliderFloat("dt", &dt, 0.001, 0.5);
    ImGui::SliderInt("steps", &steps, 1, 500);

    ImGui::Separator();
    if (ImGui::Button("Schroedinger Flow")) {
        computeSchroedingerFlow();
    }
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("geometry-central & Polyscope example project");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // Make sure a mesh name was given
    if (!inputFilename) {
        std::cerr << "Please specify a mesh file as argument" << std::endl;
        return EXIT_FAILURE;
    }

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = myCallback;

    // Load mesh
    std::tie(mesh, geometry) =
        readManifoldSurfaceMesh(args::get(inputFilename));

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
        polyscope::guessNiceNameFromPath(args::get(inputFilename)),
        geometry->inputVertexPositions, mesh->getFaceVertexList(),
        polyscopePermutations(*mesh));

    // Set vertex tangent spaces
    geometry->requireVertexTangentBasis();
    VertexData<Vector3> vBasisX(*mesh);
    for (Vertex v : mesh->vertices()) {
        vBasisX[v] = geometry->vertexTangentBasis[v][0];
    }
    psMesh->setVertexTangentBasisX(vBasisX);

    velocityField =
        geometrycentral::surface::computeSmoothestVertexDirectionField(
            *geometry);
    psMesh->addVertexIntrinsicVectorQuantity("VF", velocityField);

    /*
    std::vector<std::vector<std::array<Vector3, 2>>> gcStreamlines = traceField(
        *mesh, *geometry, toFaces(*mesh, *geometry, velocityField), 1, 1);

    std::vector<Vector3> streamlinePositions;
    std::vector<std::array<size_t, 2>> streamlineEdges;
    for (const std::vector<std::array<Vector3, 2>>& line : gcStreamlines) {

        for (size_t iS = 0; iS + 1 < line.size(); ++iS) {
            const std::array<Vector3, 2>& segment     = line[iS];
            const std::array<Vector3, 2>& nextSegment = line[iS + 1];
            streamlineEdges.push_back(
                {streamlinePositions.size(), streamlinePositions.size() + 1});
            streamlinePositions.push_back(segment[0]);
            streamlinePositions.push_back(nextSegment[0]);
        }
    }

    polyscope::registerCurveNetwork("streamlines", streamlinePositions,
                                    streamlineEdges);

    auto fField = toFaces(*mesh, *geometry, velocityField);
    psMesh->addFaceIntrinsicVectorQuantity("FF", fField);
    */


    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
