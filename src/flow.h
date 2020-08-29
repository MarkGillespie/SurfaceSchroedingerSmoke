#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/vector2.h"

#include "geometrycentral/surface/trace_geodesic.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/numerical/linear_solvers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

using Wavefunction = std::array<VertexData<std::complex<double>>, 2>;

// Convert from a vertex-based vector field to face-based by averaging
// together the vectors on all of a face's vertices
FaceData<Vector2> toFaces(ManifoldSurfaceMesh& mesh,
                          VertexPositionGeometry& geo,
                          const VertexData<Vector2>& field);

SparseMatrix<double> computeAdvectionMatrix(ManifoldSurfaceMesh& mesh,
                                            VertexPositionGeometry& geo,
                                            const FaceData<Vector2>& velocity);

// Advect with face-based velocity field
VertexData<double> advect(ManifoldSurfaceMesh& mesh,
                          VertexPositionGeometry& geo,
                          const FaceData<Vector2>& velocity,
                          const VertexData<double>& field, double t);

// Advect with vertex-based velocity field
VertexData<double> advect(ManifoldSurfaceMesh& mesh,
                          VertexPositionGeometry& geo,
                          const VertexData<Vector2>& velocity,
                          const VertexData<double>& field, double t);

// Advect with vertex-based velocity field
VertexData<double> advectSemiLagrangian(ManifoldSurfaceMesh& mesh,
                                        VertexPositionGeometry& geo,
                                        const VertexData<Vector2>& velocity,
                                        const VertexData<double>& field,
                                        double t);

SparseMatrix<std::complex<double>> connectionD0(ManifoldSurfaceMesh& mesh);

class SchroedingerSolver {
  public:
    SchroedingerSolver(ManifoldSurfaceMesh& mesh_, VertexPositionGeometry& geo_,
                       double dt, double hbar_);

    std::vector<Wavefunction> schroedingerFlow(const Wavefunction& psi,
                                               size_t steps);

    std::vector<EdgeData<double>>
    schroedingerFlowVelocities(const Wavefunction& psi, size_t steps);

    EdgeData<double> getVelocityField(const Wavefunction& psi);

    FaceData<Vector2>
    reconstructFaceVectors(const EdgeData<double>& velocityField);

    void normalize(Wavefunction& psi);
    Wavefunction schroedingerAdvect(const Wavefunction& psi);
    VertexData<double> computePressure(const EdgeData<double>& velocity);
    Wavefunction pressureProject(Wavefunction psi);
    VertexData<Vector3> computeSpin(Wavefunction psi);

  protected:
    ManifoldSurfaceMesh& mesh;
    VertexPositionGeometry& geo;
    double hbar;

    SparseMatrix<std::complex<double>> M;
    std::unique_ptr<SquareSolver<std::complex<double>>> evolutionSolver0;
    std::unique_ptr<SquareSolver<std::complex<double>>> evolutionSolver1;
    std::unique_ptr<PositiveDefiniteSolver<double>> Lsolver;
};
