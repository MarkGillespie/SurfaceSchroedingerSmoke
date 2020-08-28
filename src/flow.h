#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/vector2.h"

#include "geometrycentral/surface/trace_geodesic.h"

#include "geometrycentral/numerical/linear_algebra_utilities.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// Convert from a vertex-based vector field to face-based by averaging together
// the vectors on all of a face's vertices
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

// Advect with vertex-based velocity field
VertexData<double> advectCombined(ManifoldSurfaceMesh& mesh,
                                  VertexPositionGeometry& geo,
                                  const VertexData<Vector2>& velocity,
                                  const VertexData<double>& field, double t);
