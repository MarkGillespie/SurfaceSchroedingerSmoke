#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


// Trace lines through a vector field on a mesh.
// Return is a list of lines, each entry is (position, normal)
// Input field should be identified (raised to power), not disambiguated
// Settings 0 for nLines results in an automatically computed value
std::vector<std::vector<std::array<Vector3, 2>>>
traceField(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo,
           const FaceData<Vector2>& field, int nSym = 1, size_t nLines = 0,
           double maxLineLength_ = -1);

// Rotate in to a new basis in R3. Vector is rotated in to new tangent plane,
// then a change of basis is performed to the new basis. Basis vectors MUST be
// unit and orthogonal -- this function doesn't check.
Vector2 rotateToTangentBasis(Vector2 v, const Vector3& oldBasisX,
                             const Vector3& oldBasisY, const Vector3& newBasisX,
                             const Vector3& newBasisY);
