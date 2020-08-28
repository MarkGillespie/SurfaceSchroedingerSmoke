#include "flow.h"

// Convert from a vertex-based vector field to face-based by averaging together
// the vectors on all of a face's vertices
FaceData<Vector2> toFaces(ManifoldSurfaceMesh& mesh,
                          VertexPositionGeometry& geo,
                          const VertexData<Vector2>& field) {
    FaceData<Vector2> faceField(mesh, Vector2{0, 0});
    geo.requireHalfedgeVectorsInVertex();
    geo.requireHalfedgeVectorsInFace();

    const HalfedgeData<Vector2>& hV = geo.halfedgeVectorsInVertex;
    const HalfedgeData<Vector2>& hF = geo.halfedgeVectorsInFace;

    for (Face f : mesh.faces()) {
        double D = f.degree();
        for (Halfedge he : f.adjacentHalfedges()) {
            Vector2 v = field[he.vertex()];

            faceField[f] += v / hV[he] * hF[he] / D;
        }
    }

    return faceField;
}

SparseMatrix<double> computeAdvectionMatrix(ManifoldSurfaceMesh& mesh,
                                            VertexPositionGeometry& geo,
                                            const FaceData<Vector2>& velocity) {
    std::vector<Eigen::Triplet<double>> T;
    VertexData<size_t> vIdx = mesh.getVertexIndices();

    geo.requireFaceAreas();
    const FaceData<double>& area = geo.faceAreas;

    geo.requireHalfedgeVectorsInFace();
    const HalfedgeData<Vector2>& hF = geo.halfedgeVectorsInFace;

    Vector2 J{0, 1};
    for (Vertex v : mesh.vertices()) {
        double vArea = 0;
        for (Face f : v.adjacentFaces()) vArea += area[f];
        size_t iV = vIdx[v];
        for (Halfedge he : v.outgoingHalfedges()) {
            size_t iW    = vIdx[he.next().vertex()];
            Face f       = he.face();
            Vector2 grad = (J * hF[he.next()]).normalize();
            double entry = dot(grad, velocity[f]) * area[f] / vArea;
            T.emplace_back(iV, iW, entry);
        }
    }

    size_t n = mesh.nVertices();
    SparseMatrix<double> M(n, n);
    M.setFromTriplets(std::begin(T), std::end(T));

    return M;
}

// Advect with face-based velocity field
VertexData<double> advect(ManifoldSurfaceMesh& mesh,
                          VertexPositionGeometry& geo,
                          const FaceData<Vector2>& velocity,
                          const VertexData<double>& field, double t) {

    Vector<double> fieldVec = field.toVector();

    size_t n = mesh.nVertices();
    SparseMatrix<double> I(n, n);
    I.setIdentity();
    SparseMatrix<double> advectionMatrix =
        computeAdvectionMatrix(mesh, geo, velocity);

    SparseMatrix<double> evolution =
        I + t * advectionMatrix +
        t * t / 2. * advectionMatrix * advectionMatrix;

    Vector<double> newFieldVec = evolution * fieldVec;

    // newFieldVec *= fieldVec.norm() / newFieldVec.norm();

    return VertexData<double>(mesh, newFieldVec);
}

// Advect with vertex-based velocity field
VertexData<double> advect(ManifoldSurfaceMesh& mesh,
                          VertexPositionGeometry& geo,
                          const VertexData<Vector2>& velocity,
                          const VertexData<double>& field, double t) {

    FaceData<Vector2> faceVelocity = toFaces(mesh, geo, velocity);

    return advect(mesh, geo, faceVelocity, field, t);
}

// Advect with vertex-based velocity field
VertexData<double> advectSemiLagrangian(ManifoldSurfaceMesh& mesh,
                                        VertexPositionGeometry& geo,
                                        const VertexData<Vector2>& velocity,
                                        const VertexData<double>& field,
                                        double t) {
    VertexData<double> advectedField(mesh);

    for (Vertex v : mesh.vertices()) {
        TraceGeodesicResult path =
            traceGeodesic(geo, SurfacePoint(v), -t * velocity[v]);

        advectedField[v] = path.endPoint.interpolate(field);
    }

    return advectedField;
}

// Advect with vertex-based velocity field
VertexData<double> advectCombined(ManifoldSurfaceMesh& mesh,
                                  VertexPositionGeometry& geo,
                                  const VertexData<Vector2>& velocity,
                                  const VertexData<double>& field, double t) {
    double s = 0.95;
    return s * advect(mesh, geo, velocity, field, t) +
           (1 - s) * advectSemiLagrangian(mesh, geo, velocity, field, t);
}
