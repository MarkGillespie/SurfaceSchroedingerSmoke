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

SparseMatrix<std::complex<double>> connectionD0(ManifoldSurfaceMesh& mesh,
                                                VertexPositionGeometry& geo) {
    std::vector<Eigen::Triplet<std::complex<double>>> T;

    geo.requireTransportVectorsAlongHalfedge();
    HalfedgeData<Vector2> hedgeRot = geo.transportVectorsAlongHalfedge;

    VertexData<size_t> vIdx = mesh.getVertexIndices();
    EdgeData<size_t> eIdx   = mesh.getEdgeIndices();

    std::complex<double> i(0, 1);
    for (Edge e : mesh.edges()) {
        Vertex v = e.halfedge().vertex();
        Vertex w = e.halfedge().next().vertex();
        std::complex<double> rot =
            static_cast<std::complex<double>>(hedgeRot[e.halfedge()]);
        std::complex<double> halfRot = std::sqrt(rot);

        T.emplace_back(eIdx[e], vIdx[w], halfRot);
        T.emplace_back(eIdx[e], vIdx[w], -std::conj(halfRot));
    }

    SparseMatrix<std::complex<double>> d0(mesh.nEdges(), mesh.nVertices());
    d0.setFromTriplets(std::begin(T), std::end(T));
    return d0;
}

SchroedingerSolver::SchroedingerSolver(ManifoldSurfaceMesh& mesh_,
                                       VertexPositionGeometry& geo_, double dt,
                                       double hbar_)
    : mesh(mesh_), geo(geo_), hbar(hbar_) {

    geo.requireCotanLaplacian();
    geo.requireVertexConnectionLaplacian();
    geo.requireVertexLumpedMassMatrix();

    SparseMatrix<double> L = geo.cotanLaplacian;
    M = geo.vertexLumpedMassMatrix.cast<std::complex<double>>();

    /*
    SparseMatrix<std::complex<double>> Lconn = geo.vertexConnectionLaplacian;
    SparseMatrix<std::complex<double>> d0conn = connectionD0(mesh, geo);

    geo.requireDECOperators();
    SparseMatrix<double> DEC_L = geo.d0.transpose() * geo.hodge1 * geo.d0;
    for (size_t iV = 0; iV < mesh.nVertices(); ++iV) {
        Vector<std::complex<double>> vecI =
            Vector<std::complex<double>>::Zero(mesh.nVertices());
        vecI(iV) = 1;

        Vector<std::complex<double>> DECLV =
            d0conn.transpose() * geo.hodge1 * d0conn * vecI;
        Vector<std::complex<double>> LV = Lconn * vecI;

        if ((DECLV - LV).norm() > 1e-8) {
            std::cerr << "DEC laplacian and cotan laplacian disgree?"
                      << std::endl;
            exit(1);
        }

          Vector<double> vecI = Vector<double>::Zero(mesh.nVertices());
          vecI(iV)            = 1;

          Vector<double> DECLV = DEC_L * vecI;
          Vector<double> LV    = L * vecI;

          if ((DECLV - LV).norm() > 1e-8) {
              std::cerr << "DEC laplacian and cotan laplacian disgree?"
                        << std::endl;
              exit(1);
          }
  }
  */

    std::complex<double> i(0, 1);
    SparseMatrix<std::complex<double>> evolution0 =
        M.cast<std::complex<double>>() - i * hbar * dt / 4. * L;
    SparseMatrix<std::complex<double>> evolution1 =
        M.cast<std::complex<double>>() + i * hbar * dt / 4. * L;

    evolutionSolver0 = std::unique_ptr<SquareSolver<std::complex<double>>>(
        new SquareSolver<std::complex<double>>(evolution0));
    evolutionSolver1 = std::unique_ptr<SquareSolver<std::complex<double>>>(
        new SquareSolver<std::complex<double>>(evolution1));
    Lsolver = std::unique_ptr<PositiveDefiniteSolver<double>>(
        new PositiveDefiniteSolver<double>(L));
}

Wavefunction SchroedingerSolver::schroedingerAdvect(const Wavefunction& psi) {

    Vector<std::complex<double>> psi0 = psi[0].toVector();
    Vector<std::complex<double>> psi1 = psi[1].toVector();

    Vector<std::complex<double>> newPsi0 = evolutionSolver0->solve(M * psi0);
    Vector<std::complex<double>> newPsi1 = evolutionSolver1->solve(M * psi1);

    return Wavefunction{VertexData<std::complex<double>>(mesh, newPsi0),
                        VertexData<std::complex<double>>(mesh, newPsi1)};
}

void SchroedingerSolver::normalize(Wavefunction& psi) {
    for (Vertex v : mesh.vertices()) {
        // Apparently std::norm(std::complex) computes the *squared* norm
        double vNorm = sqrt(std::norm(psi[0][v]) + std::norm(psi[1][v]));
        psi[0][v] /= vNorm;
        psi[1][v] /= vNorm;
    }
}

VertexData<double>
SchroedingerSolver::computePressure(const EdgeData<double>& velocity) {

    geo.requireDECOperators();
    Vector<double> delV = geo.d0.transpose() * geo.hodge1 * velocity.toVector();

    Vector<double> pressureVec = Lsolver->solve(delV);

    Vector<double> updatedVelocity = velocity.toVector() - geo.d0 * pressureVec;
    Vector<double> updatedDivergence =
        geo.d0.transpose() * geo.hodge1 * updatedVelocity;

    for (size_t iV = 0; iV < mesh.nVertices(); ++iV) {
        if (abs(updatedDivergence(iV)) > 1e-8) {
            std::cerr << "Pressure projection failed (inside)" << std::endl;
            exit(1);
        }
    }

    return VertexData<double>(mesh, pressureVec);
}

Wavefunction SchroedingerSolver::pressureProject(Wavefunction psi) {
    EdgeData<double> velocityField = getVelocityField(psi);

    VertexData<double> pressure = computePressure(velocityField);
    Wavefunction oldPsi         = psi;

    std::complex<double> i(0, 1);
    for (Vertex v : mesh.vertices()) {
        psi[0][v] *= exp(-i * pressure[v] / hbar);
        psi[1][v] *= exp(-i * pressure[v] / hbar);
    }

    /*
    EdgeData<double> updatedVelocityField = getVelocityField(psi);

    geo.requireDECOperators();
    Vector<double> edgeUpdate =
        updatedVelocityField.toVector() - velocityField.toVector();
    Vector<double> gradP = geo.d0 * pressure.toVector();

    EdgeData<size_t> eIdx = mesh.getEdgeIndices();
    for (Edge e : mesh.edges()) {
        size_t iE = eIdx[e];
        Vertex v  = e.halfedge().vertex();
        Vertex w  = e.halfedge().next().vertex();
        if (abs(edgeUpdate(iE) + gradP(iE)) > 1e-8) {
            std::cerr << "Error: gradP = " << gradP(iE)
                      << "  but edge update was " << edgeUpdate(iE)
                      << "\t| new velocity: " << updatedVelocityField[e]
                      << ", old velocity: " << velocityField[e]
                      << "\t| err / hbar = "
                      << (edgeUpdate(iE) + gradP(iE)) / hbar << std::endl;
            // exit(1);
        }
    }

    Vector<double> divergence =
        geo.d0.transpose() * geo.hodge1 * updatedVelocityField.toVector();

    for (size_t iV = 0; iV < mesh.nVertices(); ++iV) {
        if (abs(divergence(iV)) > 1e-8) {
            std::cerr << "Pressure projection failed" << std::endl;
            exit(1);
        }
    }
    */

    return psi;
}

// v@s;
// @s.x = 2*(@psi1re*@psi2re - @psi1im*@psi2im);
// @s.y = 2*(@psi1re*@psi2im + @psi2re*@psi1im);
// @s.z = pow(@psi1re,2)+pow(@psi1im,2)
//                          -pow(@psi2re,2)-pow(@psi2im,2);
VertexData<Vector3> SchroedingerSolver::computeSpin(Wavefunction psi) {
    VertexData<Vector3> spins(mesh);
    for (Vertex v : mesh.vertices()) {
        std::complex<double> prod = psi[0][v] * psi[1][v];

        spins[v] = Vector3{2 * std::real(prod), 2 * std::imag(prod),
                           std::norm(psi[0][v]) - std::norm(psi[1][v])};
        spins[v] = spins[v] / 2. + Vector3{0.5, 0.5, 0.5};
    }
    return spins;
}

std::vector<Wavefunction>
SchroedingerSolver::schroedingerFlow(const Wavefunction& psi, size_t steps) {
    std::vector<Wavefunction> frames{psi};

    for (size_t iS = 0; iS < steps; ++iS) {
        const Wavefunction oldPsi = frames[frames.size() - 1];

        Wavefunction psiTmp = schroedingerAdvect(oldPsi);
        normalize(psiTmp);
        frames.push_back(pressureProject(psiTmp));
    }

    return frames;
}

std::vector<EdgeData<double>>
SchroedingerSolver::schroedingerFlowVelocities(const Wavefunction& psi,
                                               size_t steps) {
    std::vector<EdgeData<double>> velocities;
    for (Wavefunction psi_t : schroedingerFlow(psi, steps)) {
        velocities.push_back(getVelocityField(psi_t));
    }

    return velocities;
}

EdgeData<double> SchroedingerSolver::getVelocityField(const Wavefunction& psi) {
    EdgeData<double> velocity(mesh, 0);

    geo.requireTransportVectorsAlongHalfedge();
    HalfedgeData<Vector2> hedgeRot = geo.transportVectorsAlongHalfedge;
    std::complex<double> i(0, 1);
    for (Edge e : mesh.edges()) {
        Vertex v = e.halfedge().vertex();
        Vertex w = e.halfedge().next().vertex();
        // std::complex<double> rot =
        //     static_cast<std::complex<double>>(hedgeRot[e.halfedge()]);
        // std::complex<double> psi_v_0 = rot * psi[0][v];
        // std::complex<double> psi_v_1 = rot * psi[1][v];
        std::complex<double> psi_v_0 = psi[0][v];
        std::complex<double> psi_v_1 = psi[1][v];
        // std::complex<double> psi_v_0 = rot * psi[0][v];
        // std::complex<double> psi_v_1 = rot * psi[1][v];
        std::complex<double> psi_w_0 = psi[0][w];
        std::complex<double> psi_w_1 = psi[1][w];

        std::complex<double> psi_v_bar_psi_w =
            std::conj(psi_v_0) * psi_w_0 + std::conj(psi_v_1) * psi_w_1;
        velocity[e] = hbar * std::arg(psi_v_bar_psi_w);

        // double a    = std::real(psi_v_bar_psi_w);
        // double b    = std::imag(psi_v_bar_psi_w);
        // velocity[e] = hbar * atan2(b, a);

        // std::complex<double> psi_avg_0 = (psi_v_0 + psi_w_0) / 2.;
        // std::complex<double> psi_avg_1 = (psi_v_1 + psi_w_1) / 2.;
        // double r0_sq                   = std::norm(psi_avg_0);
        // double r1_sq                   = std::norm(psi_avg_1);

        // double q = r0_sq + r1_sq;

        // std::complex<double> dpsi_0 = psi_w_0 / psi_v_0;
        // std::complex<double> dpsi_1 = psi_w_0 / psi_v_0;

        // velocity[e] =
        //     hbar * (r0_sq * std::arg(dpsi_0) + r1_sq * std::arg(dpsi_1)) / q;

        // velocity[e] = hbar / q *
        //               std::real(-std::conj(psi_avg_0) * i * dpsi_0 -
        //                         std::conj(psi_avg_1) * i * dpsi_1);


        // velocity[e] = hbar * std::arg(std::conj(psi_v_0) * psi_w_0 +
        //                               std::conj(psi_v_1) * psi_w_1);
    }

    return velocity;
}

FaceData<Vector2> SchroedingerSolver::reconstructFaceVectors(
    const EdgeData<double>& velocityField) {

    geo.requireHalfedgeVectorsInFace();
    const HalfedgeData<Vector2>& faceVec = geo.halfedgeVectorsInFace;

    auto hedgeVelocity = [&](Halfedge he) {
        double sign = (he.edge().halfedge() == he) ? 1 : -1;
        return sign * velocityField[he.edge()];
    };

    FaceData<Vector2> faceVectors(mesh);
    for (Face f : mesh.faces()) {
        Vector2 v0 = faceVec[f.halfedge()];
        Vector2 v1 = faceVec[f.halfedge().next()];
        double a0  = hedgeVelocity(f.halfedge());
        double a1  = hedgeVelocity(f.halfedge().next());

        double det = v0.x * v1.y - v0.y * v1.x;
        double w0  = (v1.y * a0 - v0.y * a1) / det;
        double w1  = (-v1.x * a0 + v0.x * a1) / det;

        faceVectors[f] = {w0, w1};
    }

    return faceVectors;
}
