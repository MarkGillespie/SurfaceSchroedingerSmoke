#include "streamlines.h"
#include "polyscope/surface_mesh.h"

// Helpers for tracing
namespace {

struct PointNormal {
    Vector3 p;
    Vector3 n;
};

struct FacePoint {
    Face f;
    Vector3 baryWeights;
};

struct IntersectionResult {
    double tRay;
    double tLine;
};


IntersectionResult rayLineIntersection(Vector2 rayStart, Vector2 rayDir,
                                       Vector2 lineA, Vector2 lineB) {

    Vector2 v1 = rayStart - lineA;
    Vector2 v2 = lineB - lineA;
    Vector2 v3{-rayDir.y, rayDir.x};

    double cross21 = v2.x * v1.y - v2.y * v1.x;
    double tRay    = cross21 / dot(v2, v3);
    double tLine   = dot(v1, v3) / dot(v2, v3);

    if (tRay < 0) {
        tRay = std::numeric_limits<double>::infinity();
    }


    return IntersectionResult{tRay, tLine};
}


Vector3 unitSum(Vector3 v) {
    double sum = v.x + v.y + v.z;
    return v / sum;
}

std::array<Halfedge, 3> faceHedges(Face f) {
    return {f.halfedge(), f.halfedge().next(), f.halfedge().next().next()};
}
std::array<Vertex, 3> faceVertices(Face f) {
    return {f.halfedge().vertex(), f.halfedge().next().vertex(),
            f.halfedge().next().next().vertex()};
}


class FieldTracer {

  public:
    // Core members
    ManifoldSurfaceMesh& mesh;
    VertexPositionGeometry& geo;
    FaceData<Vector2> faceVectors; // disambiguated
    int nSym;

    // Parameters
    double maxLineLength;
    size_t maxFaceCount;

    // Cached geometric quantities
    FaceData<Vector2> vert1InFaceBasis, vert2InFaceBasis;
    double totalArea;

    // Cached connectivity data

    // Input should be identified (raised to power), not disambiguated
    FieldTracer(ManifoldSurfaceMesh& mesh_, VertexPositionGeometry& geo_,
                const FaceData<Vector2>& field, int nSym_ = 1,
                double maxLineLength_ = -1)
        : mesh(mesh_), geo(geo_), nSym(nSym_) {

        faceVectors = FaceData<Vector2>(mesh);
        // Prepare the field
        for (Face f : mesh.faces()) {
            std::complex<double> c{field[f].x, field[f].y};
            std::complex<double> cRoot = std::pow(c, 1.0f / nSym);
            faceVectors[f]             = Vector2{cRoot.real(), cRoot.imag()};
        }

        geo.requireFaceAreas();
        geo.requireFaceTangentBasis();

        vert1InFaceBasis = FaceData<Vector2>(mesh);
        vert2InFaceBasis = FaceData<Vector2>(mesh);

        // Cache geometry quantities
        totalArea = 0;
        for (Face f : mesh.faces()) {
            totalArea += geo.faceAreas[f];

            // Face basis
            Vector3 X = geo.faceTangentBasis[f][0];
            Vector3 Y = geo.faceTangentBasis[f][1];

            const VertexData<Vector3>& pos = geo.inputVertexPositions;

            // Find each of the vertices as a point in the basis
            // The first vertex is implicitly at (0,0)
            Vertex v0    = f.halfedge().vertex();
            Vertex v1    = f.halfedge().next().vertex();
            Vertex v2    = f.halfedge().next().next().vertex();
            Vector3 pos1 = pos[v1] - pos[v0];
            Vector3 pos2 = pos[v2] - pos[v0];
            Vector2 p1{dot(X, pos1), dot(Y, pos1)};
            Vector2 p2{dot(X, pos2), dot(Y, pos2)};

            // Save data
            vert1InFaceBasis[f] = p1;
            vert2InFaceBasis[f] = p2;
        }

        maxLineLength =
            (maxLineLength_ > 0) ? maxLineLength : std::sqrt(totalArea) * .5;
        maxFaceCount = static_cast<size_t>(std::ceil(std::sqrt(mesh.nFaces())));
    }

    Vector3 facePointInR3(FacePoint p) {

        const VertexData<Vector3>& pos = geo.inputVertexPositions;

        Vertex v0 = p.f.halfedge().vertex();
        Vertex v1 = p.f.halfedge().next().vertex();
        Vertex v2 = p.f.halfedge().next().next().vertex();

        return p.baryWeights[0] * pos[v0] + p.baryWeights[1] * pos[v1] +
               p.baryWeights[2] * pos[v2];
    }


    Vector2 barycentricToR2(FacePoint p) {
        return p.baryWeights[1] * vert1InFaceBasis[p.f] +
               p.baryWeights[2] * vert2InFaceBasis[p.f];
    }

    // Trace a single line through the field
    // traceSign should be 1.0 or -1.0, useful for tracing lines backwards
    // through field
    std::vector<std::array<Vector3, 2>>
    traceLine(FacePoint startPoint, Vector2 startDir, double traceSign = 1.0) {

        // Accumulate the result here
        std::vector<std::array<Vector3, 2>> points;

        // Add the initial point
        Vector3 initPoint = facePointInR3(startPoint);
        geo.requireFaceNormals();
        points.push_back({{initPoint, geo.faceNormals[startPoint.f]}});

        // Trace!
        FacePoint currPoint    = startPoint;
        Vector2 currDir        = startDir;
        size_t nFaces          = 0;
        double lengthRemaining = maxLineLength;
        Face prevFace, prevPrevFace;
        while (lengthRemaining > 0 && currPoint.f != prevPrevFace &&
               nFaces < maxFaceCount) {

            nFaces++;
            Face currFace = currPoint.f;

            // Keep track of the last two faces visited
            prevPrevFace = prevFace;
            prevFace     = currFace;

            // Get the data in the basis for this face
            Vector2 v0{0, 0};
            Vector2 v1       = vert1InFaceBasis[currFace];
            Vector2 v2       = vert2InFaceBasis[currFace];
            Vector2 pointPos = barycentricToR2(currPoint);

            // Pick the best symmetric direction in the face
            Vector2 faceDir = faceVectors[currFace];
            Vector2 traceDir{0., 0.};
            double bestAlign = -std::numeric_limits<double>::infinity();
            double deltaRot  = 2.0 * PI / nSym;
            for (int iSym = 0; iSym < nSym; iSym++) {

                double alignScore = dot(traceSign * faceDir, currDir);

                if (alignScore > bestAlign) {
                    bestAlign = alignScore;
                    traceDir  = traceSign * faceDir;
                }

                faceDir = faceDir.rotate(deltaRot);
            }

            // Find when we would exit each edge (if ever)
            IntersectionResult hit0 =
                rayLineIntersection(pointPos, traceDir, v0, v1);
            IntersectionResult hit1 =
                rayLineIntersection(pointPos, traceDir, v1, v2);
            IntersectionResult hit2 =
                rayLineIntersection(pointPos, traceDir, v2, v0);

            std::array<Halfedge, 3> fHalfedges = faceHedges(currFace);

            // Check which edge we would exit first
            Face nextFace;
            Halfedge nextHe, exitHe;
            unsigned int exitHeLocalIndex = -1;
            double tCross, tRay;
            // TODO: wrong cyclic permutation?
            if (hit0.tRay <= hit1.tRay && hit0.tRay <= hit2.tRay) {
                exitHeLocalIndex = 0;
                exitHe           = fHalfedges[0];
                nextHe           = exitHe.twin();
                tCross           = 1.0 - hit0.tLine;
                tRay             = hit0.tRay;
            } else if (hit1.tRay <= hit0.tRay && hit1.tRay <= hit2.tRay) {
                exitHeLocalIndex = 1;
                exitHe           = fHalfedges[1];
                nextHe           = exitHe.twin();
                tCross           = 1.0 - hit1.tLine;
                tRay             = hit1.tRay;
            } else if (hit2.tRay <= hit0.tRay && hit2.tRay <= hit1.tRay) {
                exitHeLocalIndex = 2;
                exitHe           = fHalfedges[2];
                nextHe           = exitHe.twin();
                tCross           = 1.0 - hit2.tLine;
                tRay             = hit2.tRay;
            } else {
                // cerr << "tracing failure :(" << endl;
                return points;
            }

            geo.requireFaceTangentBasis();
            geo.requireFaceNormals();
            // If the ray would end before exiting the face, end it
            if (tRay > lengthRemaining) {
                tRay              = lengthRemaining;
                Vector2 endingPos = pointPos + tRay * traceDir;
                Vector3 endingPosR3 =
                    geo.inputVertexPositions[currFace.halfedge().vertex()] +
                    endingPos.x * geo.faceTangentBasis[currFace][0] +
                    endingPos.y * geo.faceTangentBasis[currFace][1];
                points.push_back({{endingPosR3, geo.faceNormals[currFace]}});
                break;
            }

            // Generate a point for this intersection
            Vector2 newPointLocal = pointPos + tRay * traceDir;
            Vector3 newPointR3 =
                geo.inputVertexPositions[currFace.halfedge().vertex()] +
                newPointLocal.x * geo.faceTangentBasis[currFace][0] +
                newPointLocal.y * geo.faceTangentBasis[currFace][1];
            Vector3 newNormal = geo.faceNormals[currFace];
            if (nextHe != Halfedge()) {
                nextFace  = nextHe.face();
                newNormal = (newNormal + geo.faceNormals[nextFace]).normalize();
            }
            points.push_back({{newPointR3, newNormal}});


            // Quit if we hit a boundary
            if (!nextHe.isInterior()) {
                break;
            }


            // Subtract off the length expended in this face
            lengthRemaining -= tRay;


            // Find a direction of travel in the new face
            currDir = rotateToTangentBasis(traceDir,
                                           geo.faceTangentBasis[currFace][0],
                                           geo.faceTangentBasis[currFace][1],
                                           geo.faceTangentBasis[nextFace][0],
                                           geo.faceTangentBasis[nextFace][1]);

            std::array<Halfedge, 3> nextFaceHedges = faceHedges(nextFace);
            // Figure out which halfedge in the next face is nextHe
            unsigned int nextHeLocalIndex = 0;
            while (nextFaceHedges[nextHeLocalIndex] != nextHe) {
                nextHeLocalIndex++;
            }

            // On a non-manifold / not oriented mesh, the halfedges we're
            // transitioning between might actually point the same way. If so,
            // flip tCross.
            Vertex currTailVert = faceVertices(currFace)[exitHeLocalIndex];
            Vertex nextTailVert = faceVertices(nextFace)[nextHeLocalIndex];
            if (currTailVert == nextTailVert) {
                tCross = 1.0 - tCross;
            }

            // Clamp for safety
            tCross = std::fmin(tCross, 1.0);
            tCross = std::fmax(tCross, 0.);

            // Find barycentric coordinates in the new face
            currPoint = FacePoint{nextFace, {0, 0, 0}};
            currPoint.baryWeights[nextHeLocalIndex]           = 1.0 - tCross;
            currPoint.baryWeights[(nextHeLocalIndex + 1) % 3] = tCross;


            // Pull the result slightly towards the center of the new face to
            // minimize numerical difficulties
            currPoint.baryWeights = unitSum(10000.f * currPoint.baryWeights +
                                            Vector3{1, 1, 1} / 3.f);
        }


        return points;
    }
};


}; // namespace

// Rotate in to a new basis in R3. Vector is rotated in to new tangent plane,
// then a change of basis is performed to the new basis. Basis vectors MUST be
// unit and orthogonal -- this function doesn't check.
Vector2 rotateToTangentBasis(Vector2 v, const Vector3& oldBasisX,
                             const Vector3& oldBasisY, const Vector3& newBasisX,
                             const Vector3& newBasisY) {

    Vector3 oldNormal = cross(oldBasisX, oldBasisY);
    Vector3 newNormal = cross(newBasisX, newBasisY);

    // If the vectors are nearly in plane, no rotation is needed
    // (can't just go through the same code path as below because cross
    // yields a degenerate direction)
    double EPS    = 0.0000001;
    double dotVal = dot(oldNormal, newNormal);
    Vector3 oldBasisInPlaneX, oldBasisInPlaneY;
    if (dotVal > (1.0 - EPS)) {
        // No rotation
        oldBasisInPlaneX = oldBasisX;
        oldBasisInPlaneY = oldBasisY;
    } else if (dotVal < -(1.0 - EPS)) {
        // 180 degree rotation
        oldBasisInPlaneX = -oldBasisX;
        oldBasisInPlaneY = -oldBasisY;
    } else {
        // general rotation
        Vector3 edgeV    = cross(oldNormal, newNormal).normalize();
        double angle     = atan2(dot(edgeV, cross(oldNormal, newNormal)),
                             dot(oldNormal, newNormal));
        oldBasisInPlaneX = oldBasisX.rotateAround(edgeV, angle);
        oldBasisInPlaneY = oldBasisY.rotateAround(edgeV, angle);
        // oldBasisInPlaneX = oldBasisX.rotate_around(edgeV, angle);
        // oldBasisInPlaneY = oldBasisY.rotate_around(edgeV, angle);
    }

    // Now it's just a goood old-fashioned change of basis
    double xComp = v.x * dot(oldBasisInPlaneX, newBasisX) +
                   v.y * dot(oldBasisInPlaneY, newBasisX);
    double yComp = v.x * dot(oldBasisInPlaneX, newBasisY) +
                   v.y * dot(oldBasisInPlaneY, newBasisY);
    return Vector2{xComp, yComp};
}

std::vector<std::vector<std::array<Vector3, 2>>>
traceField(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geo,
           const FaceData<Vector2>& field, int nSym, size_t nLines,
           double maxLineLength) {

    // Only works on triangle meshes
    for (const Face& face : mesh.faces()) {
        if (face.degree() != 3) {
            polyscope::warning("field tracing only supports triangular meshes");
            return std::vector<std::vector<std::array<Vector3, 2>>>();
        }
    }

    // Preliminaries

    // Create a tracer
    FieldTracer tracer(mesh, geo, field, nSym, maxLineLength);

    // Compute a reasonable number of lines if no count was specified
    if (nLines == 0) {
        double lineFactor = 10;
        nLines            = static_cast<size_t>(
            std::ceil(lineFactor * std::sqrt(mesh.nFaces())));
    }

    geo.requireFaceAreas();
    // Shuffle the list of faces to get a reasonable distribution of starting
    // points Build a list of faces to start lines in. Unusually large faces get
    // listed multiple times so we start more lines in them. This roughly
    // approximates a uniform sampling of the mesh. Small faces get oversampled,
    // but that's much less visually striking than large faces getting
    // undersampled.
    std::vector<Face> faceQueue;
    {
        double meanArea = tracer.totalArea / mesh.nFaces();
        for (Face f : mesh.faces()) {
            faceQueue.push_back(f);
            double faceArea = geo.faceAreas[f];
            while (faceArea > meanArea) {
                faceQueue.push_back(f);
                faceArea -= meanArea;
            }
        }

        // Shuffle the list (if we're tracing fewer lines than the size of the
        // list, we want them to be distributed evenly)
        auto randomEngine = std::default_random_engine{};
        std::shuffle(faceQueue.begin(), faceQueue.end(), randomEngine);

        // Make sure the queue of faces to process is long enough by repeating
        // it
        int iAppend = 0;
        while (faceQueue.size() < nLines) {
            faceQueue.push_back(faceQueue[iAppend++]);
        }
    }

    // Unit random numbers
    std::random_device device;
    std::mt19937 mt(device());
    std::uniform_real_distribution<double> unitDist(0.0, 1.0);
    auto unitRand = [&]() { return unitDist(mt); };

    // == Trace the lines
    // cout << "Tracing lines through vector field... " << endl;
    std::vector<std::vector<std::array<Vector3, 2>>> lineList;
    for (size_t i = 0; i < nLines; i++) {

        // Get the next starting face
        Face startFace = faceQueue.back();
        faceQueue.pop_back();

        // Generate a random point in the face
        double r1 = unitRand();
        double r2 = unitRand();
        Vector3 randPoint{1.0 - std::sqrt(r1), std::sqrt(r1) * (1.0 - r2),
                          r2 * std::sqrt(r1)}; // uniform sampling in triangle
        randPoint =
            unitSum(10000.f * randPoint +
                    Vector3{1, 1, 1} / 3.f); // pull slightly towards center

        std::cout << "Starting at face " << startFace.getIndex()
                  << " at coordinates " << randPoint << std::endl;


        // trace half of lines backwards through field, reduces concentration
        // near areas of convergence
        double traceSign = unitRand() > 0.5 ? 1.0 : -1.0;

        // Generate a random direction
        // (the tracing code snaps the velocity to the best-fitting direction,
        // this just serves the role of picking a random direction in symmetric
        // fields)
        Vector2 randomDir =
            (Vector2{unitRand() - .5, unitRand() - .5}).normalize();

        // Trace
        lineList.push_back(tracer.traceLine(FacePoint{startFace, randPoint},
                                            randomDir, traceSign));
    }
    // cout << "    ... done tracing field." << endl;


    return lineList;
}
