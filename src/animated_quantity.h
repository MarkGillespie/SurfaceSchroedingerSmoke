#pragma once

#include "polyscope/affine_remapper.h"
#include "polyscope/histogram.h"
#include "polyscope/render/color_maps.h"
#include "polyscope/render/engine.h"
#include "polyscope/surface_mesh.h"

namespace polyscope {

// Assuming that our values are uniformly spaced keyframes in [0, 1],
// returns the first frame before t, and the coefficient to interpolate
// along the way to the next frame
std::tuple<size_t, double> interpolate(double t, size_t N);

class SurfaceAnimatedScalarQuantity : public SurfaceMeshQuantity {
  public:
    SurfaceAnimatedScalarQuantity(std::string name, SurfaceMesh& mesh_,
                                  std::string definedOn, DataType dataType);

    virtual void draw() override;
    virtual void buildCustomUI() override;
    virtual std::string niceName() override;
    virtual void geometryChanged() override;

    virtual void writeToFile(std::string filename = "");

    // === Members
    const DataType dataType;

    // === Get/set visualization parameters

    // The color map
    SurfaceAnimatedScalarQuantity* setColorMap(std::string val);
    std::string getColorMap();

    // Data limits mapped in to colormap
    SurfaceAnimatedScalarQuantity* setMapRange(std::pair<double, double> val);
    std::pair<double, double> getMapRange();
    SurfaceAnimatedScalarQuantity* resetMapRange(); // reset to full range

    SurfaceAnimatedScalarQuantity* setTime(double t);
    double getTime();

  protected:
    // === Visualization parameters

    // Affine data maps and limits
    std::pair<float, float> vizRange;
    std::pair<double, double> dataRange;
    Histogram hist;

    // UI internals
    PersistentValue<float> time;
    PersistentValue<std::string> cMap;
    const std::string definedOn;
    std::shared_ptr<render::ShaderProgram> program;

    // Helpers
    virtual void fillColorBuffers(render::ShaderProgram& p,
                                  bool update = false) = 0;
    virtual void createProgram()                       = 0;
    void setProgramUniforms(render::ShaderProgram& program);
};


// ========================================================
// ==========           Vertex Scalar            ==========
// ========================================================

class SurfaceAnimatedVertexScalarQuantity
    : public SurfaceAnimatedScalarQuantity {
  public:
    SurfaceAnimatedVertexScalarQuantity(
        std::string name, std::vector<std::vector<double>> values_,
        SurfaceMesh& mesh_, DataType dataType_ = DataType::STANDARD);

    virtual void createProgram() override;

    virtual void fillColorBuffers(render::ShaderProgram& p,
                                  bool update = false) override;

    void buildVertexInfoGUI(size_t vInd) override;
    virtual void writeToFile(std::string filename = "") override;


    // === Members
    std::vector<std::vector<double>> values;
};

SurfaceAnimatedVertexScalarQuantity* addAnimatedVertexScalarQuantityImpl(
    SurfaceMesh& mesh, std::string name,
    const std::vector<std::vector<double>>& data,
    DataType type = DataType::STANDARD);


template <typename T>
SurfaceAnimatedVertexScalarQuantity*
addAnimatedVertexScalarQuantity(SurfaceMesh& mesh, std::string name,
                                const std::vector<T>& data,
                                DataType type = DataType::STANDARD) {
    std::vector<std::vector<double>> standardizedData;
    for (const T& frame : data) {
        validateSize(frame, mesh.vertexDataSize,
                     "vertex scalar quantity " + name);
        standardizedData.push_back(standardizeArray<double, T>(frame));
    }

    return addAnimatedVertexScalarQuantityImpl(mesh, name, standardizedData,
                                               type);
}

// ==== Common base class

// Represents a general vector field associated with a surface mesh, including
// R3 fields in the ambient space and R2 fields embedded in the surface
class SurfaceAnimatedVectorQuantity : public SurfaceMeshQuantity {
  public:
    SurfaceAnimatedVectorQuantity(
        std::string name, SurfaceMesh& mesh_, MeshElement definedOn_,
        VectorType vectorType_ = VectorType::STANDARD);


    virtual void draw() override;
    virtual void buildCustomUI() override;

    // Allow children to append to the UI
    virtual void drawSubUI();

    // === Members

    // Note: these vectors are not the raw vectors passed in by the user, but
    // have been rescaled such that the longest has length 1 (unless type is
    // VectorType::Ambient)
    const VectorType vectorType;
    std::vector<glm::vec3> vectorRoots;
    std::vector<glm::vec3> vectors;

    // === Option accessors

    //  The vectors will be scaled such that the longest vector is this long
    SurfaceAnimatedVectorQuantity* setVectorLengthScale(double newLength,
                                                        bool isRelative = true);
    double getVectorLengthScale();

    // The radius of the vectors
    SurfaceAnimatedVectorQuantity* setVectorRadius(double val,
                                                   bool isRelative = true);
    double getVectorRadius();

    // The color of the vectors
    SurfaceAnimatedVectorQuantity* setVectorColor(glm::vec3 color);
    glm::vec3 getVectorColor();

    // Material
    SurfaceAnimatedVectorQuantity* setMaterial(std::string name);
    std::string getMaterial();

    // Enable the ribbon visualization
    SurfaceAnimatedVectorQuantity* setRibbonEnabled(bool newVal);
    bool isRibbonEnabled();

    virtual void recomputeVectors() = 0;

    SurfaceAnimatedVectorQuantity* setTime(double t);
    double getTime();

  protected:
    // === Visualization options
    PersistentValue<ScaledValue<float>> vectorLengthMult;
    PersistentValue<ScaledValue<float>> vectorRadius;
    PersistentValue<glm::vec3> vectorColor;
    PersistentValue<std::string> material;
    PersistentValue<float> time;

    // The map that takes values to [0,1] for drawing
    AffineRemapper<glm::vec3> mapper;

    MeshElement definedOn;

    // A ribbon viz that is appropriate for some fields
    std::unique_ptr<RibbonArtist> ribbonArtist;
    PersistentValue<bool> ribbonEnabled;

    // GL things
    void prepareProgram();
    std::shared_ptr<render::ShaderProgram> program;

    // Set up the mapper for vectors
    void prepareVectorMapper();
};

// ==== Intrinsic vectors at faces

class SurfaceAnimatedFaceIntrinsicVectorQuantity
    : public SurfaceAnimatedVectorQuantity {
  public:
    SurfaceAnimatedFaceIntrinsicVectorQuantity(
        std::string name, std::vector<std::vector<glm::vec2>> vectors_,
        SurfaceMesh& mesh_, int nSym = 1,
        VectorType vectorType_ = VectorType::STANDARD);

    int nSym;
    std::vector<std::vector<glm::vec2>> vectorField;

    virtual void recomputeVectors() override;

    virtual void draw() override;

    void drawSubUI() override;

    virtual std::string niceName() override;
    void buildFaceInfoGUI(size_t fInd) override;
};

SurfaceAnimatedFaceIntrinsicVectorQuantity*
addAnimatedFaceIntrinsicVectorQuantityImpl(
    SurfaceMesh& mesh, std::string name,
    const std::vector<std::vector<glm::vec2>>& data, int nSym = 1,
    VectorType type = VectorType::STANDARD);


template <typename T>
SurfaceAnimatedFaceIntrinsicVectorQuantity*
addAnimatedFaceIntrinsicVectorQuantity(SurfaceMesh& mesh, std::string name,
                                       const std::vector<T>& data, int nSym = 1,
                                       VectorType type = VectorType::STANDARD) {
    std::vector<std::vector<glm::vec2>> standardizedData;
    for (const T& frame : data) {
        validateSize(frame, mesh.faceDataSize,
                     "face intrinsic vector quantity " + name);
        standardizedData.push_back(standardizeVectorArray<glm::vec2, 2>(frame));
    }

    return addAnimatedFaceIntrinsicVectorQuantityImpl(
        mesh, name, standardizedData, nSym, type);
}

class SurfaceAnimatedColorQuantity : public SurfaceMeshQuantity {
  public:
    SurfaceAnimatedColorQuantity(std::string name, SurfaceMesh& mesh_,
                                 std::string definedOn);

    virtual void draw() override;
    virtual void buildCustomUI() override;
    virtual std::string niceName() override;

    virtual void geometryChanged() override;

    SurfaceAnimatedColorQuantity* setTime(double t);
    double getTime();

  protected:
    // UI internals
    PersistentValue<float> time;
    const std::string definedOn;
    std::shared_ptr<render::ShaderProgram> program;

    // Helpers
    virtual void createProgram() = 0;

    virtual void fillColorBuffers(render::ShaderProgram& p) = 0;
};

// ========================================================
// ==========           Vertex Color             ==========
// ========================================================

class SurfaceAnimatedVertexColorQuantity : public SurfaceAnimatedColorQuantity {
  public:
    SurfaceAnimatedVertexColorQuantity(
        std::string name, std::vector<std::vector<glm::vec3>> values_,
        SurfaceMesh& mesh_);

    virtual void createProgram() override;
    virtual void fillColorBuffers(render::ShaderProgram& p) override;

    void buildVertexInfoGUI(size_t vInd) override;

    // === Members
    std::vector<std::vector<glm::vec3>> values;
};

SurfaceAnimatedVertexColorQuantity* addAnimatedVertexColorQuantityImpl(
    SurfaceMesh& mesh, std::string name,
    const std::vector<std::vector<glm::vec3>>& data);


template <typename T>
SurfaceAnimatedVertexColorQuantity*
addAnimatedVertexColorQuantity(SurfaceMesh& mesh, std::string name,
                               const std::vector<T>& data) {
    std::vector<std::vector<glm::vec3>> standardizedData;
    for (const T& frame : data) {
        validateSize(frame, mesh.vertexDataSize,
                     "vertex color quantity " + name);
        standardizedData.push_back(standardizeVectorArray<glm::vec3, 3>(frame));
    }

    return addAnimatedVertexColorQuantityImpl(mesh, name, standardizedData);
}

} // namespace polyscope
