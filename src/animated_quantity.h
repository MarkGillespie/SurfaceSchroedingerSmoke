#pragma once

#include "polyscope/affine_remapper.h"
#include "polyscope/histogram.h"
#include "polyscope/render/color_maps.h"
#include "polyscope/render/engine.h"
#include "polyscope/surface_mesh.h"

namespace polyscope {

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

class SurfaceVertexAnimatedScalarQuantity
    : public SurfaceAnimatedScalarQuantity {
  public:
    SurfaceVertexAnimatedScalarQuantity(
        std::string name, std::vector<std::vector<double>> values_,
        SurfaceMesh& mesh_, DataType dataType_ = DataType::STANDARD);

    virtual void createProgram() override;

    virtual void fillColorBuffers(render::ShaderProgram& p,
                                  bool update = false) override;

    void buildVertexInfoGUI(size_t vInd) override;
    virtual void writeToFile(std::string filename = "") override;


    // Assuming that our values are uniformly spaced keyframes in [0, 1],
    // returns the first frame before t, and the coefficient to interpolate
    // along the way to the next frame
    std::tuple<size_t, double> interpolate(double t);

    // === Members
    std::vector<std::vector<double>> values;
};

SurfaceVertexAnimatedScalarQuantity* addAnimatedVertexScalarQuantityImpl(
    SurfaceMesh& mesh, std::string name,
    const std::vector<std::vector<double>>& data,
    DataType type = DataType::STANDARD);


template <typename T>
SurfaceVertexAnimatedScalarQuantity*
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

} // namespace polyscope
