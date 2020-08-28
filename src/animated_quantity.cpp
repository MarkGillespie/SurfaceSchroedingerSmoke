// Copyright 2017-2019, Nicholas Sharp and the Polyscope contributors.
// http://polyscope.run.
#include "animated_quantity.h"

#include "polyscope/file_helpers.h"
#include "polyscope/polyscope.h"
#include "polyscope/render/engine.h"
#include "polyscope/render/shaders.h"

#include "imgui.h"

using std::cout;
using std::endl;

namespace polyscope {

SurfaceAnimatedScalarQuantity::SurfaceAnimatedScalarQuantity(
    std::string name, SurfaceMesh& mesh_, std::string definedOn_,
    DataType dataType_)
    : SurfaceMeshQuantity(name, mesh_, true), dataType(dataType_),
      time(uniquePrefix() + name + "#time", 0),
      cMap(uniquePrefix() + name + "#cmap", defaultColorMap(dataType)),
      definedOn(definedOn_) {}

void SurfaceAnimatedScalarQuantity::draw() {
    if (!isEnabled()) return;

    if (program == nullptr) {
        createProgram();
    }

    // Set uniforms
    parent.setTransformUniforms(*program);
    setProgramUniforms(*program);

    program->draw();
}

void SurfaceAnimatedScalarQuantity::writeToFile(std::string filename) {
    polyscope::warning("Writing to file not yet implemented for this datatype");
}


// Update range uniforms
void SurfaceAnimatedScalarQuantity::setProgramUniforms(
    render::ShaderProgram& program) {
    program.setUniform("u_rangeLow", vizRange.first);
    program.setUniform("u_rangeHigh", vizRange.second);
}


SurfaceAnimatedScalarQuantity* SurfaceAnimatedScalarQuantity::resetMapRange() {
    switch (dataType) {
    case DataType::STANDARD:
        vizRange = dataRange;
        break;
    case DataType::SYMMETRIC: {
        double absRange =
            std::max(std::abs(dataRange.first), std::abs(dataRange.second));
        vizRange = std::make_pair(-absRange, absRange);
    } break;
    case DataType::MAGNITUDE:
        vizRange = std::make_pair(0., dataRange.second);
        break;
    }

    requestRedraw();
    return this;
}

void SurfaceAnimatedScalarQuantity::buildCustomUI() {
    ImGui::SameLine();

    // == Options popup
    if (ImGui::Button("Options")) {
        ImGui::OpenPopup("OptionsPopup");
    }
    if (ImGui::BeginPopup("OptionsPopup")) {

        if (ImGui::MenuItem("Write to file")) writeToFile();
        if (ImGui::MenuItem("Reset colormap range")) resetMapRange();

        ImGui::EndPopup();
    }

    if (render::buildColormapSelector(cMap.get())) {
        program.reset();
        setColorMap(getColorMap());
    }

    // Draw the histogram of values
    hist.colormapRange = vizRange;
    hist.buildUI();

    // Data range
    // Note: %g specifies are generally nicer than %e, but here we don't
    // acutally have a choice. ImGui (for somewhat valid reasons) links the
    // resolution of the slider to the decimal width of the formatted number.
    // When %g formats a number with few decimal places, sliders can break.
    // There is no way to set a minimum number of decimal places with %g,
    // unfortunately.
    {
        switch (dataType) {
        case DataType::STANDARD:
            ImGui::DragFloatRange2("", &vizRange.first, &vizRange.second,
                                   (dataRange.second - dataRange.first) / 100.,
                                   dataRange.first, dataRange.second,
                                   "Min: %.3e", "Max: %.3e");
            break;
        case DataType::SYMMETRIC: {
            float absRange =
                std::max(std::abs(dataRange.first), std::abs(dataRange.second));
            ImGui::DragFloatRange2("##range_symmetric", &vizRange.first,
                                   &vizRange.second, absRange / 100., -absRange,
                                   absRange, "Min: %.3e", "Max: %.3e");
        } break;
        case DataType::MAGNITUDE: {
            ImGui::DragFloatRange2("##range_mag", &vizRange.first,
                                   &vizRange.second, vizRange.second / 100.,
                                   0.0, dataRange.second, "Min: %.3e",
                                   "Max: %.3e");
        } break;
        }
    }

    if (ImGui::SliderFloat("Time", &time.get(), 0.0, 1., "%.5f", 1.)) {
        time.manuallyChanged();
        fillColorBuffers(*program, true);
        requestRedraw();
    }
}

void SurfaceAnimatedScalarQuantity::geometryChanged() { program.reset(); }

SurfaceAnimatedScalarQuantity*
SurfaceAnimatedScalarQuantity::setColorMap(std::string name) {
    cMap = name;
    hist.updateColormap(cMap.get());
    requestRedraw();
    return this;
}
std::string SurfaceAnimatedScalarQuantity::getColorMap() { return cMap.get(); }

SurfaceAnimatedScalarQuantity*
SurfaceAnimatedScalarQuantity::setMapRange(std::pair<double, double> val) {
    vizRange = val;
    requestRedraw();
    return this;
}
std::pair<double, double> SurfaceAnimatedScalarQuantity::getMapRange() {
    return vizRange;
}

SurfaceAnimatedScalarQuantity*
SurfaceAnimatedScalarQuantity::setTime(double t) {
    time.set(t);
    requestRedraw();
    return this;
}
double SurfaceAnimatedScalarQuantity::getTime() { return time.get(); }

std::string SurfaceAnimatedScalarQuantity::niceName() {
    return name + " (" + definedOn + " scalar)";
}

// ========================================================
// ==========           Vertex Scalar            ==========
// ========================================================

SurfaceVertexAnimatedScalarQuantity::SurfaceVertexAnimatedScalarQuantity(
    std::string name, std::vector<std::vector<double>> values_,
    SurfaceMesh& mesh_, DataType dataType_)
    : SurfaceAnimatedScalarQuantity(name, mesh_, "vertex", dataType_),
      values(std::move(values_))

{
    hist.updateColormap(cMap.get());

    // TODO: update histogram over time?
    hist.buildHistogram(values[0], parent.vertexAreas);

    dataRange = robustMinMax(values[0], 1e-5);
    resetMapRange();
}

void SurfaceVertexAnimatedScalarQuantity::createProgram() {
    // Create the program to draw this quantity
    program = render::engine->generateShaderProgram(
        {render::VERTCOLOR_SURFACE_VERT_SHADER,
         render::VERTCOLOR_SURFACE_FRAG_SHADER},
        DrawMode::Triangles);

    // Fill color buffers
    parent.fillGeometryBuffers(*program);
    fillColorBuffers(*program);
    render::engine->setMaterial(*program, parent.getMaterial());
}


void SurfaceVertexAnimatedScalarQuantity::fillColorBuffers(
    render::ShaderProgram& p, bool update) {
    std::vector<double> colorval;
    colorval.reserve(3 * parent.nFacesTriangulation());

    size_t iPrev, iNext;
    double sPrev, sNext;

    std::tie(iPrev, sPrev) = interpolate(getTime());
    iNext                  = iPrev + 1;
    sNext                  = 1 - sPrev;

    auto interpolatedValue = [&](size_t iV) {
        double valPrev = values[iPrev][iV];
        double valNext = values[iNext][iV];

        return sPrev * valPrev + sNext * valNext;
    };

    for (size_t iF = 0; iF < parent.nFaces(); iF++) {
        auto& face = parent.faces[iF];
        size_t D   = face.size();

        // implicitly triangulate from root
        size_t vRoot = face[0];
        for (size_t j = 1; (j + 1) < D; j++) {
            size_t vB = face[j];
            size_t vC = face[(j + 1) % D];

            colorval.push_back(interpolatedValue(vRoot));
            colorval.push_back(interpolatedValue(vB));
            colorval.push_back(interpolatedValue(vC));
        }
    }

    // Store data in buffers
    p.setAttribute("a_colorval", colorval);

    if (!update) p.setTextureFromColormap("t_colormap", cMap.get());
}

void SurfaceVertexAnimatedScalarQuantity::writeToFile(std::string filename) {

    throw std::runtime_error("not implemented");

    /* TODO
    if (filename == "") {
      filename = promptForFilename();
      if (filename == "") {
        return;
      }
    }

    // For now, just always write scalar to U texture coordinate

    cout << "Writing vertex value to file " << filename << " in U coordinate of
    texture map" << endl;

    HalfedgeMesh* mesh = parent.mesh;
    CornerData<Vector2> scalarVal(mesh, Vector2{0.0, 0.0});
    for (CornerPtr c : mesh->corners()) {
      scalarVal[c].x = values[c.vertex()];
    }

    WavefrontOBJ::write(filename, *parent.geometry, scalarVal);
    */
}

void SurfaceVertexAnimatedScalarQuantity::buildVertexInfoGUI(size_t vInd) {
    ImGui::TextUnformatted(name.c_str());
    ImGui::NextColumn();

    size_t iPrev, iNext;
    double sPrev, sNext;
    std::tie(iPrev, sPrev) = interpolate(getTime());
    iNext                  = iPrev + 1;
    sNext                  = 1 - sPrev;

    auto interpolatedValue = [&](size_t iV) {
        double valPrev = values[iPrev][iV];
        double valNext = values[iNext][iV];

        return sPrev * valPrev + sNext * valNext;
    };
    ImGui::Text("%g", interpolatedValue(vInd));
    ImGui::NextColumn();
}

std::tuple<size_t, double>
SurfaceVertexAnimatedScalarQuantity::interpolate(double t) {
    if (t >= 1 - 1e-8) {
        return std::make_tuple(values.size() - 2, 0);
    } else if (t <= 1e-8) {
        return std::make_tuple(0, 1);
    }
    double N  = values.size();
    size_t iF = floor(t * (N - 1));

    double remainder = t - ((double)iF) / N;

    // std::cout << "interpolating time " << t << " computed iF = " << iF
    //           << "   and   remainder = " << remainder << std::endl;

    return std::make_tuple(iF, 1 - remainder * (N - 1));
}

SurfaceVertexAnimatedScalarQuantity* addAnimatedVertexScalarQuantityImpl(
    SurfaceMesh& mesh, std::string name,
    const std::vector<std::vector<double>>& data, DataType type) {

    std::vector<std::vector<double>> permutedData;

    for (const std::vector<double>& frame : data) {
        permutedData.push_back(applyPermutation(frame, mesh.vertexPerm));
    }


    SurfaceVertexAnimatedScalarQuantity* q =
        new SurfaceVertexAnimatedScalarQuantity(name, permutedData, mesh, type);
    mesh.addQuantity(q);
    return q;
}
} // namespace polyscope
