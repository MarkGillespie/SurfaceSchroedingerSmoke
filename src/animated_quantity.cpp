// Copyright 2017-2019, Nicholas Sharp and the Polyscope contributors.
// http://polyscope.run.
#include "animated_quantity.h"

#include "polyscope/file_helpers.h"
#include "polyscope/polyscope.h"
#include "polyscope/render/engine.h"
#include "polyscope/render/shaders.h"

#include "polyscope/trace_vector_field.h"

#include "imgui.h"

#include <complex>
#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

namespace polyscope {

std::tuple<size_t, double> interpolate(double t, size_t N) {
    if (t >= 1 - 1e-8) {
        return std::make_tuple(N - 2, 0);
    } else if (t <= 1e-8) {
        return std::make_tuple(0, 1);
    }
    size_t iF = floor(t * (double)(N - 1));

    double remainder = t - ((double)iF) / ((double)(N - 1));

    // std::cout << "interpolating time " << t << " computed iF = " << iF
    //           << "   and   remainder = " << remainder << std::endl;

    return std::make_tuple(iF, 1 - remainder * (double)(N - 1));
}


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

SurfaceAnimatedVertexScalarQuantity::SurfaceAnimatedVertexScalarQuantity(
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

void SurfaceAnimatedVertexScalarQuantity::createProgram() {
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


void SurfaceAnimatedVertexScalarQuantity::fillColorBuffers(
    render::ShaderProgram& p, bool update) {
    std::vector<double> colorval;
    colorval.reserve(3 * parent.nFacesTriangulation());

    size_t iPrev, iNext;
    double sPrev, sNext;

    std::tie(iPrev, sPrev) = interpolate(getTime(), values.size());
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

void SurfaceAnimatedVertexScalarQuantity::writeToFile(std::string filename) {

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

void SurfaceAnimatedVertexScalarQuantity::buildVertexInfoGUI(size_t vInd) {
    ImGui::TextUnformatted(name.c_str());
    ImGui::NextColumn();

    size_t iPrev, iNext;
    double sPrev, sNext;
    std::tie(iPrev, sPrev) = interpolate(getTime(), values.size());
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

SurfaceAnimatedVertexScalarQuantity* addAnimatedVertexScalarQuantityImpl(
    SurfaceMesh& mesh, std::string name,
    const std::vector<std::vector<double>>& data, DataType type) {

    std::vector<std::vector<double>> permutedData;

    for (const std::vector<double>& frame : data) {
        permutedData.push_back(applyPermutation(frame, mesh.vertexPerm));
    }


    SurfaceAnimatedVertexScalarQuantity* q =
        new SurfaceAnimatedVertexScalarQuantity(name, permutedData, mesh, type);
    mesh.addQuantity(q);
    return q;
}

// ===========================================================================
//                      Face Vector Quantity
// ===========================================================================
SurfaceAnimatedVectorQuantity::SurfaceAnimatedVectorQuantity(
    std::string name, SurfaceMesh& mesh_, MeshElement definedOn_,
    VectorType vectorType_)
    : SurfaceMeshQuantity(name, mesh_), vectorType(vectorType_),
      vectorLengthMult(uniquePrefix() + name + "#vectorLengthMult",
                       vectorType == VectorType::AMBIENT ? absoluteValue(1.0)
                                                         : relativeValue(0.02)),
      vectorRadius(uniquePrefix() + name + "#vectorRadius",
                   relativeValue(0.0025)),
      vectorColor(uniquePrefix() + "#vectorColor", getNextUniqueColor()),
      material(uniquePrefix() + "#material", "clay"),
      time(uniquePrefix() + "#time", 0.0), definedOn(definedOn_),
      ribbonEnabled(uniquePrefix() + "#ribbonEnabled", false) {}

void SurfaceAnimatedVectorQuantity::prepareVectorMapper() {

    // Create a mapper (default mapper is identity)
    if (vectorType == VectorType::AMBIENT) {
        mapper.setMinMax(vectors);
    } else {
        mapper = AffineRemapper<glm::vec3>(vectors, DataType::MAGNITUDE);
    }
}

void SurfaceAnimatedVectorQuantity::draw() {
    if (!isEnabled()) return;

    if (program == nullptr) prepareProgram();

    // Set uniforms
    parent.setTransformUniforms(*program);

    program->setUniform("u_radius", getVectorRadius());
    program->setUniform("u_baseColor", getVectorColor());

    if (vectorType == VectorType::AMBIENT) {
        program->setUniform("u_lengthMult", 1.0);
    } else {
        program->setUniform("u_lengthMult", getVectorLengthScale());
    }

    glm::mat4 P    = view::getCameraPerspectiveMatrix();
    glm::mat4 Pinv = glm::inverse(P);
    program->setUniform("u_invProjMatrix", glm::value_ptr(Pinv));
    program->setUniform("u_viewport", render::engine->getCurrentViewport());

    program->draw();
}

void SurfaceAnimatedVectorQuantity::prepareProgram() {

    program = render::engine->generateShaderProgram(
        {render::PASSTHRU_VECTOR_VERT_SHADER, render::VECTOR_GEOM_SHADER,
         render::VECTOR_FRAG_SHADER},
        DrawMode::Points);

    // Fill buffers
    std::vector<glm::vec3> mappedVectors;
    for (glm::vec3& v : vectors) {
        mappedVectors.push_back(mapper.map(v));
    }

    program->setAttribute("a_vector", mappedVectors);
    program->setAttribute("a_position", vectorRoots);

    render::engine->setMaterial(*program, getMaterial());
}

void SurfaceAnimatedVectorQuantity::buildCustomUI() {
    ImGui::SameLine();
    if (ImGui::ColorEdit3("Color", &vectorColor.get()[0],
                          ImGuiColorEditFlags_NoInputs)) {
        setVectorColor(getVectorColor());
    }
    ImGui::SameLine();


    // === Options popup
    if (ImGui::Button("Options")) {
        ImGui::OpenPopup("OptionsPopup");
    }
    if (ImGui::BeginPopup("OptionsPopup")) {
        if (render::buildMaterialOptionsGui(material.get())) {
            material.manuallyChanged();
            setMaterial(
                material
                    .get()); // trigger the other updates that happen on set()
        }
        ImGui::EndPopup();
    }


    // Only get to set length for non-ambient vectors
    if (vectorType != VectorType::AMBIENT) {
        if (ImGui::SliderFloat("Length", vectorLengthMult.get().getValuePtr(),
                               0.0, .1, "%.5f", 3.)) {
            vectorLengthMult.manuallyChanged();
            requestRedraw();
        }
    }

    if (ImGui::SliderFloat("Radius", vectorRadius.get().getValuePtr(), 0.0, .1,
                           "%.5f", 3.)) {
        vectorRadius.manuallyChanged();
        requestRedraw();
    }

    { // Draw max and min magnitude
        ImGui::TextUnformatted(mapper.printBounds().c_str());
    }


    if (ImGui::SliderFloat("Time", &time.get(), 0.0, 1., "%.5f", 1.)) {
        time.manuallyChanged();
        recomputeVectors();
        prepareProgram();
        requestRedraw();
    }

    drawSubUI();
}

void SurfaceAnimatedVectorQuantity::drawSubUI() {}

SurfaceAnimatedVectorQuantity*
SurfaceAnimatedVectorQuantity::setVectorLengthScale(double newLength,
                                                    bool isRelative) {
    vectorLengthMult = ScaledValue<double>(newLength, isRelative);
    requestRedraw();
    return this;
}
double SurfaceAnimatedVectorQuantity::getVectorLengthScale() {
    return vectorLengthMult.get().asAbsolute();
}

SurfaceAnimatedVectorQuantity*
SurfaceAnimatedVectorQuantity::setVectorRadius(double val, bool isRelative) {
    vectorRadius = ScaledValue<double>(val, isRelative);
    requestRedraw();
    return this;
}
double SurfaceAnimatedVectorQuantity::getVectorRadius() {
    return vectorRadius.get().asAbsolute();
}

SurfaceAnimatedVectorQuantity*
SurfaceAnimatedVectorQuantity::setVectorColor(glm::vec3 color) {
    vectorColor = color;
    requestRedraw();
    return this;
}
glm::vec3 SurfaceAnimatedVectorQuantity::getVectorColor() {
    return vectorColor.get();
}

SurfaceAnimatedVectorQuantity*
SurfaceAnimatedVectorQuantity::setMaterial(std::string m) {
    material = m;
    if (program) render::engine->setMaterial(*program, getMaterial());
    if (ribbonArtist && ribbonArtist->program)
        render::engine->setMaterial(*ribbonArtist->program, material.get());
    requestRedraw();
    return this;
}
std::string SurfaceAnimatedVectorQuantity::getMaterial() {
    return material.get();
}

SurfaceAnimatedVectorQuantity*
SurfaceAnimatedVectorQuantity::setRibbonEnabled(bool val) {
    ribbonEnabled = val;
    requestRedraw();
    return this;
}
bool SurfaceAnimatedVectorQuantity::isRibbonEnabled() {
    return ribbonEnabled.get();
}


SurfaceAnimatedVectorQuantity*
SurfaceAnimatedVectorQuantity::setTime(double t) {
    time.set(t);
    requestRedraw();
    return this;
}
double SurfaceAnimatedVectorQuantity::getTime() { return time.get(); }


// ========================================================
// ==========        Intrinsic Face Vector       ==========
// ========================================================


SurfaceAnimatedFaceIntrinsicVectorQuantity::
    SurfaceAnimatedFaceIntrinsicVectorQuantity(
        std::string name, std::vector<std::vector<glm::vec2>> vectors_,
        SurfaceMesh& mesh_, int nSym_, VectorType vectorType_)
    : SurfaceAnimatedVectorQuantity(name, mesh_, MeshElement::FACE,
                                    vectorType_),
      nSym(nSym_), vectorField(vectors_) {

    // Set vector roots
    for (size_t iF = 0; iF < parent.nFaces(); iF++) {
        // Face center
        auto& face           = parent.faces[iF];
        glm::vec3 faceCenter = parent.faceCenter(iF);

        for (int iRot = 0; iRot < nSym; iRot++) {
            vectorRoots.push_back(faceCenter);
        }
    }

    // set vectors
    recomputeVectors();

    prepareVectorMapper();
}

void SurfaceAnimatedFaceIntrinsicVectorQuantity::recomputeVectors() {
    parent.ensureHaveFaceTangentSpaces();

    double rotAngle = 2.0 * PI / nSym;
    Complex rot     = std::exp(Complex(0, 1) * rotAngle);

    size_t iPrev, iNext;
    float sPrev, sNext;

    std::tie(iPrev, sPrev) = interpolate(getTime(), vectorField.size());
    iNext                  = iPrev + 1;
    sNext                  = 1 - sPrev;

    auto interpolatedValue = [&](size_t iF) {
        glm::vec2 valPrev = vectorField[iPrev][iF];
        glm::vec2 valNext = vectorField[iNext][iF];

        return sPrev * valPrev + sNext * valNext;
    };

    // Copy the vectors
    vectors.clear();
    for (size_t iF = 0; iF < parent.nFaces(); iF++) {

        glm::vec3 normal = parent.faceNormals[iF];
        glm::vec3 basisX = parent.faceTangentSpaces[iF][0];
        glm::vec3 basisY = parent.faceTangentSpaces[iF][1];

        glm::vec2 vec = interpolatedValue(iF);
        Complex angle = std::pow(Complex(vec.x, vec.y), 1.0 / nSym);

        // Face center
        auto& face           = parent.faces[iF];
        size_t D             = face.size();
        glm::vec3 faceCenter = parent.faceCenter(iF);

        for (int iRot = 0; iRot < nSym; iRot++) {
            glm::vec3 vec =
                basisX * (float)angle.real() + basisY * (float)angle.imag();
            vectors.push_back(vec);

            angle *= rot;
        }
    }
}

void SurfaceAnimatedFaceIntrinsicVectorQuantity::buildFaceInfoGUI(size_t iF) {
    ImGui::TextUnformatted(name.c_str());
    ImGui::NextColumn();

    size_t iPrev, iNext;
    float sPrev, sNext;

    std::tie(iPrev, sPrev) = interpolate(getTime(), vectorField.size());
    iNext                  = iPrev + 1;
    sNext                  = 1 - sPrev;

    auto interpolatedValue = [&](size_t iF) {
        glm::vec2 valPrev = vectorField[iPrev][iF];
        glm::vec2 valNext = vectorField[iNext][iF];

        return sPrev * valPrev + sNext * valNext;
    };

    std::stringstream buffer;
    buffer << "<" << interpolatedValue(iF).x << "," << interpolatedValue(iF).y
           << ">";
    ImGui::TextUnformatted(buffer.str().c_str());

    ImGui::NextColumn();
    ImGui::NextColumn();
    ImGui::Text("magnitude: %g", glm::length(interpolatedValue(iF)));
    ImGui::NextColumn();
}

void SurfaceAnimatedFaceIntrinsicVectorQuantity::draw() {
    SurfaceAnimatedVectorQuantity::draw();

    if (ribbonEnabled.get() && isEnabled()) {

        size_t iPrev, iNext;
        float sPrev, sNext;

        std::tie(iPrev, sPrev) = interpolate(getTime(), vectorField.size());
        iNext                  = iPrev + 1;
        sNext                  = 1 - sPrev;

        auto interpolatedValue = [&](size_t iF) {
            glm::vec2 valPrev = vectorField[iPrev][iF];
            glm::vec2 valNext = vectorField[iNext][iF];

            return sPrev * valPrev + sNext * valNext;
        };
        std::vector<glm::vec2> currentVectorField;
        for (size_t iF = 0; iF < vectorField[0].size(); ++iF) {
            currentVectorField.push_back(interpolatedValue(iF));
        }


        // Make sure we have a ribbon artist
        if (ribbonArtist == nullptr) {
            // Warning: expensive... Creates noticeable UI lag
            ribbonArtist.reset(new RibbonArtist(
                parent, traceField(parent, currentVectorField, nSym, 2500)));
            render::engine->setMaterial(*ribbonArtist->program, material.get());
        }

        // Update transform matrix from parent
        ribbonArtist->objectTransform = parent.objectTransform;
        ribbonArtist->draw();
    }
}

void SurfaceAnimatedFaceIntrinsicVectorQuantity::drawSubUI() {

    if (ImGui::Checkbox("Draw ribbon", &ribbonEnabled.get()))
        setRibbonEnabled(isRibbonEnabled());
    if (ribbonEnabled.get() && ribbonArtist != nullptr) {
        ribbonArtist->buildParametersGUI();
    }
}


std::string SurfaceAnimatedFaceIntrinsicVectorQuantity::niceName() {
    return name + " (face intrinsic vector)";
}

SurfaceAnimatedFaceIntrinsicVectorQuantity*
addAnimatedFaceIntrinsicVectorQuantityImpl(
    SurfaceMesh& mesh, std::string name,
    const std::vector<std::vector<glm::vec2>>& data, int nSym,
    VectorType type) {

    std::vector<std::vector<glm::vec2>> permutedData;

    for (const std::vector<glm::vec2>& frame : data) {
        permutedData.push_back(applyPermutation(frame, mesh.facePerm));
    }


    SurfaceAnimatedFaceIntrinsicVectorQuantity* q =
        new SurfaceAnimatedFaceIntrinsicVectorQuantity(name, permutedData, mesh,
                                                       nSym, type);
    mesh.addQuantity(q);
    return q;
}

// ========================================================
// ==========                                   ==========
// ==========          Animated Color           ==========
// ==========                                   ==========
// ========================================================

SurfaceAnimatedColorQuantity::SurfaceAnimatedColorQuantity(
    std::string name, SurfaceMesh& mesh_, std::string definedOn_)
    : SurfaceMeshQuantity(name, mesh_, true),
      time(uniquePrefix() + name + "#time", 0), definedOn(definedOn_) {}

void SurfaceAnimatedColorQuantity::draw() {
    if (!isEnabled()) return;

    if (program == nullptr) {
        createProgram();
    }

    // Set uniforms
    parent.setTransformUniforms(*program);

    program->draw();
}

void SurfaceAnimatedColorQuantity::buildCustomUI() {
    if (ImGui::SliderFloat("Time", &time.get(), 0.0, 1., "%.5f", 1.)) {
        time.manuallyChanged();
        fillColorBuffers(*program);
        requestRedraw();
    }
}

SurfaceAnimatedColorQuantity* SurfaceAnimatedColorQuantity::setTime(double t) {
    time.set(t);
    requestRedraw();
    return this;
}
double SurfaceAnimatedColorQuantity::getTime() { return time.get(); }

// ========================================================
// ==========           Vertex Color            ==========
// ========================================================

SurfaceAnimatedVertexColorQuantity::SurfaceAnimatedVertexColorQuantity(
    std::string name, std::vector<std::vector<glm::vec3>> values_,
    SurfaceMesh& mesh_)
    : SurfaceAnimatedColorQuantity(name, mesh_, "vertex"),
      values(std::move(values_)) {}

void SurfaceAnimatedVertexColorQuantity::createProgram() {
    // Create the program to draw this quantity
    program = render::engine->generateShaderProgram(
        {render::VERTCOLOR3_SURFACE_VERT_SHADER,
         render::VERTCOLOR3_SURFACE_FRAG_SHADER},
        DrawMode::Triangles);

    // Fill color buffers
    parent.fillGeometryBuffers(*program);
    fillColorBuffers(*program);
    render::engine->setMaterial(*program, parent.getMaterial());
}

void SurfaceAnimatedVertexColorQuantity::fillColorBuffers(
    render::ShaderProgram& p) {
    std::vector<glm::vec3> colorval;
    colorval.reserve(3 * parent.nFacesTriangulation());

    size_t iPrev, iNext;
    float sPrev, sNext;

    std::tie(iPrev, sPrev) = interpolate(getTime(), values.size());
    iNext                  = iPrev + 1;
    sNext                  = 1 - sPrev;

    auto interpolatedValue = [&](size_t iF) {
        glm::vec3 valPrev = values[iPrev][iF];
        glm::vec3 valNext = values[iNext][iF];

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
}

void SurfaceAnimatedVertexColorQuantity::buildVertexInfoGUI(size_t vInd) {
    ImGui::TextUnformatted(name.c_str());
    ImGui::NextColumn();

    size_t iPrev, iNext;
    float sPrev, sNext;

    std::tie(iPrev, sPrev) = interpolate(getTime(), values.size());
    iNext                  = iPrev + 1;
    sNext                  = 1 - sPrev;

    auto interpolatedValue = [&](size_t iF) {
        glm::vec3 valPrev = values[iPrev][iF];
        glm::vec3 valNext = values[iNext][iF];

        return sPrev * valPrev + sNext * valNext;
    };

    glm::vec3 tempColor = interpolatedValue(vInd);
    ImGui::ColorEdit3("", &tempColor[0],
                      ImGuiColorEditFlags_NoInputs |
                          ImGuiColorEditFlags_NoPicker);
    ImGui::SameLine();
    std::string colorStr = to_string_short(tempColor);
    ImGui::TextUnformatted(colorStr.c_str());
    ImGui::NextColumn();
}

std::string SurfaceAnimatedColorQuantity::niceName() {
    return name + " (" + definedOn + " color)";
}

void SurfaceAnimatedColorQuantity::geometryChanged() { program.reset(); }

SurfaceAnimatedVertexColorQuantity* addAnimatedVertexColorQuantityImpl(
    SurfaceMesh& mesh, std::string name,
    const std::vector<std::vector<glm::vec3>>& data) {

    std::vector<std::vector<glm::vec3>> permutedData;

    for (const std::vector<glm::vec3>& frame : data) {
        permutedData.push_back(applyPermutation(frame, mesh.vertexPerm));
    }


    SurfaceAnimatedVertexColorQuantity* q =
        new SurfaceAnimatedVertexColorQuantity(name, permutedData, mesh);
    mesh.addQuantity(q);
    return q;
}

} // namespace polyscope
