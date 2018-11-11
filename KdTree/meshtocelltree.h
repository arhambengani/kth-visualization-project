/*********************************************************************
 *  Author  : Arham
 *  Init    : Thursday, December 14, 2017 - 13:08:55
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <kxcelltree/kxcelltreemoduledefine.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/ports/meshport.h>
#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <ctime>
//#include <inviwo/core/ports/volumeport.h>
//#include <inviwo/core/properties/boolcompositeproperty.h>
//#include <inviwo/core/properties/boolproperty.h>
//#include <inviwo/core/properties/buttonproperty.h>
//#include <inviwo/core/properties/compositeproperty.h>
//#include <inviwo/core/properties/fileproperty.h>
//#include <inviwo/core/properties/minmaxproperty.h>
//#include <inviwo/core/properties/optionproperty.h>
//#include <inviwo/core/properties/ordinalproperty.h>
//#include <inviwo/core/properties/stringproperty.h>
//#include <inviwo/core/properties/transferfunctionproperty.h>

namespace inviwo
{

/** \docpage{org.inviwo.KxMeshToCellTree, Mesh To Cell Tree}
    ![](org.inviwo.KxMeshToCellTree.png?classIdentifier=org.inviwo.KxMeshToCellTree)

    Explanation of how to use the processor.
    
    ### Inports
      * __<Inport1>__ <description>.
    
    ### Outports
      * __<Outport1>__ <description>.
    
    ### Properties
      * __<Prop1>__ <description>.
      * __<Prop2>__ <description>
*/


/** \class KxMeshToCellTree
    \brief VERY_BRIEFLY_DESCRIBE_THE_PROCESSOR
    
    DESCRIBE_THE_PROCESSOR_FROM_A_DEVELOPER_PERSPECTIVE

    @author Arham
*/
    
struct node
{
    std::vector<dvec3> key_value;
    float median;
    node *left;
    node *right;
};
    
class IVW_MODULE_KXCELLTREE_API KxMeshToCellTree : public Processor
{ 
//Friends
//Types
public:

//Construction / Deconstruction
public:
    KxMeshToCellTree();
    virtual ~KxMeshToCellTree() = default;

//Methods
public:
    virtual const ProcessorInfo getProcessorInfo() const override;
    static const ProcessorInfo processorInfo_;

protected:
    ///Our main computation function
    virtual void process() override;

    void ShowTriangleOutlines(const std::vector<dvec3>& Vertices, const std::vector<uvec3>& TriangleIndices,
                              std::vector<BasicMesh::Vertex>& OutVertices,
                              IndexBufferRAM* pIndexBufferLine);
    
    void makeTree(std::vector<dvec3> Vertices, int NumVertices, int depth, node *leaf);
    
    void DrawBoundingBox(const size_t NumVertices, const std::vector<dvec3>& Vertices,
                         std::vector<BasicMesh::Vertex>& OutVertices,
                         IndexBufferRAM* pIndexBufferLine);
    
    float median(std::vector<float>& array, const size_t NumVertices);
    
    void insert(int node_side, node *leaf, std::vector<dvec3> key);
    
    void makeTree_Dimension(std::vector<float> array, std::vector<dvec3> Vertices, int NumVertices, int depth, node *leaf);
    
    void boundingBox(std::vector<BasicMesh::Vertex>& OutVertices, IndexBufferRAM* pIndexBufferLine, node *leaf);
    
    node *search(node *leaf, vec3 value, int depth);

//Ports
public:
    ///Input of the triangle mesh
    MeshInport portInMesh;

    ///Output of the outline of the kD-tree or cell tree
    MeshOutport portOutMesh;

//Properties
public:

//Attributes
private:

};

} // namespace
