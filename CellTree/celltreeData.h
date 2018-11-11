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
#include <math.h>
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

const int Nb = 5;
    
struct celltree_node
{
    unsigned int dim;
    unsigned int child;
    
    union {
        struct { float Lmax, Rmin; } node;
        struct { unsigned int start, size; } leaf;
    };
    
    celltree_node *left;
    celltree_node *right;
};

class IVW_MODULE_KXCELLTREE_API KxCellTree : public Processor
{
//Friends
//Types
public:

//Construction / Deconstruction
public:
    KxCellTree();
    virtual ~KxCellTree() = default;

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

    void DrawBoundingBox(const size_t NumVertices, const std::vector<dvec3>& Vertices,
                         std::vector<BasicMesh::Vertex>& OutVertices,
                         IndexBufferRAM* pIndexBufferLine);
    
    void CellTreeDatastructure(const std::vector<dvec3>& Vertices, const size_t NumVertices, const std::vector<uvec3>& TriangleIndices, std::vector<BasicMesh::Vertex>& OutVertices, IndexBufferRAM* pIndexBufferLine, struct celltree_node *leaf, const std::vector<dvec3>& allVert);
    
    void insert(int node_side, celltree_node *leaf, int dim, int child, float Lmax, float Rmin, int size, int start);
    
    celltree_node *search(celltree_node *leaf, vec3 value);

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
