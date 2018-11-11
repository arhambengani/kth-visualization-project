/*********************************************************************
 *  Author  : Arham
 *  Init    : Thursday, December 14, 2017 - 13:08:55
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <kxcelltree/meshtocelltree.h>

namespace inviwo
{

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo KxMeshToCellTree::processorInfo_
{
    "org.inviwo.KxMeshToCellTree",      // Class identifier
    "Mesh To Cell Tree",                // Display name
    "Undefined",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};

const ProcessorInfo KxMeshToCellTree::getProcessorInfo() const
{
    return processorInfo_;
}


KxMeshToCellTree::KxMeshToCellTree()
    :Processor()
    ,portInMesh("InMesh")
    ,portOutMesh("OutMesh")
{
    addPort(portInMesh);
    addPort(portOutMesh);
    //addProperty();
}

//Function to show the triangles
void KxMeshToCellTree::ShowTriangleOutlines(const std::vector<dvec3>& InVertices, const std::vector<uvec3>& TriangleIndices,
                                            std::vector<BasicMesh::Vertex>& OutVertices, IndexBufferRAM* pIndexBufferLine)
{
    //Prepare memory
    const size_t NumVerticesBefore = OutVertices.size();
    OutVertices.reserve(OutVertices.size() + NumVerticesBefore);

    //Add all vertices
    for(const auto& InputVertex : InVertices)
    {
        OutVertices.emplace_back();
        OutVertices.back().pos = InputVertex;
    }

    //const size_t NumIndicesBefore = pIndexBufferLine->getSize();
    for(const uvec3& Triangle : TriangleIndices)
    {
        pIndexBufferLine->add((uint32_t)NumVerticesBefore + Triangle[0]); //from
        pIndexBufferLine->add((uint32_t)NumVerticesBefore + Triangle[1]); //to

        pIndexBufferLine->add((uint32_t)NumVerticesBefore + Triangle[0]); //from
        pIndexBufferLine->add((uint32_t)NumVerticesBefore + Triangle[2]); //to

        pIndexBufferLine->add((uint32_t)NumVerticesBefore + Triangle[1]); //from
        pIndexBufferLine->add((uint32_t)NumVerticesBefore + Triangle[2]); //to
    }
}
   
//Function to find the median
float KxMeshToCellTree::median(std::vector<float>& array, const size_t NumVertices)
{
    // Allocate an array of the same size and sort it.
    float* sortArray = new float[NumVertices];
    for(size_t i(0);i<NumVertices;i++) {
        sortArray[i] = array[i];
    }
    for (int i = NumVertices - 1; i > 0; --i) {
        for (int j = 0; j < i; ++j) {
            if (sortArray[j] > sortArray[j+1]) {
                float temp = sortArray[j];
                sortArray[j] = sortArray[j+1];
                sortArray[j+1] = temp;
            }
        }
    }
    
    // Middle or average of middle values in the sorted array.
    float median = 0.0;
    if ((NumVertices % 2) == 0) {
        median = (sortArray[NumVertices/2] + sortArray[(NumVertices/2) - 1])/2.0;
    } else {
        median = sortArray[NumVertices/2];
    }
    delete [] sortArray;
    return median;

}
    
//Function to draw the bounding box
void KxMeshToCellTree::DrawBoundingBox(const size_t NumVertices, const std::vector<dvec3>& Vertices, std::vector<BasicMesh::Vertex>& OutVertices, IndexBufferRAM* pIndexBufferLine)
{
    //Compute the bounding box
    
    //Iterate over all vertices: Get Bounding Box
    vec3 BoundingBoxMin(std::numeric_limits<float>::max());
    vec3 BoundingBoxMax(std::numeric_limits<float>::min());
    for(size_t i(0);i<NumVertices;i++)
    {
        //Shorthand for the currently visited vertex
        const dvec3& CurrentVertex = Vertices[i];
        
        if (BoundingBoxMax.x < CurrentVertex.x) BoundingBoxMax.x = CurrentVertex.x;
        if (BoundingBoxMax.y < CurrentVertex.y) BoundingBoxMax.y = CurrentVertex.y;
        if (BoundingBoxMax.z < CurrentVertex.z) BoundingBoxMax.z = CurrentVertex.z;
        
        if (BoundingBoxMin.x > CurrentVertex.x) BoundingBoxMin.x = CurrentVertex.x;
        if (BoundingBoxMin.y > CurrentVertex.y) BoundingBoxMin.y = CurrentVertex.y;
        if (BoundingBoxMin.z > CurrentVertex.z) BoundingBoxMin.z = CurrentVertex.z;
    }
    
//    LogInfo("Computed Bounding Box: " << BoundingBoxMin << " to " << BoundingBoxMax);
    
    //Draw the bounding box
    // - 8 vertices
    const size_t NumVerticesBefore = OutVertices.size();
    OutVertices.resize(NumVerticesBefore + 8);
    OutVertices[NumVerticesBefore + 0].pos = {BoundingBoxMin.x, BoundingBoxMin.y, BoundingBoxMin.z};
    OutVertices[NumVerticesBefore + 1].pos = {BoundingBoxMax.x, BoundingBoxMin.y, BoundingBoxMin.z};
    OutVertices[NumVerticesBefore + 2].pos = {BoundingBoxMin.x, BoundingBoxMax.y, BoundingBoxMin.z};
    OutVertices[NumVerticesBefore + 3].pos = {BoundingBoxMax.x, BoundingBoxMax.y, BoundingBoxMin.z};
    OutVertices[NumVerticesBefore + 4].pos = {BoundingBoxMin.x, BoundingBoxMin.y, BoundingBoxMax.z};
    OutVertices[NumVerticesBefore + 5].pos = {BoundingBoxMax.x, BoundingBoxMin.y, BoundingBoxMax.z};
    OutVertices[NumVerticesBefore + 6].pos = {BoundingBoxMin.x, BoundingBoxMax.y, BoundingBoxMax.z};
    OutVertices[NumVerticesBefore + 7].pos = {BoundingBoxMax.x, BoundingBoxMax.y, BoundingBoxMax.z};
    // - 12 edges
    const size_t NumIndicesBefore = pIndexBufferLine->getSize();
    pIndexBufferLine->setSize(NumIndicesBefore + 12*2);
    pIndexBufferLine->set(NumIndicesBefore +  0, (uint32_t)NumVerticesBefore + 0); //from
    pIndexBufferLine->set(NumIndicesBefore +  1, (uint32_t)NumVerticesBefore + 1); //to
    pIndexBufferLine->set(NumIndicesBefore +  2, (uint32_t)NumVerticesBefore + 0); //from
    pIndexBufferLine->set(NumIndicesBefore +  3, (uint32_t)NumVerticesBefore + 2); //to
    pIndexBufferLine->set(NumIndicesBefore +  4, (uint32_t)NumVerticesBefore + 0); //from
    pIndexBufferLine->set(NumIndicesBefore +  5, (uint32_t)NumVerticesBefore + 4); //to
    
    pIndexBufferLine->set(NumIndicesBefore +  6, (uint32_t)NumVerticesBefore + 1); //from
    pIndexBufferLine->set(NumIndicesBefore +  7, (uint32_t)NumVerticesBefore + 5); //to
    pIndexBufferLine->set(NumIndicesBefore +  8, (uint32_t)NumVerticesBefore + 1); //from
    pIndexBufferLine->set(NumIndicesBefore +  9, (uint32_t)NumVerticesBefore + 3); //to
    
    pIndexBufferLine->set(NumIndicesBefore + 10, (uint32_t)NumVerticesBefore + 2); //from
    pIndexBufferLine->set(NumIndicesBefore + 11, (uint32_t)NumVerticesBefore + 3); //to
    pIndexBufferLine->set(NumIndicesBefore + 12, (uint32_t)NumVerticesBefore + 2); //from
    pIndexBufferLine->set(NumIndicesBefore + 13, (uint32_t)NumVerticesBefore + 6); //to
    
    pIndexBufferLine->set(NumIndicesBefore + 14, (uint32_t)NumVerticesBefore + 3); //from
    pIndexBufferLine->set(NumIndicesBefore + 15, (uint32_t)NumVerticesBefore + 7); //to
    
    pIndexBufferLine->set(NumIndicesBefore + 16, (uint32_t)NumVerticesBefore + 4); //from
    pIndexBufferLine->set(NumIndicesBefore + 17, (uint32_t)NumVerticesBefore + 6); //to
    pIndexBufferLine->set(NumIndicesBefore + 18, (uint32_t)NumVerticesBefore + 4); //from
    pIndexBufferLine->set(NumIndicesBefore + 19, (uint32_t)NumVerticesBefore + 5); //to
    
    pIndexBufferLine->set(NumIndicesBefore + 20, (uint32_t)NumVerticesBefore + 5); //from
    pIndexBufferLine->set(NumIndicesBefore + 21, (uint32_t)NumVerticesBefore + 7); //to
    
    pIndexBufferLine->set(NumIndicesBefore + 22, (uint32_t)NumVerticesBefore + 6); //from
    pIndexBufferLine->set(NumIndicesBefore + 23, (uint32_t)NumVerticesBefore + 7); //to
    
}

//Function to insert into the tree
void KxMeshToCellTree::insert(int node_side, node *leaf, std::vector<dvec3> key)
{
    if(node_side == 0)
    {
        leaf->left=new node;
        leaf->left->key_value=key;
        leaf->left->median=NULL;
        leaf->left->left=NULL;    //Sets the left child of the child node to null
        leaf->left->right=NULL;
    }
    else if(node_side == 1)
    {
        leaf->right=new node;
        leaf->right->key_value=key;
        leaf->right->median=NULL;
        leaf->right->left=NULL;  //Sets the left child of the child node to null
        leaf->right->right=NULL; //Sets the right child of the child node to null
    }
}
    
void KxMeshToCellTree::makeTree_Dimension(std::vector<float> array, std::vector<dvec3> Vertices, int NumVertices, int depth, node *leaf)
{
    float medianValue = median(array, NumVertices);
    int count_left = 0, count_right = 0;
    std::vector<dvec3> left_vertices(NumVertices);
    std::vector<dvec3> right_vertices(NumVertices);
    int dimenstion = depth % 3;
    if (dimenstion == 0)
    {
        for(int j(0);j<NumVertices;j++)
        {
            const dvec3& CurrentVertex = Vertices[j];
            if (CurrentVertex.x <= medianValue) {
                left_vertices[count_left] = CurrentVertex;
                count_left++;
            }
            else if (CurrentVertex.x >= medianValue) {
                right_vertices[count_right] = CurrentVertex;
                count_right++;
            }
        }
    }
    else if (dimenstion == 1)
    {
        for(int j(0);j<NumVertices;j++)
        {
            const dvec3& CurrentVertex = Vertices[j];
            if (CurrentVertex.y <= medianValue) {
                left_vertices[count_left] = CurrentVertex;
                count_left++;
            }
            else if (CurrentVertex.y >= medianValue) {
                right_vertices[count_right] = CurrentVertex;
                count_right++;
            }
        }
    }
    else if (dimenstion == 2)
    {
        for(int j(0);j<NumVertices;j++)
        {
            const dvec3& CurrentVertex = Vertices[j];
            if (CurrentVertex.z <= medianValue) {
                left_vertices[count_left] = CurrentVertex;
                count_left++;
            }
            else if (CurrentVertex.z >= medianValue) {
                right_vertices[count_right] = CurrentVertex;
                count_right++;
            }
        }
    }
    left_vertices.resize(count_left);
    right_vertices.resize(count_right);
    leaf->median=medianValue;
    insert(0, leaf, left_vertices);
    insert(1, leaf, right_vertices);
    depth = depth + 1;
    makeTree(left_vertices, left_vertices.size(), depth, leaf->left);
    makeTree(right_vertices, right_vertices.size(), depth, leaf->right);
}
    
//Function to find the left and right node of the tree
void KxMeshToCellTree::makeTree(std::vector<dvec3> Vertices, int NumVertices, int depth, node *leaf)
{
    int dimenstion = depth % 3;
    if (dimenstion == 0 && depth < 15)
    {
        std::vector<float> array(NumVertices);
        for(int i(0);i<NumVertices;i++)
        {
            const dvec3& CurrentVertex = Vertices[i];
            array[i] = CurrentVertex.x;
        }
        makeTree_Dimension(array, Vertices, NumVertices, depth, leaf);
    }
    else if (dimenstion == 1 && depth < 15)
    {
        std::vector<float> array(NumVertices);
        for(int i(0);i<NumVertices;i++)
        {
            const dvec3& CurrentVertex = Vertices[i];
            array[i] = CurrentVertex.y;
        }
        makeTree_Dimension(array, Vertices, NumVertices, depth, leaf);
    }
    else if (dimenstion == 2 && depth < 15)
    {
        std::vector<float> array(NumVertices);
        for(int i(0);i<NumVertices;i++)
        {
            const dvec3& CurrentVertex = Vertices[i];
            array[i] = CurrentVertex.z;
        }
        makeTree_Dimension(array, Vertices, NumVertices, depth, leaf);
    }
}
    
void KxMeshToCellTree::boundingBox(std::vector<BasicMesh::Vertex>& OutVertices, IndexBufferRAM* pIndexBufferLine, node *leaf)
{
    std::vector<dvec3> np_left = leaf->left->key_value;
    int size_is_left = np_left.size();
    DrawBoundingBox(size_is_left, np_left, OutVertices, pIndexBufferLine);
    std::vector<dvec3> np_right = leaf->right->key_value;
    int size_is_right = np_right.size();
    DrawBoundingBox(size_is_right, np_right, OutVertices, pIndexBufferLine);
    if(leaf->left->left != NULL && leaf->left->right != NULL)
    {
        boundingBox(OutVertices, pIndexBufferLine, leaf->left);
    }
    if(leaf->right->left != NULL && leaf->right->right != NULL)
    {
        boundingBox(OutVertices, pIndexBufferLine, leaf->right);
    }
}
    
//Function to search the tree
node *KxMeshToCellTree::search(node *leaf, vec3 value, int depth)
{
    int dimenstion = depth % 3;
    struct node *return_value = new node;
    return_value = leaf;
    if (dimenstion == 0)
    {
        if (value.x <= leaf->median && leaf->left!=NULL) {
            return_value = search(leaf->left, value, depth+1);
        }else if (value.x > leaf->median && leaf->right!=NULL) {
            return_value = search(leaf->right, value, depth+1);
        }
    }
    else if (dimenstion == 1)
    {
        if (value.y <= leaf->median && leaf->left!=NULL) {
            return_value = search(leaf->left, value, depth+1);
        }else if (value.y > leaf->median && leaf->right!=NULL) {
            return_value = search(leaf->right, value, depth+1);
        }
    }
    else if (dimenstion == 2)
    {
        if (value.z <= leaf->median && leaf->left!=NULL) {
            return_value = search(leaf->left, value, depth+1);
        }else if (value.z > leaf->median && leaf->right!=NULL) {
            return_value = search(leaf->right, value, depth+1);
        }
    }
    
    return return_value;
}

void KxMeshToCellTree::process()
{
    //////////////////////////
    //Get input mesh
    
    std::shared_ptr<const Mesh> data = portInMesh.getData();
    if (!data) return;

    //////////////////////////
    // Get Vertices
    
    auto pit = util::find_if(data->getBuffers(), [](const auto& buf)
    {
        return buf.first.type == BufferType::PositionAttrib;
    });

    if (pit == data->getBuffers().end())
    {
        return;
    }
    const auto posRam = pit->second->getRepresentation<BufferRAM>();
    if (!posRam)
    {
        return;
    }
    if (posRam->getDataFormat()->getComponents() != 3)
    {
        return;
    }

    const size_t NumVertices = posRam->getSize();
    LogInfo("Number of vertices: " << NumVertices);

    std::vector<dvec3> Vertices(NumVertices);
    for(size_t i(0);i<NumVertices;i++)
    {
        Vertices[i] = posRam->getAsDVec3(i);
    }

    //////////////////////////
    // Get Triangles

    const auto& FirstIndexBuffer = data->getIndexBuffers()[0];
    if (FirstIndexBuffer.first.dt != DrawType::Triangles) return;
    if (FirstIndexBuffer.first.ct != ConnectivityType::None) return;
    const auto& Tris = FirstIndexBuffer.second->getRAMRepresentation()->getDataContainer();

    //Number of Triangles
    const size_t NumTriangles = Tris.size() / 3;
    LogInfo("Number of triangles: " << NumTriangles);

    std::vector<uvec3> TriangleIndices(NumTriangles);
    for (size_t i = 0; i < Tris.size(); i += 3)
    {
        TriangleIndices[i/3] = {Tris[i], Tris[i + 1], Tris[i + 2]};
    }
    
    //Get the output
    std::shared_ptr<BasicMesh> OutMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> OutVertices;
    IndexBufferRAM* pIndexBufferLine = OutMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    
    DrawBoundingBox(NumVertices, Vertices, OutVertices, pIndexBufferLine);
    
    struct node *leaf = new node;
    leaf->median=NULL;
    leaf->left=NULL;
    leaf->right=NULL;
    int depth = 0;
    makeTree(Vertices, NumVertices, depth, leaf);

    boundingBox(OutVertices, pIndexBufferLine, leaf);
    
    //We get a list of small number of candidate triangles that can be checked for interpolation point.
    std::clock_t c_start = std::clock();
    vec3 value = Vertices[7];
    struct node *new_leaf = new node;
    new_leaf = search(leaf, value, 0);
    std::clock_t c_end = std::clock();
    
    LogInfo("Computed time: " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC);
    
    OutMesh->addVertices(OutVertices);
    portOutMesh.setData(OutMesh);

//    //Draw the outlines of the triangles - just for fun
//    ShowTriangleOutlines(Vertices, TriangleIndices, OutVertices, pIndexBufferLine);

}

} // namespace

