/*********************************************************************
 *  Author  : Arham
 *  Init    : Thursday, December 14, 2017 - 13:08:55
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <kxcelltree/celltreeData.h>

namespace inviwo
{

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo KxCellTree::processorInfo_
{
    "org.inviwo.KxCellTree",      // Class identifier
    "CellTree Data Structure",    // Display name
    "Undefined",              // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};

const ProcessorInfo KxCellTree::getProcessorInfo() const
{
    return processorInfo_;
}


KxCellTree::KxCellTree()
    :Processor()
    ,portInMesh("InMesh")
    ,portOutMesh("OutMesh")
{
    addPort(portInMesh);
    addPort(portOutMesh);
    //addProperty();
}
    
//Function to insert value into the tree
void KxCellTree::insert(int node_side, celltree_node *leaf, int dim, int child, float Lmax, float Rmin, int size, int start)
{
    if(node_side == 0)
    {
        leaf->left=new celltree_node;
        leaf->left->dim = dim;
        leaf->left->child = child;
        leaf->left->node.Lmax = Lmax;
        leaf->left->node.Rmin = Rmin;
        leaf->left->leaf.start = start;
        leaf->left->leaf.size = size;
        leaf->left->left = NULL;    //Sets the left child of the child node to null
        leaf->left->right = NULL;
    }
    else if(node_side == 1)
    {
        leaf->right = new celltree_node;
        leaf->right->dim = dim;
        leaf->right->child = child;
        leaf->right->node.Lmax = Lmax;
        leaf->right->node.Rmin = Rmin;
        leaf->right->leaf.start = start;
        leaf->right->leaf.size = size;
        leaf->right->left = NULL;    //Sets the left child of the child node to null
        leaf->right->right = NULL;
    }
}

//Function to show the triangles
void KxCellTree::ShowTriangleOutlines(const std::vector<dvec3>& InVertices, const std::vector<uvec3>& TriangleIndices,
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
    
//Function to draw the bounding box
void KxCellTree::DrawBoundingBox(const size_t NumVertices, const std::vector<dvec3>& Vertices, std::vector<BasicMesh::Vertex>& OutVertices, IndexBufferRAM* pIndexBufferLine)
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
    
//Function to compute the nodes based on cell tree data structure
void KxCellTree::CellTreeDatastructure(const std::vector<dvec3>& Vertices, const size_t NumVertices, const std::vector<uvec3>& TriangleIndices, std::vector<BasicMesh::Vertex>& OutVertices, IndexBufferRAM* pIndexBufferLine, struct celltree_node *leaf, const std::vector<dvec3>& allVert)
{
    size_t size = TriangleIndices.size();
    size_t child = leaf->child;
    
    //Limiting the number of child ( can be inclreased or decreased depending on the depth of the tree desired )
    if (size > 0 && child < 15)
    {
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
    
    std::vector<dvec3> plane(5);
    for (size_t i(0.0); i<4.0; i++) {
        plane[i].x = (BoundingBoxMin.x + ((i+1.0)/5.0) * (BoundingBoxMax.x - BoundingBoxMin.x));
        plane[i].y = BoundingBoxMin.y + ((i+1.0)/5.0) * (BoundingBoxMax.y - BoundingBoxMin.y);
        plane[i].z = BoundingBoxMin.z + ((i+1.0)/5.0) * (BoundingBoxMax.z - BoundingBoxMin.z);
    }
    plane[4].x = BoundingBoxMax.x;
    plane[4].y = BoundingBoxMax.y;
    plane[4].z = BoundingBoxMax.z;
    
    std::vector<vec3> bucket(5);
    std::vector<vec3> lmaxValue(5);
    std::vector<vec3> rminValue(5);
    for (int i(0); i<5; i++) {
        lmaxValue[i].x = BoundingBoxMax.x;
        lmaxValue[i].y = BoundingBoxMax.y;
        lmaxValue[i].z = BoundingBoxMax.z;
        
        rminValue[i].x = BoundingBoxMin.x;
        rminValue[i].y = BoundingBoxMin.y;
        rminValue[i].z = BoundingBoxMin.z;
    }
    
    for(const uvec3& Triangle : TriangleIndices)
    {
        std::vector<dvec3> CurrentVertices(3);
        CurrentVertices[0] = allVert[Triangle[0]];
        CurrentVertices[1] = allVert[Triangle[1]];
        CurrentVertices[2] = allVert[Triangle[2]];
        vec3 cmax(std::numeric_limits<float>::min());
        vec3 cmin(std::numeric_limits<float>::max());
        vec3 midValue;
        for(size_t i(0);i<3;i++)
        {
            if (cmax.x < CurrentVertices[i].x) cmax.x = CurrentVertices[i].x;
            if (cmax.y < CurrentVertices[i].y) cmax.y = CurrentVertices[i].y;
            if (cmax.z < CurrentVertices[i].z) cmax.z = CurrentVertices[i].z;
            
            if (cmin.x > CurrentVertices[i].x) cmin.x = CurrentVertices[i].x;
            if (cmin.y > CurrentVertices[i].y) cmin.y = CurrentVertices[i].y;
            if (cmin.z > CurrentVertices[i].z) cmin.z = CurrentVertices[i].z;
        }
        midValue.x = (cmin.x + cmax.x)/2.0;
        midValue.y = (cmin.y + cmax.y)/2.0;
        midValue.z = (cmin.z + cmax.z)/2.0;
        if (midValue.x <= plane[0].x) {
            bucket[0].x += 1;
            if (cmax.x > lmaxValue[0].x) {
                lmaxValue[0].x = cmax.x;
            }
            if (cmin.x < rminValue[0].x) {
                rminValue[0].x = cmin.x;
            }
        } else if (midValue.x <= plane[1].x) {
            bucket[1].x += 1;
            if (cmax.x > lmaxValue[1].x) {
                lmaxValue[1].x = cmax.x;
            }
            if (cmin.x < rminValue[1].x) {
                rminValue[1].x = cmin.x;
            }
        } else if (midValue.x <= plane[2].x) {
            bucket[2].x += 1;
            if (cmax.x > lmaxValue[2].x) {
                lmaxValue[2].x = cmax.x;
            }
            if (cmin.x < rminValue[2].x) {
                rminValue[2].x = cmin.x;
            }
        } else if (midValue.x <= plane[3].x) {
            bucket[3].x += 1;
            if (cmax.x > lmaxValue[3].x) {
                lmaxValue[3].x = cmax.x;
            }
            if (cmin.x < rminValue[3].x) {
                rminValue[3].x = cmin.x;
            }
        } else if (midValue.x <= plane[4].x) {
            bucket[4].x += 1;
            if (cmax.x > lmaxValue[4].x) {
                lmaxValue[4].x = cmax.x;
            }
            if (cmin.x < rminValue[4].x) {
                rminValue[4].x = cmin.x;
            }
        }
        
        if (midValue.y <= plane[0].y) {
            bucket[0].y += 1;
            if (cmax.y > lmaxValue[0].y) {
                lmaxValue[0].y = cmax.y;
            }
            if (cmin.y< rminValue[0].y) {
                rminValue[0].y = cmin.y;
            }
        } else if (midValue.y <= plane[1].y) {
            bucket[1].y += 1;
            if (cmax.y > lmaxValue[1].y) {
                lmaxValue[1].y = cmax.y;
            }
            if (cmin.y < rminValue[1].y) {
                rminValue[1].y = cmin.y;
            }
        } else if (midValue.y <= plane[2].y) {
            bucket[2].y += 1;
            if (cmax.y > lmaxValue[2].y) {
                lmaxValue[2].y = cmax.y;
            }
            if (cmin.y < rminValue[2].y) {
                rminValue[2].y = cmin.y;
            }
        } else if (midValue.y <= plane[3].y) {
            bucket[3].y += 1;
            if (cmax.y > lmaxValue[3].y) {
                lmaxValue[3].y = cmax.y;
            }
            if (cmin.y < rminValue[3].y) {
                rminValue[3].y = cmin.y;
            }
        } else if (midValue.y <= plane[4].y) {
            bucket[4].y += 1;
            if (cmax.y > lmaxValue[4].y) {
                lmaxValue[4].y = cmax.y;
            }
            if (cmin.y < rminValue[4].y) {
                rminValue[4].y = cmin.y;
            }
        }
        
        if (midValue.z <= plane[0].z) {
            bucket[0].z += 1;
            if (cmax.z > lmaxValue[0].z) {
                lmaxValue[0].z = cmax.z;
            }
            if (cmin.z< rminValue[0].z) {
                rminValue[0].z = cmin.z;
            }
        } else if (midValue.z <= plane[1].z) {
            bucket[1].z += 1;
            if (cmax.z > lmaxValue[1].z) {
                lmaxValue[1].z = cmax.z;
            }
            if (cmin.z < rminValue[1].z) {
                rminValue[1].z = cmin.z;
            }
        } else if (midValue.z <= plane[2].z) {
            bucket[2].z += 1;
            if (cmax.z > lmaxValue[2].z) {
                lmaxValue[2].z = cmax.z;
            }
            if (cmin.z < rminValue[2].z) {
                rminValue[2].z = cmin.z;
            }
        } else if (midValue.z <= plane[3].z) {
            bucket[3].z += 1;
            if (cmax.z > lmaxValue[3].z) {
                lmaxValue[3].z = cmax.z;
            }
            if (cmin.z < rminValue[3].z) {
                rminValue[3].z = cmin.z;
            }
        } else if (midValue.z <= plane[4].z) {
            bucket[4].z += 1;
            if (cmax.z > lmaxValue[4].z) {
                lmaxValue[4].z = cmax.z;
            }
            if (cmin.z < rminValue[4].z) {
                rminValue[4].z = cmin.z;
            }
        }
    }
    
    std::vector<dvec3> optimal_plane(4);
    
    optimal_plane[0].x = (lmaxValue[0].x * (bucket[0].x)) - (rminValue[0].x * (bucket[1].x + bucket[2].x + bucket[3].x + bucket[4].x));
    optimal_plane[1].x = (lmaxValue[1].x * (bucket[0].x + bucket[1].x)) - (rminValue[1].x * (bucket[2].x + bucket[3].x + bucket[4].x));
    optimal_plane[2].x = (lmaxValue[2].x * (bucket[0].x + bucket[1].x + bucket[2].x)) - (rminValue[2].x * (bucket[3].x + bucket[4].x));
    optimal_plane[3].x = (lmaxValue[3].x * (bucket[0].x + bucket[1].x + bucket[2].x + bucket[3].x)) - (rminValue[3].x * (bucket[4].x));
    
    optimal_plane[0].y = (lmaxValue[0].y * (bucket[0].y)) - (rminValue[0].y * (bucket[1].y + bucket[2].y + bucket[3].y + bucket[4].y));
    optimal_plane[1].y = (lmaxValue[1].y * (bucket[0].y + bucket[1].y)) - (rminValue[1].y * (bucket[2].y + bucket[3].y + bucket[4].y));
    optimal_plane[2].y = (lmaxValue[2].y * (bucket[0].y + bucket[1].y + bucket[2].y)) - (rminValue[2].y * (bucket[3].y + bucket[4].y));
    optimal_plane[3].y = (lmaxValue[3].y * (bucket[0].y + bucket[1].y + bucket[2].y + bucket[3].y)) - (rminValue[3].y * (bucket[4].y));
    
    optimal_plane[0].z = (lmaxValue[0].z * (bucket[0].z)) - (rminValue[0].z * (bucket[1].z + bucket[2].z + bucket[3].z + bucket[4].z));
    optimal_plane[1].z = (lmaxValue[1].z * (bucket[0].z + bucket[1].z)) - (rminValue[1].z * (bucket[2].z + bucket[3].z + bucket[4].z));
    optimal_plane[2].z = (lmaxValue[2].z * (bucket[0].z + bucket[1].z + bucket[2].z)) - (rminValue[2].z * (bucket[3].z + bucket[4].z));
    optimal_plane[3].z = (lmaxValue[3].z * (bucket[0].z + bucket[1].z + bucket[2].z + bucket[3].z)) - (rminValue[3].z * (bucket[4].z));
    
    double optimal_plane_x = optimal_plane[0].x;
    double optimal_plane_y = optimal_plane[0].y;
    double optimal_plane_z = optimal_plane[0].z;
    int optimalPlane_x = 0, optimalPlane_y = 0, optimalPlane_z = 0;
    for (int i(0); i<4; i++) {
        if (optimal_plane_x >= optimal_plane[i].x) {
            optimal_plane_x = optimal_plane[i].x;
            optimalPlane_x = i;
        }
        if (optimal_plane_y >= optimal_plane[i].y) {
            optimal_plane_y = optimal_plane[i].y;
            optimalPlane_y = i;
        }
        if (optimal_plane_z >= optimal_plane[i].z) {
            optimal_plane_z = optimal_plane[i].z;
            optimalPlane_z = i;
        }
    }
    
    double optimal_dimension = optimal_plane_x;
    int optimalNum = optimalPlane_x;
    int dim = 0;
    if (optimal_dimension > optimal_plane_y) {
        optimal_dimension = optimal_plane_y;
        optimalNum = optimalPlane_y;
        dim = 1;
    }
    if (optimal_dimension > optimal_plane_z) {
        optimal_dimension = optimal_plane_z;
        optimalNum = optimalPlane_z;
        dim = 2;
    }
        
    //Till this part of the program we have the optimal plane or the plane along the which we will divide the node.
    
    std::vector<dvec3> vertices_left(NumVertices);
    std::vector<dvec3> vertices_right(NumVertices);
    std::vector<uvec3> left_TriangleIndices(TriangleIndices.size());
    std::vector<uvec3> right_TriangleIndices(TriangleIndices.size());
    int count_left = 0, count_right = 0;
    int countTriangle_left = 0, countTriangle_right = 0;
    if (dim == 0) {
        for (size_t i(0); i<NumVertices; i++) {
            const dvec3& CurrentVertex = Vertices[i];
            if (CurrentVertex.x <= plane[optimalNum].x) {
                vertices_left[count_left] = CurrentVertex;
                count_left++;
            }
            if (CurrentVertex.x >= plane[optimalNum].x) {
                vertices_right[count_right] = CurrentVertex;
                count_right++;
            }
        }
        int triangle_start_left = 0, triangle_start_right = 0;
        int count = 0;
        for(const uvec3& Triangle : TriangleIndices)
        {
            std::vector<dvec3> CurrentVertices(3);
            CurrentVertices[0] = allVert[Triangle[0]];
            CurrentVertices[1] = allVert[Triangle[1]];
            CurrentVertices[2] = allVert[Triangle[2]];
            vec3 cmax(std::numeric_limits<float>::min());
            vec3 cmin(std::numeric_limits<float>::max());
            for(size_t i(0);i<3;i++)
            {
                if (cmax.x < CurrentVertices[i].x) cmax.x = CurrentVertices[i].x;
                if (cmin.x > CurrentVertices[i].x) cmin.x = CurrentVertices[i].x;
            }
            if (cmin.x <= plane[optimalNum].x) {
                left_TriangleIndices[countTriangle_left] = Triangle;
                if (countTriangle_left == 0) {
                    triangle_start_left = count;
                }
                countTriangle_left++;
            }
            if (cmax.x >= plane[optimalNum].x) {
                right_TriangleIndices[countTriangle_right] = Triangle;
                if (countTriangle_right == 0) {
                    triangle_start_right = count;
                }
                countTriangle_right++;
            }
            count++;
        }
        vertices_left.resize(count_left);
        vertices_right.resize(count_right);
        left_TriangleIndices.resize(countTriangle_left);
        right_TriangleIndices.resize(countTriangle_right);
        insert(0, leaf, 0, leaf->child + 1, lmaxValue[optimalNum].x, rminValue[optimalNum].x, countTriangle_left, triangle_start_left);
        insert(1, leaf, 0, leaf->child + 1, lmaxValue[optimalNum].x, rminValue[optimalNum].x, countTriangle_right, triangle_start_right);
        DrawBoundingBox(count_left, vertices_left, OutVertices, pIndexBufferLine);
        DrawBoundingBox(count_right, vertices_right, OutVertices, pIndexBufferLine);
        CellTreeDatastructure(vertices_left, count_left, left_TriangleIndices, OutVertices, pIndexBufferLine, leaf->left, allVert);
        CellTreeDatastructure(vertices_right, count_right, right_TriangleIndices, OutVertices, pIndexBufferLine, leaf->right, allVert);
    }else if (dim == 1) {
        for (size_t i(0); i<NumVertices; i++) {
            const dvec3& CurrentVertex = Vertices[i];
            if (CurrentVertex.y <= plane[optimalNum].y) {
                vertices_left[count_left] = CurrentVertex;
                count_left++;
            }
            if (CurrentVertex.y >= plane[optimalNum].y) {
                vertices_right[count_right] = CurrentVertex;
                count_right++;
            }
        }
        int triangle_start_left = 0, triangle_start_right = 0;
        int count = 0;
        for(const uvec3& Triangle : TriangleIndices)
        {
            std::vector<dvec3> CurrentVertices(3);
            CurrentVertices[0] = allVert[Triangle[0]];
            CurrentVertices[1] = allVert[Triangle[1]];
            CurrentVertices[2] = allVert[Triangle[2]];
            vec3 cmax(std::numeric_limits<float>::min());
            vec3 cmin(std::numeric_limits<float>::max());
            for(size_t i(0);i<3;i++)
            {
                if (cmax.y < CurrentVertices[i].y) cmax.y = CurrentVertices[i].y;
                if (cmin.y > CurrentVertices[i].y) cmin.y = CurrentVertices[i].y;
            }
            if (cmin.y <= plane[optimalNum].y) {
                left_TriangleIndices[countTriangle_left] = Triangle;
                if (countTriangle_left == 0) {
                    triangle_start_left = count;
                }
                countTriangle_left++;
            }
            if (cmax.y >= plane[optimalNum].y) {
                right_TriangleIndices[countTriangle_right] = Triangle;
                if (countTriangle_right == 0) {
                    triangle_start_right = count;
                }
                countTriangle_right++;
            }
        }
        vertices_left.resize(count_left);
        vertices_right.resize(count_right);
        left_TriangleIndices.resize(countTriangle_left);
        right_TriangleIndices.resize(countTriangle_right);
        insert(0, leaf, 1, leaf->child + 1, lmaxValue[optimalNum].x, rminValue[optimalNum].x, countTriangle_left, triangle_start_left);
        insert(1, leaf, 1, leaf->child + 1, lmaxValue[optimalNum].x, rminValue[optimalNum].x, countTriangle_right, triangle_start_right);
        DrawBoundingBox(count_left, vertices_left, OutVertices, pIndexBufferLine);
        DrawBoundingBox(count_right, vertices_right, OutVertices, pIndexBufferLine);
        CellTreeDatastructure(vertices_left, count_left, left_TriangleIndices, OutVertices, pIndexBufferLine, leaf->left, allVert);
        CellTreeDatastructure(vertices_right, count_right, right_TriangleIndices, OutVertices, pIndexBufferLine, leaf->right, allVert);
    }else if (dim == 2) {
        for (size_t i(0); i<NumVertices; i++) {
            const dvec3& CurrentVertex = Vertices[i];
            if (CurrentVertex.z <= plane[optimalNum].z) {
                vertices_left[count_left] = CurrentVertex;
                count_left++;
            }
            if (CurrentVertex.z >= plane[optimalNum].z) {
                vertices_right[count_right] = CurrentVertex;
                count_right++;
            }
        }
        int triangle_start_left = 0, triangle_start_right = 0;
        int count = 0;
        for(const uvec3& Triangle : TriangleIndices)
        {
            std::vector<dvec3> CurrentVertices(3);
            CurrentVertices[0] = allVert[Triangle[0]];
            CurrentVertices[1] = allVert[Triangle[1]];
            CurrentVertices[2] = allVert[Triangle[2]];
            vec3 cmax(std::numeric_limits<float>::min());
            vec3 cmin(std::numeric_limits<float>::max());
            for(size_t i(0);i<3;i++)
            {
                if (cmax.z < CurrentVertices[i].z) cmax.z = CurrentVertices[i].z;
                if (cmin.z > CurrentVertices[i].z) cmin.z = CurrentVertices[i].z;
            }
            if (cmin.z <= plane[optimalNum].z) {
                left_TriangleIndices[countTriangle_left] = Triangle;
                if (countTriangle_left == 0) {
                    triangle_start_left = count;
                }
                countTriangle_left++;
            }
            if (cmax.z >= plane[optimalNum].z) {
                right_TriangleIndices[countTriangle_right] = Triangle;
                if (countTriangle_right == 0) {
                    triangle_start_right = count;
                }
                countTriangle_right++;
            }
        }
        vertices_left.resize(count_left);
        vertices_right.resize(count_right);
        left_TriangleIndices.resize(countTriangle_left);
        right_TriangleIndices.resize(countTriangle_right);
        insert(0, leaf, 2, leaf->child + 1, lmaxValue[optimalNum].x, rminValue[optimalNum].x, countTriangle_left, triangle_start_left);
        insert(1, leaf, 2, leaf->child + 1, lmaxValue[optimalNum].x, rminValue[optimalNum].x, countTriangle_right, triangle_start_right);
        DrawBoundingBox(count_left, vertices_left, OutVertices, pIndexBufferLine);
        DrawBoundingBox(count_right, vertices_right, OutVertices, pIndexBufferLine);
        CellTreeDatastructure(vertices_left, count_left, left_TriangleIndices, OutVertices, pIndexBufferLine, leaf->left, allVert);
        CellTreeDatastructure(vertices_right, count_right, right_TriangleIndices, OutVertices, pIndexBufferLine, leaf->right, allVert);
    }
    }
}
    
//Function to search the tree
celltree_node *KxCellTree::search(celltree_node *leaf, vec3 value)
{
    
    double val = value.x;
    if (leaf->dim == 0) {
        val = value.x;
    }else if (leaf->dim == 1) {
        val = value.y;
    }else if (leaf->dim == 2) {
        val = value.z;
    }
    struct celltree_node *return_value = new celltree_node;
    return_value = leaf;
    if (leaf->node.Lmax >= val) {
        if (leaf->left != NULL && leaf->left->leaf.size != 0) {
            return_value = search(leaf->left, value);
        }
    } else if (leaf->node.Rmin <= val) {
        if (leaf->right != NULL && leaf->right->leaf.size != 0) {
            return_value = search(leaf->right, value);
        }
    }
    return return_value;
}

void KxCellTree::process()
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
    
    float Lmax = std::numeric_limits<float>::min();
    float Rmin = std::numeric_limits<float>::max();
    for(size_t i(0);i<NumVertices;i++)
    {
        //Shorthand for the currently visited vertex
        const dvec3& CurrentVertex = Vertices[i];
        
        if (Lmax < CurrentVertex.x) Lmax = CurrentVertex.x;
        if (Rmin > CurrentVertex.x) Rmin = CurrentVertex.x;
    }
    
    struct celltree_node *leaf = new celltree_node;
    leaf->left = NULL;
    leaf->right = NULL;
    leaf->dim = 0;
    leaf->child = 0;
    leaf->node.Lmax = Lmax;
    leaf->node.Rmin = Rmin;
    leaf->leaf.start = 0;
    leaf->leaf.size = TriangleIndices.size();
    
    CellTreeDatastructure(Vertices, NumVertices, TriangleIndices, OutVertices, pIndexBufferLine, leaf, Vertices);
    
    //Draw the outlines of the triangles - just for fun
    ShowTriangleOutlines(Vertices, TriangleIndices, OutVertices, pIndexBufferLine);
    
    //We get a list of small number of candidate triangles that can be checked for interpolation point.
    std::clock_t c_start = std::clock();
    vec3 value = Vertices[7];
    struct celltree_node *node = new celltree_node;
    node = search(leaf, value);
    std::clock_t c_end = std::clock();
    
    LogInfo("Computed time: " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC);
    
    OutMesh->addVertices(OutVertices);
    portOutMesh.setData(OutMesh);
}

} // namespace
