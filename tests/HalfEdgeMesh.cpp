#ifndef IN_IDE_PARSER
#define CATCH_CONFIG_MAIN
#endif

#include "catch.hpp"

#include "HalfEdgeMesh.h"

#include <vector>
#include <iterator>
#include <unordered_map>
#include <set>
#include <random>

template< typename T >
struct InstanceCounter {
    int Count = 0;
    std::vector<const T*> instances;
    
    void add( const T *ptr ) {
        const T *b = ptr,
                *e = ptr + 1;
        for (const T *p : instances) {
            REQUIRE((e <= p || p < b));
        }
        
        instances.push_back(ptr);
    }
    void remove( const T *ptr ) {
        auto iter = std::find(instances.begin(), instances.end(), ptr);
        
        REQUIRE(iter != instances.end());
        std::swap(*iter, instances.back());
        instances.pop_back();
    }
};

// Helper to make sure that memory allocation, constructor calls and that the destructor 
// all works as intended
template< typename T >
struct InstanceCheck {
    InstanceCheck *this_;
    static InstanceCounter<T> Instances;
    
    InstanceCheck( const InstanceCheck& ) = delete;
    InstanceCheck() :
        this_(this)
    {
        Instances.Count++;
        Instances.add((T*)this);
    }
    ~InstanceCheck()
    {
        REQUIRE(this_ == this);
        this_ = nullptr;
        Instances.remove((T*)this);
        Instances.Count--;
    }
};
template< typename T >
InstanceCounter<T> InstanceCheck<T>::Instances;


struct Vertex :
    public InstanceCheck<Vertex>
{
    int data = 0;
};

struct HalfEdge  :
    public InstanceCheck<HalfEdge>
{
    int data = 0;
};

struct Edge  :
    public InstanceCheck<Edge>
{
    int data = 0;
};

struct Face  :
    public InstanceCheck<Face>
{
    int data = 0;
};

using Mesh = HalfEdgeMesh::HalfEdgeMesh<Face, Edge, HalfEdge, Vertex>;


using VertexHandle = Mesh::VertexHandle;
using HEdgeHandle = Mesh::HalfEdgeHandle;
using EdgeHandle = Mesh::EdgeHandle;
using FaceHandle = Mesh::FaceHandle;


using CVertexHandle = Mesh::CVertexHandle;
using CHEdgeHandle = Mesh::CHalfEdgeHandle;
using CEdgeHandle = Mesh::CEdgeHandle;
using CFaceHandle = Mesh::CFaceHandle;

struct EdgeInfo {
    CEdgeHandle handle;
    CVertexHandle v1, v2;
};
struct FaceInfo {
    CFaceHandle handle;
    std::vector<CVertexHandle> vertexes;
};

// This functions tries to make sure that the internal structure of the HalfEdgeMesh
// is valid and consistent
// It does NOT check that it actually contains the mesh it was given
void verifyInvariants( const Mesh &mesh,
    std::vector<CVertexHandle> vertexes,
    std::vector<CEdgeHandle> edges,
    std::vector<CFaceHandle> faces )
{
    std::vector<CHEdgeHandle> hedges;
    
    auto hasHEdge = [&]( CHEdgeHandle hedge ) {
        for (CHEdgeHandle e : hedges) {
            if (e == hedge) return true;
        }
        return false;
    };
    
    auto hasVertex = [&]( CVertexHandle vertex ) {
        for (CVertexHandle v : vertexes) {
            if (v == vertex) return true;
        }
        return false;
    };
    
    auto hasFace = [&]( CFaceHandle face ) {
        for (CFaceHandle f : faces) {
            if (f == face) return true;
        }
        return false;
    };
    
    
    // Verify that all half edge pointers are valid
    for (CEdgeHandle edge : edges) {
        INFO("Edge " << edge);
        CHEdgeHandle hedge = mesh.getEdgeHalfEdge(edge);
        REQUIRE(hedge);
        
        CHEdgeHandle pair = mesh.getHalfEdgePair(hedge);
        REQUIRE(pair);
        
        REQUIRE(mesh.getHalfEdgePair(pair) == hedge);
        REQUIRE(mesh.getHalfEdgeEdge(hedge) == edge);
        REQUIRE(mesh.getHalfEdgeEdge(pair) == edge);
        
        CVertexHandle v1 = mesh.getHalfEdgeVertex(hedge),
                     v2 = mesh.getHalfEdgeVertex(pair);
                     
        REQUIRE(v1);
        REQUIRE(v2);
        
        REQUIRE(hasVertex(v1));
        REQUIRE(hasVertex(v2));
        
        
        CFaceHandle f1 = mesh.getHalfEdgeFace(hedge),
                   f2 = mesh.getHalfEdgeFace(pair);
                   
        if (f1) {
            REQUIRE(hasFace(f1));
        }
        if (f2) {
            REQUIRE(hasFace(f2));
        }
        
        
        hedges.push_back(hedge);
        hedges.push_back(pair);
    }
    
    for (CHEdgeHandle hedge : hedges) {
        INFO("Half Edge " << hedge);
        CHEdgeHandle next = mesh.getHalfEdgeNext(hedge);
        REQUIRE(next);
        REQUIRE(hasHEdge(next));
        
        REQUIRE(mesh.getHalfEdgePrev(next) == hedge);   
    }
    
    // Verify that you can walk around the vertex
    for (CVertexHandle vertex : vertexes ) {
        INFO("Vertex " << vertex);
        CHEdgeHandle hedge = mesh.getVertexHalfEdge(vertex);
        
        if (hedge == nullptr) continue;
        
        // Walk around the vertex
        CHEdgeHandle root = hedge;
        
        std::set<CHEdgeHandle> visited;
        do {
            // Make sure that we only visit each half edge once
            REQUIRE(visited.count(hedge) == 0);
            visited.insert(hedge);
            
            REQUIRE(hedge);
            REQUIRE(hasHEdge(hedge));
            
            CHEdgeHandle pair = mesh.getHalfEdgePair(hedge);
            REQUIRE(pair);
            
            REQUIRE(mesh.getHalfEdgeVertex(pair) == vertex);
            
            hedge = mesh.getHalfEdgeNext(pair);
        } while (hedge != root);
    }
    
    // Check that all half edges can be reached from its vertex
    for (CHEdgeHandle hedge : hedges) {
        CHEdgeHandle pair = mesh.getHalfEdgePair(hedge);
        CVertexHandle vertex = mesh.getHalfEdgeVertex(pair);
        
        CHEdgeHandle root = mesh.getVertexHalfEdge(vertex);
        CHEdgeHandle e = root;
        do {
            if (e == hedge) break;
            
            e = mesh.getHalfEdgePair(e);
            e = mesh.getHalfEdgeNext(e);
        } while (root != e);
        REQUIRE(e == hedge);
    }
    
    // Verify that you can walk around the faces
    for (CFaceHandle face : faces) {
        CHEdgeHandle hedge = mesh.getFaceHalfEdge(face);
        REQUIRE(hedge);
        
        CHEdgeHandle root = hedge;
        std::set<CHEdgeHandle> visited;
        do {
            // Make sure that we only visit each half edge once
            REQUIRE(visited.count(hedge) == 0);
            visited.insert(hedge);
            
            REQUIRE(mesh.getHalfEdgeFace(hedge) == face);
            hedge = mesh.getHalfEdgeNext(hedge);
        } while (hedge != root);
    }
    
    // Verify that all half edges can be reaches from its face
    for (CHEdgeHandle hedge : hedges) {
        CFaceHandle face = mesh.getHalfEdgeFace(hedge);
        if (face == nullptr) continue;
        
        CHEdgeHandle root = mesh.getFaceHalfEdge(face);
        CHEdgeHandle e = root;
        do {
            if (e == hedge) break;
            
            e = mesh.getHalfEdgeNext(e);
        } while (root != e);
        REQUIRE(e == hedge);
    }
}

void dumpState( const Mesh &mesh,
    std::vector<CVertexHandle> vertexes,
    std::vector<CEdgeHandle> edges,
    std::vector<CFaceHandle> faces )
{
    std::vector<CHEdgeHandle> hedges;    
    
    for (CEdgeHandle edge : edges) {
        CHEdgeHandle hedge = mesh.getEdgeHalfEdge(edge);
        CHEdgeHandle pair = mesh.getHalfEdgePair(hedge);
        
        hedges.push_back(hedge);
        hedges.push_back(pair);
    }
    
    std::unordered_map<CVertexHandle, int> vertexLookup;
    std::unordered_map<CEdgeHandle, int> edgeLookup;
    std::unordered_map<CHEdgeHandle, int> hedgeLookup;
    std::unordered_map<CFaceHandle, int> faceLookup;
    
    vertexLookup.emplace(nullptr, -1);
    edgeLookup.emplace(nullptr, -1);
    hedgeLookup.emplace(nullptr, -1);
    faceLookup.emplace(nullptr, -1);
    
    {
        printf("=== Vertex Translation Table ===\n");
        int idx = 0;
        for (CVertexHandle v : vertexes) {
            printf("%p\t%i\n", v, idx);
            vertexLookup.emplace(v, idx++);
        }
        idx = 0;
        printf("=== Edge Translation Table ===\n");
        for (CEdgeHandle e : edges) {
            printf("%p\t%i\n", e, idx);
            edgeLookup.emplace(e, idx++);
        }
        idx = 0;
        printf("=== Half Edge Translation Table ===\n");
        for (CHEdgeHandle e : hedges) {
            printf("%p\t%i\n", e, idx);
            hedgeLookup.emplace(e, idx++);
        }
        idx = 0;
        printf("=== Face Translation Table ===\n");
        for (CFaceHandle f : faces) {
            printf("%p\t%i\n", f, idx);
            faceLookup.emplace(f, idx++);
        }
    }
    
    
    printf("************************\n");
    for (CVertexHandle v : vertexes) {
        CHEdgeHandle e = mesh.getVertexHalfEdge(v);
        
        printf("Vertex %i:\n"
               "   hedge: %i\n",
               vertexLookup[v],
               hedgeLookup[e]
        );
    }
    
    printf("************************\n");
    for (CEdgeHandle e : edges) {
        CHEdgeHandle h = mesh.getEdgeHalfEdge(e);
        printf("Edge %i:\n"
               "   hedge: %i\n",
               edgeLookup[e],
               hedgeLookup[h]
        );
    }
    
    printf("************************\n");
    for (CFaceHandle f : faces) {
        CHEdgeHandle h = mesh.getFaceHalfEdge(f);
        printf("Face %i:\n"
               "   hedge: %i\n",
               faceLookup[f],
               hedgeLookup[h]
        );
    }
    
    printf("************************\n");
    for (CHEdgeHandle e : hedges) {
        CHEdgeHandle pair = mesh.getHalfEdgePair(e), 
                    next = mesh.getHalfEdgeNext(e), 
                    prev = mesh.getHalfEdgePrev(e);
        CVertexHandle vert = mesh.getHalfEdgeVertex(e);
        CEdgeHandle edge = mesh.getHalfEdgeEdge(e);
        CFaceHandle face = mesh.getHalfEdgeFace(e);
        
        printf("HalfEdge %i:\n"
               "   pair: %i\n"
               "   next: %i\n"
               "   prev: %i\n"
               "   vert: %i\n"
               "   edge: %i\n"
               "   face: %i\n",
               hedgeLookup[e],
               hedgeLookup[pair],
               hedgeLookup[next],
               hedgeLookup[prev],
               vertexLookup[vert],
               edgeLookup[edge],
               faceLookup[face]
        );
    }
    
    
    
}

// Verifies that the internal structure of the HalfEdgeMesh reprecents the expected mesh
void verifyMesh( const Mesh &mesh, 
    std::vector<CVertexHandle> vertexes,
    std::vector<EdgeInfo> edges,
    std::vector<FaceInfo> faces
)
{
    for (EdgeInfo &edge : edges) {
        auto handle = mesh.findEdge(edge.v1, edge.v2);
        if (!edge.handle) {
            REQUIRE(handle);
            edge.handle = handle;
        }
        else {
            REQUIRE(handle == edge.handle);
        }
    }

    for (const FaceInfo &face : faces) {
        CHEdgeHandle edge = mesh.findHalfEdge(
            face.vertexes[0],   
            face.vertexes[1]
        );
        REQUIRE(edge);

        for (int i=1; i <= face.vertexes.size(); ++i) {
            REQUIRE(mesh.getHalfEdgeFace(edge) == face.handle);
            CVertexHandle vertex = face.vertexes[i % face.vertexes.size()];

            CVertexHandle edgeVert = mesh.getHalfEdgeVertex(edge);
            REQUIRE(edgeVert == vertex);

            edge = mesh.getHalfEdgeNext(edge);
            REQUIRE(edge);
        }
    }
    SECTION("Iterators")
    {
        { // Vertexes
            std::set<CVertexHandle> visited;

            for (CVertexHandle handle : mesh.vertexes()) {
                REQUIRE(visited.count(handle) == 0);
                visited.insert(handle);
            }

            REQUIRE(visited.size() == vertexes.size());
            for (CVertexHandle handle : vertexes) {
                REQUIRE(visited.count(handle) == 1);
            }
        }
        { // HEdges
            std::set<CHEdgeHandle> visited;

            for (CHEdgeHandle handle : mesh.halfEdges()) {
                REQUIRE(visited.count(handle) == 0);
                visited.insert(handle);
            }

            REQUIRE(visited.size() == edges.size()*2);
            for (const auto &edge : edges) {
                CHEdgeHandle hedge = mesh.getEdgeHalfEdge(edge.handle);
                CHEdgeHandle pair = mesh.getHalfEdgePair(hedge);

                REQUIRE(visited.count(hedge) == 1);
                REQUIRE(visited.count(pair) == 1);
            }
        }
        { // Edges
            std::set<CEdgeHandle> visited;

            for (CEdgeHandle handle : mesh.edges()) {
                REQUIRE(visited.count(handle) == 0);
                visited.insert(handle);
            }

            REQUIRE(visited.size() == edges.size());
            for (const auto &edge : edges) {
                REQUIRE(visited.count(edge.handle) == 1);
            }
        }
        { // Faces
            std::set<CFaceHandle> visited;

            for (CFaceHandle handle : mesh.faces()) {
                REQUIRE(visited.count(handle) == 0);
                visited.insert(handle);
            }

            REQUIRE(visited.size() == faces.size());
            for (const auto &face : faces) {
                REQUIRE(visited.count(face.handle) == 1);
            }
        }
        { // VertexVertex
            for (CVertexHandle vertex : vertexes) {
                std::vector<CVertexHandle> handles;

                for (const auto &edge : edges) {
                    if (edge.v1 == vertex) {
                        handles.push_back(edge.v2);
                    }
                    if (edge.v2 == vertex) {
                        handles.push_back(edge.v1);
                    }
                }

                std::set<CVertexHandle> visited;
                for (CVertexHandle handle : mesh.vertexVertexes(vertex)) {
                    REQUIRE(visited.count(handle) == 0);
                    visited.insert(handle);
                }
                
                REQUIRE(visited.size() == handles.size());
                for (const auto &vertex : handles) {
                    REQUIRE(visited.count(vertex) == 1);
                }
            }
        }
        { // VertexOutHalfEdge
            for (CVertexHandle vertex : vertexes) {
                std::vector<CHEdgeHandle> handles;

                for (const auto &edge : edges) {
                    if (edge.v1 != vertex && edge.v2 != vertex) continue;

                    CHEdgeHandle hedge = mesh.getEdgeHalfEdge(edge.handle);
                    if (mesh.getHalfEdgeVertex(hedge) != vertex) {
                        handles.push_back(hedge);
                    }
                    else {
                        handles.push_back(mesh.getHalfEdgePair(hedge));
                    }
                }

                std::set<CHEdgeHandle> visited;
                for (auto handle : mesh.vertexOutHalfEdges(vertex)) {
                    REQUIRE(visited.count(handle) == 0);
                    visited.insert(handle);
                }
                
                REQUIRE(visited.size() == handles.size());
                for (const auto &handle : handles) {
                    REQUIRE(visited.count(handle) == 1);
                }
            }
        }
        { // VertexInHalfEdge
            for (CVertexHandle vertex : vertexes) {
                std::vector<CHEdgeHandle> handles;

                for (const auto &edge : edges) {
                    if (edge.v1 != vertex && edge.v2 != vertex) continue;

                    CHEdgeHandle hedge = mesh.getEdgeHalfEdge(edge.handle);
                    if (mesh.getHalfEdgeVertex(hedge) == vertex) {
                        handles.push_back(hedge);
                    }
                    else {
                        handles.push_back(mesh.getHalfEdgePair(hedge));
                    }
                }

                std::set<CHEdgeHandle> visited;
                for (auto handle : mesh.vertexInHalfEdges(vertex)) {
                    REQUIRE(visited.count(handle) == 0);
                    visited.insert(handle);
                }
                
                REQUIRE(visited.size() == handles.size());
                for (const auto &handle : handles) {
                    REQUIRE(visited.count(handle) == 1);
                }
            }
        }
        { // VertexEdge
            for (CVertexHandle vertex : vertexes) {
                std::vector<CEdgeHandle> handles;

                for (const auto &edge : edges) {
                    if (edge.v1 == vertex || edge.v2 == vertex) {
                        handles.push_back(edge.handle);
                    }
                }

                std::set<CEdgeHandle> visited;
                for (auto handle : mesh.vertexEdges(vertex)) {
                    REQUIRE(visited.count(handle) == 0);
                    visited.insert(handle);
                }
                
                REQUIRE(visited.size() == handles.size());
                for (const auto &handle : handles) {
                    REQUIRE(visited.count(handle) == 1);
                }
            }
        }
        { // VertexFace
            for (CVertexHandle vertex : vertexes) {
                std::vector<CFaceHandle> handles;

                for (const auto &face : faces) {
                    for (CVertexHandle handle : face.vertexes) {
                        if (vertex == handle) {
                            handles.push_back(face.handle);
                            break;
                        }
                    }
                }

                std::set<CFaceHandle> visited;
                for (auto handle : mesh.vertexFaces(vertex)) {
                    REQUIRE(visited.count(handle) == 0);
                    visited.insert(handle);
                }
                
                REQUIRE(visited.size() == handles.size());
                for (const auto &handle : handles) {
                    REQUIRE(visited.count(handle) == 1);
                }
            }
        }
        { // FaceVertex
            for (const FaceInfo &info : faces) {
                std::set<CVertexHandle> visited;
                for (auto handle : mesh.faceVertexes(info.handle)) {
                    REQUIRE(visited.count(handle) == 0);
                    visited.insert(handle);
                }

                REQUIRE(visited.size() == info.vertexes.size());
                for (const auto &handle : info.vertexes) {
                    REQUIRE(visited.count(handle) == 1);
                }
            }
        }
        { // FaceHEdge
            for (const FaceInfo &info : faces) {
                std::set<CHEdgeHandle> visited;
                for (auto handle : mesh.faceHalfEdges(info.handle)) {
                    REQUIRE(visited.count(handle) == 0);
                    visited.insert(handle);
                }

                REQUIRE(visited.size() == info.vertexes.size());
                for (size_t i=0, s=info.vertexes.size(); i < s; ++i) {
                    CHEdgeHandle handle = mesh.findHalfEdge(info.vertexes[i], info.vertexes[(i+1)%s]);
                    REQUIRE (handle);
                    REQUIRE(visited.count(handle) == 1);
                }
            }
        }
        { // FaceEdge
            for (const FaceInfo &info : faces) {
                std::set<CEdgeHandle> visited;
                for (auto handle : mesh.faceEdges(info.handle)) {
                    REQUIRE(visited.count(handle) == 0);
                    visited.insert(handle);
                }

                REQUIRE(visited.size() == info.vertexes.size());
                for (size_t i=0, s=info.vertexes.size(); i < s; ++i) {
                    CEdgeHandle handle = mesh.findEdge(info.vertexes[i], info.vertexes[(i+1)%s]);
                    REQUIRE (handle);
                    REQUIRE(visited.count(handle) == 1);
                }
            }
        }
        { // FaceFace 
            for (const FaceInfo &info : faces) {
                std::set<CFaceHandle> neighbors;
                
                for (size_t i=0, s=info.vertexes.size(); i < s; ++i) {
                    CVertexHandle v1 = info.vertexes[i],
                                 v2 = info.vertexes[(i+1)%s];

                    for (const FaceInfo &face : faces) {
                        if (face.handle == info.handle) continue;

                        for (size_t a=0,b=face.vertexes.size(); a < b; ++a) {
                            if (face.vertexes[a] == v2 && face.vertexes[(a+1)%b] == v1) {
                                REQUIRE(neighbors.count(face.handle) == 0);
                                neighbors.insert(face.handle);
                            }
                        }
                    }
                }


                std::set<CFaceHandle> visited;
                for (CFaceHandle handle : mesh.faceFaces(info.handle)) {
                    REQUIRE(handle);
                    REQUIRE(visited.count(handle) == 0);
                    visited.insert(handle);
                }
                
                REQUIRE(visited.size() == neighbors.size());
                for (auto handle : neighbors) {
                    REQUIRE(visited.count(handle) == 1);
                }
            }
        }
    }
}

TEST_CASE( "HalfEdgeMesh", "[Common][HalfEdgeMesh]" )
{
    SECTION("SlabAllocator")
    {
        using SlabAllocator = HalfEdgeMesh::internal::SlabAllocator<Vertex, 8>;
        
        
        SlabAllocator allocator;
        
        std::vector<Vertex*> vertexes;
        
        const int Target = 128,
                  Dev = 32;
        
        std::ranlux24 randEng;
        using Dist = std::bernoulli_distribution;
        Dist rand;
        for (int i=0; i < 4048; ++i) {
            float p = (Target - vertexes.size()) / (float)Dev;
            if (p < 0.f) p = 0.f;
            if (p > 1.f) p = 1.f;
            
            if (rand(randEng,Dist::param_type(p))) {
                int c = Vertex::Instances.Count;
                Vertex *v = allocator.allocate();
                REQUIRE((c+1) == Vertex::Instances.Count);
                
                REQUIRE(std::find(vertexes.begin(), vertexes.end(),v) == vertexes.end());
                vertexes.push_back(v);
            }
            else {
                int idx = std::uniform_int_distribution<>(0, vertexes.size()-1)(randEng);
                Vertex *v = vertexes[idx];
                
                int c = Vertex::Instances.Count;
                allocator.free(v);
                REQUIRE((c-1) == Vertex::Instances.Count);
                
                vertexes.erase(vertexes.begin()+idx);
            }
        }
        
        // Now the allocator should have a good mix of used and free slabs over multiple blocks
        auto beg = allocator.begin(),
             end = allocator.end();
             
        REQUIRE(beg != end);
        
        for (auto iter=beg; iter != end; ++iter) {
            auto pos = std::find(vertexes.begin(), vertexes.end(), &*iter);
            REQUIRE(pos != vertexes.end());
            
            REQUIRE(iter->data == 0);
            iter->data = 1;
        }
        
        // Make sure all vertexes was visited
        for (auto v : vertexes) {
            REQUIRE(v->data == 1);
        }
        
        // in reverse
        auto iter = end;
        do {
            --iter;
            
            REQUIRE(iter->data == 1);
            iter->data = 2;
            
        } while(iter != beg);
        
        // Make sure all vertexes was visited
        for (auto v : vertexes) {
            REQUIRE(v->data == 2);
        }
    }
    
    SECTION("Edges") 
    {
        Mesh mesh;
        /*
            1 -- 2
            |  /
            | /
            3    4
        */

        auto v1 = mesh.createVertex(),
             v2 = mesh.createVertex(),
             v3 = mesh.createVertex(),
             v4 = mesh.createVertex();

        auto e1 = mesh.createEdge(v1, v2),
             e2 = mesh.createEdge(v2, v3),
             e3 = mesh.createEdge(v3, v1);

        REQUIRE(e1);
        REQUIRE(e2);
        REQUIRE(e3);

        
        verifyInvariants(mesh,
            {v1,v2,v3,v4},
            {e1,e2,e3},
            {}
        );
        verifyMesh(mesh,
            {v1,v2,v3,v4},
            {{e1,v1,v2}, {e2,v2,v3}, {e3,v3,v1}},
            {}
        );
    }

    SECTION("Faces")
    {
        Mesh mesh;
        /*
                  3      4
                  | \  /   \ <=f2
              f1=>|  1------5
                  |/  \     
                  2    6
        */

        auto v1 = mesh.createVertex(),
             v2 = mesh.createVertex(),
             v3 = mesh.createVertex(),
             v4 = mesh.createVertex(),
             v5 = mesh.createVertex(),
             v6 = mesh.createVertex();

        REQUIRE(v1);
        REQUIRE(v2);
        REQUIRE(v3);
        REQUIRE(v4);
        REQUIRE(v5);
        REQUIRE(v6);

        auto e1 = mesh.createEdge(v1, v2),
             e2 = mesh.createEdge(v1, v6);

        REQUIRE(e1);
        REQUIRE(e2);
        
        auto f1Vertexes = {v1, v2, v3};
        auto f1 = mesh.createFace(f1Vertexes.begin(), f1Vertexes.size());
        REQUIRE(f1);
        
        auto e3 = mesh.findEdge(v2,v3),
             e4 = mesh.findEdge(v1,v3);
             
        REQUIRE(e3);
        REQUIRE(e4);

        verifyInvariants(mesh,
            {v1,v2,v3,v4,v5,v6},
            {e1,e2,e3,e4},
            {f1}
        );
        
        auto f2Vertexes = {v1, v4, v5};
        auto f2 = mesh.createFace(f2Vertexes.begin(), f2Vertexes.size());
        REQUIRE(f2);

        auto e5 = mesh.findEdge(v1,v4),
             e6 = mesh.findEdge(v4,v5),
             e7 = mesh.findEdge(v5,v1);
             
        REQUIRE(e5);
        REQUIRE(e6);
        REQUIRE(e7);
        
        verifyInvariants(mesh,
            {v1,v2,v3,v4,v5,v6},
            {e1,e2,e3,e4,e5,e6,e7},
            {f1,f2}
        );
        
        /*
                  3   6  4
                  | \ | /   \ <=f2
              f1=>|  1------5
                  |/    ___/ <=f3
                  2----/
        */

        auto f3Vertexes = {v2,v1,v5};
        auto f3 = mesh.createFace(f3Vertexes.begin(), f3Vertexes.size());

        auto e8 = mesh.findEdge(v2, v5);
        
        auto f4Vertexes = {v1,v3,v4};
        auto f4 = mesh.createFace(f4Vertexes.begin(), f4Vertexes.size());
        
        REQUIRE(f4 == nullptr);
        
        auto e9 = mesh.findEdge(v3,v4);
        
        verifyInvariants(mesh,
            {v1,v2,v3,v4,v5,v6},
            {e1,e2,e3,e4,e5,e6,e7,e8,e9},
            {f1,f2,f3}
        );
        verifyMesh(mesh,
            {v1, v2, v3, v4, v5, v6},
            {{e1,v1, v2}, {{}, v2, v3}, {{}, v3,v1}, {e2, v1, v6}, {{}, v1, v4}, {{}, v4, v5}, {{}, v1, v5}, {{}, v5, v2}, {{}, v3, v4}},
            {{f1, {v1, v2, v3}}, {f2, {v1, v4, v5}}, {f3, {v2, v1, v5}}}
        );
    }
    
    REQUIRE(Vertex::Instances.Count == 0);
    REQUIRE(HalfEdge::Instances.Count == 0);
    REQUIRE(Edge::Instances.Count == 0);
    REQUIRE(Face::Instances.Count == 0);
}
