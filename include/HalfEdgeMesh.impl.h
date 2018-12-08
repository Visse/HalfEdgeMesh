#pragma once

#ifndef HALF_EDGE_MESH_H
// mostly for the ide
#include "HalfEdgeMesh.h"
#endif

namespace HalfEdgeMesh
{
# define HALF_EDGE_MESH_T template< typename FaceData, typename EdgeData, typename HalfEdgeData, typename VertexData >
# define HALF_EDGE_MESH_C HalfEdgeMesh<FaceData, EdgeData, HalfEdgeData, VertexData>

    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::VertexHandle 
    HALF_EDGE_MESH_C::createVertex()
    {
        VertexEntry *v = mVertexes.allocate();
        if (v == nullptr) {
            HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createVertex - Failed to allocate vertex.\n");
        }
        return v;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::EdgeHandle 
    HALF_EDGE_MESH_C::createEdge( CVertexHandle v1_, CVertexHandle v2_ )
    {
        VertexEntry *v1 = asEntry(v1_),
                    *v2 = asEntry(v2_);
        
        if (v1 == nullptr || v2 == nullptr || v1 == v2) {
            HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createEdge - Either v1 or v2 was null.\n");
            return nullptr;
        }
        
        if (EdgeHandle edge = findEdge(v1, v2)) {
            return edge;
        }
        
        HalfEdgeEntry *v1_insertion_point = findFreeHEdge(v1->hedge),
                      *v2_insertion_point = findFreeHEdge(v2->hedge);
        
        // If we fail to find a free edge, that means creating a edge would lead to a
        // non-manifold mesh
        if ((v1->hedge && v1_insertion_point == nullptr) ||
            (v2->hedge && v2_insertion_point == nullptr)) 
        {
            HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createEdge - Failed to find a insertion point.\n");
            return nullptr;
        }
        
        EdgeEntry *edge = mEdges.allocate();
        if (!edge) {
            HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createEdge - Failed to allocate edge.\n");
            return nullptr;
        }
        
        HalfEdgeEntry *h1 = mHalfEdges.allocate(),
                      *h2 = mHalfEdges.allocate();
        
        if (h1 == nullptr || h2 == nullptr) {
            mEdges.free(edge);
            mHalfEdges.free(h1);
            mHalfEdges.free(h2);
            HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createEdge - Failed to allocate half edges.\n");
            return nullptr;
        }
        
        edge->hedge = h1;
        
        { // init h1
            h1->edge = edge;
            h1->vertex = v2;
            h1->pair = h2;
            h1->next = h2;
            h1->prev = h2;
        }
        { // init h2
            h2->edge = edge;
            h2->vertex = v1;
            h2->pair = h1;
            h2->next = h1;
            h2->prev = h1;
        }
        
        if (v1->hedge == nullptr) {
            v1->hedge = h1;
        }
        else {
            HalfEdgeEntry *e1 = v1_insertion_point,
                          *e2 = v1_insertion_point->prev;
                          
            e2->next = h1;
            h1->prev = e2;

            h2->next = e1;
            e1->prev = h2;
        }
        
        if (v2->hedge == nullptr) {
            v2->hedge = h2;
        }
        else {
            HalfEdgeEntry *e1 = v2_insertion_point,
                          *e2 = v2_insertion_point->prev;
            
            h1->next = e1;
            e1->prev = h1;

            e2->next = h2;
            h2->prev = e2;
        }
        
        return edge;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::FaceHandle
    HALF_EDGE_MESH_C::createFace( const CVertexHandle *vertexes, size_t count )
    {
        if (count >= MAX_VERTEXES_PER_FACE) {
            HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createFace - To many vertexes.\n");
        }
        
        CHalfEdgeHandle halfEdges[MAX_VERTEXES_PER_FACE] = {};
        
        for (size_t i=0; i < count; ++i) {
            if (vertexes[i] == nullptr) {
                HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createFace - Invalid vertexes.\n");
                return nullptr;
            }
        }
        for (size_t i=0; i < count; ++i) {
            CVertexHandle v1 = vertexes[i],
                          v2 = vertexes[(i+1) % count];
            
            EdgeEntry *edge = asEntry(createEdge(v1, v2));
            if (edge == nullptr) {
                HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createFace - Failed to create edge.\n");
                return nullptr;
            }
            
            HalfEdgeEntry *hedge = edge->hedge;
            if (hedge->vertex != v2) {
                hedge = hedge->pair;
                HALF_EDGE_MESH_ASSERT(hedge->vertex == v2);
            }
            halfEdges[i] = hedge;
        }
        
        return createFace(halfEdges, count);
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::FaceHandle
    HALF_EDGE_MESH_C::createFace( const CHalfEdgeHandle *halfEdges_, size_t count )
    {
        if (count >= MAX_VERTEXES_PER_FACE) {
            HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createFace - To many vertexes.\n");
            return nullptr;
        }
        
        HalfEdgeEntry *halfEdges[MAX_VERTEXES_PER_FACE];
        
        for (size_t i=0; i < count; ++i) {
            halfEdges[i] = asEntry(halfEdges_[i]);
            if (halfEdges[i] == nullptr) {
                HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createFace - Invalid half edges.\n");
                return nullptr; 
            }
        }
        
        // check that the half edges can form a valid face
        for (size_t i=0; i < count; ++i) {
            HalfEdgeEntry *e1 = halfEdges[i],
                          *e2 = halfEdges[(i+1) % count];
                
            // Do they form a loop?
            if (e1->vertex != e2->pair->vertex) {
                HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createFace - The half edges don't form a loop.\n");
                return nullptr;
            }
            
            // Do they have a existing face?
            if (e1->face != nullptr) {
                HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createFace - The half edges have a existing face.\n");
                return nullptr;
            }
        }
        
        // Make all edges adjacent to each other (that is next points to the corrent half edge)
        for (size_t i=0; i < count; ++i) {
            HalfEdgeEntry *e1 = halfEdges[i],
                          *e2 = halfEdges[(i+1) % count];
                          
            if (makeAdjacent(e1, e2) == false) {
                HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createFace - Failed to make half edges adjacent.\n");
                return nullptr;
            }
        }
        
        // Create the face
        FaceEntry *face = mFaces.allocate();
        if (!face) {
            HALF_EDGE_MESH_PRINTF ("HalfEdgeMesh::createFace - Failed to allocate face.\n");
            return nullptr;
        }
        
        for (size_t i=0; i < count; ++i) {
            halfEdges[i]->face = face;
        }
        face->hedge = halfEdges[0];
        
        return face;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::CEdgeHandle 
    HALF_EDGE_MESH_C::findEdge( CVertexHandle v1, CVertexHandle v2 ) const
    {
        const HalfEdgeEntry *hedge = asEntry(findHalfEdge(v1, v2));
        if (hedge) {
            return hedge->edge;
        }
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::CHalfEdgeHandle 
    HALF_EDGE_MESH_C::findHalfEdge( CVertexHandle v1_, CVertexHandle v2 ) const
    {
        const VertexEntry *v1 = asEntry(v1_);
        
        if (v1 == nullptr || v2 == nullptr) {
            return nullptr;
        }
        
        const HalfEdgeEntry *head = v1->hedge;
        if (head == nullptr) return nullptr;
        
        const HalfEdgeEntry *current = head;
        do {
            if (current->vertex == v2) {
                return current;
            }
            
            current = current->pair->next;
        } while (current != head);
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::HalfEdgeEntry* 
    HALF_EDGE_MESH_C::findFreeHEdge( HalfEdgeEntry *hedge )
    {
        if (hedge == nullptr) {
            return nullptr;
        }
        
        const HalfEdgeEntry *root = hedge;
        do {
            if (hedge->face == nullptr) {
                return hedge;
            }
            hedge = hedge->pair->next;
        } while (root != hedge);
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::HalfEdgeEntry* 
    HALF_EDGE_MESH_C::findFreeHEdge( HalfEdgeEntry *before, HalfEdgeEntry *after )
    {
        HALF_EDGE_MESH_ASSERT(after->vertex == before->vertex);
        
        while (true) {
            after = after->next->pair;
            if (after == before) {
                return nullptr;
            }
            if (after->face == nullptr) {
                return after;
            }
        }
    }
    
    HALF_EDGE_MESH_T
    bool 
    HALF_EDGE_MESH_C::makeAdjacent( HalfEdgeEntry *in, HalfEdgeEntry *out )
    {
        HALF_EDGE_MESH_ASSERT(in != nullptr && out != nullptr);
        if (in->next == out) {
            // they are already adjacent
            return true;
        }
        
        HalfEdgeEntry *a = in->next,
                      *b = out->prev;

        HalfEdgeEntry *c = findFreeHEdge(in, out->pair);
        if (c == nullptr) return false;

        HalfEdgeEntry *d = c->next;

        in->next = out;
        out->prev = in;

        c->next = a;
        a->prev = c;

        b->next = d;
        d->prev = b;

        return true;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::CHalfEdgeHandle 
    HALF_EDGE_MESH_C::getFaceHalfEdge( CFaceHandle face_ ) const
    {
        const FaceEntry *face = asEntry(face_);
        if (face) {
            return face->hedge;
        }
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::CHalfEdgeHandle 
    HALF_EDGE_MESH_C::getEdgeHalfEdge( CEdgeHandle edge_ ) const
    {
        const EdgeEntry *edge = asEntry(edge_);
        if (edge) {
            return edge->hedge;
        }
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::CHalfEdgeHandle 
    HALF_EDGE_MESH_C::getVertexHalfEdge( CVertexHandle vertex_ ) const
    {
        const VertexEntry *vertex = asEntry(vertex_);
        if (vertex) {
            return vertex->hedge;
        }
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::CHalfEdgeHandle 
    HALF_EDGE_MESH_C::getHalfEdgePair( CHalfEdgeHandle edge_ ) const
    {
        const HalfEdgeEntry *edge = asEntry(edge_);
        if (edge) {
            return edge->pair;
        }
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::CHalfEdgeHandle 
    HALF_EDGE_MESH_C::getHalfEdgeNext( CHalfEdgeHandle edge_ ) const
    {
        const HalfEdgeEntry *edge = asEntry(edge_);
        if (edge) {
            return edge->next;
        }
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::CHalfEdgeHandle 
    HALF_EDGE_MESH_C::getHalfEdgePrev( CHalfEdgeHandle edge_ ) const
    {
        const HalfEdgeEntry *edge = asEntry(edge_);
        if (edge) {
            return edge->prev;
        }
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::CVertexHandle 
    HALF_EDGE_MESH_C::getHalfEdgeVertex( CHalfEdgeHandle edge_ ) const
    {
        const HalfEdgeEntry *edge = asEntry(edge_);
        if (edge) {
            return edge->vertex;
        }
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::CVertexHandle
    HALF_EDGE_MESH_C::getHalfEdgeVertexOrigin( CHalfEdgeHandle edge_ ) const
    {
        const HalfEdgeEntry *edge = asEntry(edge_);
        if (edge) {
            return edge->pair->vertex;
        }
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::CEdgeHandle 
    HALF_EDGE_MESH_C::getHalfEdgeEdge( CHalfEdgeHandle edge_ ) const
    {
        const HalfEdgeEntry *edge = asEntry(edge_);
        if (edge) {
            return edge->edge;
        }
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    typename HALF_EDGE_MESH_C::CFaceHandle 
    HALF_EDGE_MESH_C::getHalfEdgeFace( CHalfEdgeHandle edge_ ) const
    {
        const HalfEdgeEntry *edge = asEntry(edge_);
        if (edge) {
            return edge->face;
        }
        return nullptr;
    }
    
    HALF_EDGE_MESH_T
    std::pair<typename HALF_EDGE_MESH_C::CVertexHandle, typename HALF_EDGE_MESH_C::CVertexHandle>
    HALF_EDGE_MESH_C::getEdgeVertexes( CEdgeHandle edge_ ) const
    {
        const EdgeEntry *edge = asEntry(edge_);
        if (edge) {
            const HalfEdgeEntry *e1 = edge->hedge;
            const HalfEdgeEntry *e2 = e1->pair;
            
            return std::make_pair(
                e1->vertex,
                e2->vertex
            );
        }
        return std::make_pair(nullptr, nullptr);
    }
    
    HALF_EDGE_MESH_T
    std::pair<typename HALF_EDGE_MESH_C::CFaceHandle, typename HALF_EDGE_MESH_C::CFaceHandle> 
    HALF_EDGE_MESH_C::getEdgeFaces( CEdgeHandle edge_ ) const
    {
        const EdgeEntry *edge = asEntry(edge_);
        if (edge) {
            const HalfEdgeEntry *e1 = edge->hedge;
            const HalfEdgeEntry *e2 = e1->pair;
            
            return std::make_pair(
                e1->face,
                e2->face
            );
        }
        return std::make_pair(nullptr, nullptr);
        
    }
    
    HALF_EDGE_MESH_T
    bool 
    HALF_EDGE_MESH_C::isBorderFace( CFaceHandle face_ ) const
    {
        const FaceEntry *face = asEntry(face_);
        
        const HalfEdgeEntry *edge = face->hedge;
        do {
            if (edge->pair->face == nullptr) {
                return true;
            }
            edge = edge->next;
        } while (edge != face->hedge);
        return false;
    }
    
#ifndef HALF_EDGE_MESH_NO_ITERATORS

#define HALF_EDGE_MESH_CREATE_CORE_ITER( Name, FuncName, Member )           \
    HALF_EDGE_MESH_T                                                        \
    typename HALF_EDGE_MESH_C::Name##Range                                  \
    HALF_EDGE_MESH_C::FuncName() {                                          \
        return Name##Range(                                                 \
            Name##Iterator(Member.begin()),                                 \
            Name##Iterator(Member.end())                                    \
        );                                                                  \
    }                                                                       \
    HALF_EDGE_MESH_T                                                        \
    typename HALF_EDGE_MESH_C::C##Name##Range                               \
    HALF_EDGE_MESH_C::FuncName() const {                                    \
        return C##Name##Range(                                              \
            C##Name##Iterator(Member.begin()),                              \
            C##Name##Iterator(Member.end())                                 \
        );                                                                  \
    }
    
    HALF_EDGE_MESH_CREATE_CORE_ITER(Vertex, vertexes, mVertexes)
    HALF_EDGE_MESH_CREATE_CORE_ITER(Face, faces, mFaces)
    HALF_EDGE_MESH_CREATE_CORE_ITER(Edge, edges, mEdges)
    HALF_EDGE_MESH_CREATE_CORE_ITER(HalfEdge, halfEdges, mHalfEdges)
    
#undef HALF_EDGE_MESH_CREATE_CORE_ITER
#define HALF_EDGE_MESH_CREATE_ITER( Name, FuncName, Handle )                \
    HALF_EDGE_MESH_T                                                        \
    typename HALF_EDGE_MESH_C::Name##Range                                  \
    HALF_EDGE_MESH_C::FuncName( Handle handle ) {                           \
        return Name##Range(                                                 \
            Name##Iterator(asEntry(handle),0),                              \
            Name##Iterator(asEntry(handle),1)                               \
        );                                                                  \
    }                                                                       \
    HALF_EDGE_MESH_T                                                        \
    typename HALF_EDGE_MESH_C::C##Name##Range                               \
    HALF_EDGE_MESH_C::FuncName( Handle handle ) const {                     \
        return C##Name##Range(                                              \
            C##Name##Iterator(asEntry(handle),0),                           \
            C##Name##Iterator(asEntry(handle),1)                            \
        );                                                                  \
    }
    
    
    HALF_EDGE_MESH_CREATE_ITER(VertexVertex, vertexVertexes, CVertexHandle)
    HALF_EDGE_MESH_CREATE_ITER(VertexOutHalfEdge, vertexOutHalfEdges, CVertexHandle)
    HALF_EDGE_MESH_CREATE_ITER(VertexInHalfEdge, vertexInHalfEdges, CVertexHandle)
    HALF_EDGE_MESH_CREATE_ITER(VertexEdge, vertexEdges, CVertexHandle)
    HALF_EDGE_MESH_CREATE_ITER(VertexFace, vertexFaces, CVertexHandle)
    
    HALF_EDGE_MESH_CREATE_ITER(FaceVertex, faceVertexes, CFaceHandle)
    HALF_EDGE_MESH_CREATE_ITER(FaceHalfEdge, faceHalfEdges, CFaceHandle)
    HALF_EDGE_MESH_CREATE_ITER(FaceEdge, faceEdges, CFaceHandle)
    HALF_EDGE_MESH_CREATE_ITER(FaceFace, faceFaces, CFaceHandle)
    
#undef HALF_EDGE_MESH_CREATE_ITER

#endif


    
#undef HALF_EDGE_MESH_T
#undef HALF_EDGE_MESH_C
}
