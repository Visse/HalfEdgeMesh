#pragma once
#define HALF_EDGE_MESH_H

#include <type_traits>
#include <bitset>
#include <new>




/*
 *
 */


// is we compiling with debugging and haven't defined HALF_EDGE_MESH_DEBUG
// enable HALF_EDGE_MESH_DEBUG
#if !(defined(NDEBUG) || defined(HALF_EDGE_MESH_DEBUG))
#   define HALF_EDGE_MESH_DEBUG 1
#endif

#if HALF_EDGE_MESH_DEBUG
#   ifndef HALF_EDGE_MESH_ASSERT
#       include <cassert>
#       define HALF_EDGE_MESH_ASSERT(X) assert(X)
#   endif
#   ifndef HALF_EDGE_MESH_PRINTF
#       define HALF_EDGE_MESH_PRINTF(  ... ) printf( __VA_ARGS__ )
#   endif
#else
#   ifndef HALF_EDGE_MESH_ASSERT
#       define HALF_EDGE_MESH_ASSERT(X)
#   endif
#   ifndef HALF_EDGE_MESH_PRINTF
#       define HALF_EDGE_MESH_PRINTF( ... )
#   endif
#endif

#include <cassert>

namespace HalfEdgeMesh
{
    namespace internal
    {
        template< typename T, size_t BlockSize >
        class SlabAllocator;
    }
    
    template< typename FaceData, typename EdgeData, typename HalfEdgeData, typename VertexData >
    class HalfEdgeMesh {
    private:
        // Private handle types to offer better type protection
        // We also derive from the data type, instead of having it as a member
        // to be able to collapse them if they are empty (members must have a unique address)
        // 
        struct FaceHandle_t : 
            public FaceData 
        {};
        struct VertexHandle_t : 
            public VertexData 
        {};
        struct EdgeHandle_t :
            public EdgeData
        {};
        struct HalfEdgeHandle_t :
            public HalfEdgeData
        {};
        
        static const int MAX_VERTEXES_PER_FACE = 16;
        
        
    public:
        // public typedef's for the handle types
        using FaceHandle     = FaceHandle_t*;
        using VertexHandle   = VertexHandle_t*;
        using EdgeHandle     = EdgeHandle_t*;
        using HalfEdgeHandle = HalfEdgeHandle_t*;
        
        
        using CFaceHandle     = const FaceHandle_t*;
        using CVertexHandle   = const VertexHandle_t*;
        using CEdgeHandle     = const EdgeHandle_t*;
        using CHalfEdgeHandle = const HalfEdgeHandle_t*;

    public:
        HalfEdgeMesh() = default;
        ~HalfEdgeMesh() = default;
        
        HalfEdgeMesh( HalfEdgeMesh&& move ) = default;
        
        // for now no copy
        HalfEdgeMesh( const HalfEdgeMesh& ) = delete;
        HalfEdgeMesh& operator = ( const HalfEdgeMesh& ) = delete;
        
    public:
        VertexHandle createVertex();
        EdgeHandle createEdge( CVertexHandle v1, CVertexHandle v2 );
        
        FaceHandle createFace( const CVertexHandle *vertexes, size_t count );
       // FaceHandle createFace( EdgeHandle *edges, size_t count );
        FaceHandle createFace( const CHalfEdgeHandle *halfEdges, size_t count );
        
        EdgeHandle findEdge( CVertexHandle v1, CVertexHandle v2 ) {return const_cast<EdgeHandle>(asConst()->findEdge(v1,v2));}
        CEdgeHandle findEdge( CVertexHandle v1, CVertexHandle v2 ) const;
        
        HalfEdgeHandle findHalfEdge( CVertexHandle v1, CVertexHandle v2 ) {return const_cast<HalfEdgeHandle>(asConst()->findHalfEdge(v1,v2));}
        CHalfEdgeHandle findHalfEdge( CVertexHandle v1, CVertexHandle v2 ) const;
        
        
        HalfEdgeHandle  getFaceHalfEdge( CFaceHandle face ) {return const_cast<HalfEdgeHandle>(asConst()->getFaceHalfEdge(face));}
        CHalfEdgeHandle getFaceHalfEdge( CFaceHandle face ) const;
        
        HalfEdgeHandle  getEdgeHalfEdge( CEdgeHandle edge ) {return const_cast<HalfEdgeHandle>(asConst()->getEdgeHalfEdge(edge));}
        CHalfEdgeHandle getEdgeHalfEdge( CEdgeHandle edge ) const;
        
        HalfEdgeHandle  getVertexHalfEdge( CVertexHandle vertex ) {return const_cast<HalfEdgeHandle>(asConst()->getVertexHalfEdge(vertex));}
        CHalfEdgeHandle getVertexHalfEdge( CVertexHandle vertex ) const;
        
        HalfEdgeHandle  getHalfEdgePair( CHalfEdgeHandle edge ) {return const_cast<HalfEdgeHandle>(asConst()->getHalfEdgePair(edge));}
        CHalfEdgeHandle getHalfEdgePair( CHalfEdgeHandle edge ) const;
        
        HalfEdgeHandle  getHalfEdgeNext( CHalfEdgeHandle edge ) {return const_cast<HalfEdgeHandle>(asConst()->getHalfEdgeNext(edge));}
        CHalfEdgeHandle getHalfEdgeNext( CHalfEdgeHandle edge ) const;
        
        HalfEdgeHandle  getHalfEdgePrev( CHalfEdgeHandle edge ) {return const_cast<HalfEdgeHandle>(asConst()->getHalfEdgePrev(edge));}
        CHalfEdgeHandle getHalfEdgePrev( CHalfEdgeHandle edge ) const;
        
        VertexHandle  getHalfEdgeVertex( CHalfEdgeHandle edge ) {return const_cast<VertexHandle>(asConst()->getHalfEdgeVertex(edge));}
        CVertexHandle getHalfEdgeVertex( CHalfEdgeHandle edge ) const;
        
        VertexHandle  getHalfEdgeVertexOrigin( CHalfEdgeHandle edge ) {return const_cast<VertexHandle>(asConst()->getHalfEdgeVertexOrigin(edge));}
        CVertexHandle getHalfEdgeVertexOrigin( CHalfEdgeHandle edge ) const;
        
        EdgeHandle  getHalfEdgeEdge( CHalfEdgeHandle edge ) {return const_cast<EdgeHandle>(asConst()->getHalfEdgeEdge(edge));}
        CEdgeHandle getHalfEdgeEdge( CHalfEdgeHandle edge ) const;
        
        FaceHandle  getHalfEdgeFace( CHalfEdgeHandle edge ) {return const_cast<FaceHandle>(asConst()->getHalfEdgeFace(edge));}
        CFaceHandle getHalfEdgeFace( CHalfEdgeHandle edge ) const;
        
        std::pair<VertexHandle, VertexHandle> getEdgeVertexes( CEdgeHandle edge ) {
            auto pair = asConst()->getEdgeVertexes(edge);
            return std::make_pair(
                const_cast<VertexHandle>(pair.first),
                const_cast<VertexHandle>(pair.second)
            );
        }
        std::pair<CVertexHandle, CVertexHandle> getEdgeVertexes( CEdgeHandle edge ) const;
        
        // WARNING: the faces may be null
        std::pair<FaceHandle, FaceHandle> getEdgeFaces( CEdgeHandle edge ) {
            auto pair = asConst()->getEdgeFaces(edge);
            return std::make_pair(
                const_cast<FaceHandle>(pair.first),
                const_cast<FaceHandle>(pair.second)
            );
        }
        std::pair<CFaceHandle, CFaceHandle> getEdgeFaces( CEdgeHandle edge ) const;
        
        // A face is on the border if it has an edge without a neighbor
        bool isBorderFace( CFaceHandle face ) const;
        
#ifndef HALF_EDGE_MESH_NO_ITERATORS
        class VertexIterator;
        class EdgeIterator;
        class HalfEdgeIterator;
        class FaceIterator;
        
        class VertexVertexIterator;
        class VertexOutHalfEdgeIterator;
        class VertexInHalfEdgeIterator;
        class VertexEdgeIterator;
        class VertexFaceIterator;
        
        class FaceVertexIterator;
        class FaceHalfEdgeIterator;
        class FaceEdgeIterator;
        class FaceFaceIterator;
        
        class VertexRange;
        class HalfEdgeRange;
        class EdgeRange;
        class FaceRange;
        
        class VertexVertexRange;
        class VertexOutHalfEdgeRange;
        class VertexInHalfEdgeRange;
        class VertexEdgeRange;
        class VertexFaceRange;
        
        class FaceVertexRange;
        class FaceHalfEdgeRange;
        class FaceEdgeRange;
        class FaceFaceRange;
    
        class CVertexIterator;
        class CEdgeIterator;
        class CHalfEdgeIterator;
        class CFaceIterator;
        
        class CVertexVertexIterator;
        class CVertexOutHalfEdgeIterator;
        class CVertexInHalfEdgeIterator;
        class CVertexEdgeIterator;
        class CVertexFaceIterator;
        
        class CFaceVertexIterator;
        class CFaceHalfEdgeIterator;
        class CFaceEdgeIterator;
        class CFaceFaceIterator;
        
        class CVertexRange;
        class CHalfEdgeRange;
        class CEdgeRange;
        class CFaceRange;
        
        class CVertexVertexRange;
        class CVertexOutHalfEdgeRange;
        class CVertexInHalfEdgeRange;
        class CVertexEdgeRange;
        class CVertexFaceRange;
        
        class CFaceVertexRange;
        class CFaceHalfEdgeRange;
        class CFaceEdgeRange;
        class CFaceFaceRange;
        
        VertexRange vertexes();
        CVertexRange vertexes() const;
        
        HalfEdgeRange halfEdges();
        CHalfEdgeRange halfEdges() const;
        
        EdgeRange edges();
        CEdgeRange edges() const;
        
        FaceRange faces();
        CFaceRange faces() const;
        
        VertexVertexRange vertexVertexes( CVertexHandle vertex );
        CVertexVertexRange vertexVertexes( CVertexHandle vertex ) const;
        
        VertexOutHalfEdgeRange vertexOutHalfEdges( CVertexHandle vertex );
        CVertexOutHalfEdgeRange vertexOutHalfEdges( CVertexHandle vertex ) const;
        
        VertexInHalfEdgeRange vertexInHalfEdges( CVertexHandle vertex );
        CVertexInHalfEdgeRange vertexInHalfEdges( CVertexHandle vertex ) const;
        
        VertexEdgeRange vertexEdges( CVertexHandle vertex );
        CVertexEdgeRange vertexEdges( CVertexHandle vertex ) const;
        
        VertexFaceRange vertexFaces( CVertexHandle vertex );
        CVertexFaceRange vertexFaces( CVertexHandle vertex ) const;
        
        FaceVertexRange faceVertexes( CFaceHandle face );
        CFaceVertexRange faceVertexes( CFaceHandle face ) const;
        
        FaceHalfEdgeRange faceHalfEdges( CFaceHandle face );
        CFaceHalfEdgeRange faceHalfEdges( CFaceHandle face ) const;
        
        FaceEdgeRange faceEdges( CFaceHandle face );
        CFaceEdgeRange faceEdges( CFaceHandle face ) const;
        
        FaceFaceRange faceFaces( CFaceHandle face );
        CFaceFaceRange faceFaces( CFaceHandle face ) const;
        
    private:
        template< typename Iterator, typename Entry, typename Handle, bool const_ >
        class BaseHalfEdgeIterator;
        
#endif // !HALF_EDGE_MESH_NO_ITERATORS
        
        
    private:
        struct HalfEdgeEntry;
        struct FaceEntry :
            public FaceHandle_t
        {
            HalfEdgeEntry *hedge = nullptr;
        };
        
        struct EdgeEntry :
            public EdgeHandle_t
        {
            HalfEdgeEntry *hedge = nullptr;
        };
        
        struct VertexEntry :
            public VertexHandle_t
        {
            HalfEdgeEntry *hedge = nullptr;
        };
        
        struct HalfEdgeEntry :
            public HalfEdgeHandle_t
        {
            HalfEdgeEntry *pair = nullptr,
                          *next = nullptr,
                          *prev = nullptr;
            // The vertex this halfedge points to
            //    (use pair->vertex) for the vertex that this edge originates from
            VertexEntry *vertex = nullptr;
            EdgeEntry *edge = nullptr;
            FaceEntry *face = nullptr;
            
            // Used in non-constant iterators when they need to return a this half edge (taken as a const pointer)
            HalfEdgeEntry* asNonConst() const {return const_cast<HalfEdgeEntry*>(this);}
        };
        
    private:
        // Non const to prevent const function calling it
        const HalfEdgeMesh* asConst() /*const*/ {
            return this;
        }
        VertexEntry* asEntry( CVertexHandle handle ) const {
            return static_cast<VertexEntry*>(const_cast<VertexHandle>(handle));
        }
        HalfEdgeEntry* asEntry( CHalfEdgeHandle handle ) const {
            return static_cast<HalfEdgeEntry*>(const_cast<HalfEdgeHandle>(handle));
        }
        EdgeEntry* asEntry( CEdgeHandle handle ) const {
            return static_cast<EdgeEntry*>(const_cast<EdgeHandle>(handle));
        }
        FaceEntry* asEntry( CFaceHandle handle ) const {
            return static_cast<FaceEntry*>(const_cast<FaceHandle>(handle));
        }
        
        
        HalfEdgeEntry* findFreeHEdge( HalfEdgeEntry *hedge );
        HalfEdgeEntry* findFreeHEdge( HalfEdgeEntry *before, HalfEdgeEntry *after );
        bool makeAdjacent( HalfEdgeEntry *in, HalfEdgeEntry *out );
        
    private:
        static const int FACE_BLOCK_SIZE = 64,
                         EDGE_BLOCK_SIZE = 64,
                         VERTEX_BLOCK_SIZE = 64,
                         HALF_EDGE_BLOCK_SIZE = 128;
        internal::SlabAllocator<FaceEntry, FACE_BLOCK_SIZE> mFaces;
        internal::SlabAllocator<EdgeEntry, EDGE_BLOCK_SIZE> mEdges;
        internal::SlabAllocator<VertexEntry, VERTEX_BLOCK_SIZE> mVertexes;
        internal::SlabAllocator<HalfEdgeEntry, HALF_EDGE_BLOCK_SIZE> mHalfEdges;
    };
    
    namespace internal
    {
        template< typename T, size_t BlockSize >
        class SlabAllocator {
            union Slab {
                T data;
            };
            struct SlabBlock {
                typedef typename std::aligned_storage<sizeof(T), alignof(T)>::type StorageType;
                
                SlabBlock *next = nullptr;
                // Needed by bidirectional iterators
                SlabBlock *prev = nullptr;
                std::bitset<BlockSize> freeSlabs;
                
                StorageType data[BlockSize];
                
                SlabBlock() {
                    // Mark all slabs as free
                    freeSlabs.set();
                }
                ~SlabBlock() {
                    // No need to do anything on trivially destructible types
                    if (std::is_trivially_destructible<T>::value) return;
                    // No need to do anything if no allocations was made
                    if (freeSlabs.all()) return;
                    
                    for (size_t i=0; i < BlockSize; ++i) {
                        if (freeSlabs.test(i)) continue;
                        
                        // destruct object
                        get(i)->~T();
                    }
                }
                
                T* get( size_t i ) {
                    HALF_EDGE_MESH_ASSERT(i < BlockSize);
                    HALF_EDGE_MESH_ASSERT(isFree(i) == false);
                    return (T*) &data[i];
                }
                T* allocate() {
                    size_t idx = BlockSize;
                #if __GNUC__ && !__clang__
                    idx = freeSlabs._Find_first();
                #else  /// @todo figure out if visual studio provides any better way to find the first bit set
                    if (freeSlabs.none()) return nullptr;
                    
                    for (size_t i=0; i < BlockSize; ++i) {
                        if (freeSlabs.test(i)) {
                            idx = i;
                            break;
                        }
                    }
                #endif
                    
                    if (idx == BlockSize) return nullptr;
                    
                    freeSlabs.set(idx, false);
                    
                    T *ptr = get(idx);
                    new (ptr) T;
                    
                    return ptr;
                }
            
                bool contains( T *ptr ) {
                    return (ptr >= (const T*)std::begin(data) && ptr < (const T*)std::end(data));
                }
                
                void free( T *ptr ) {
                    HALF_EDGE_MESH_ASSERT (contains(ptr));
                    
                    size_t idx = ptr - (const T*)std::begin(data);
                    HALF_EDGE_MESH_ASSERT (isFree(idx) == false);
                    
                    ptr->~T();
                    freeSlabs.set(idx, true);
                }
            
                bool isFree( size_t i ) {
                    HALF_EDGE_MESH_ASSERT(i < BlockSize);
                    return freeSlabs.test(i);
                }
            };
            
            SlabBlock *mFirstBlock = nullptr,
                      *mLastBlock = nullptr;
        public:
            SlabAllocator() 
            {}
            ~SlabAllocator() {
                SlabBlock *block = mFirstBlock;
                while (block) {
                    SlabBlock *next = block->next;
                    delete block;
                    block = next;
                }
            }
            
            SlabAllocator( SlabAllocator &&move ) {
                std::swap(mFirstBlock, move.mFirstBlock);
                std::swap(mLastBlock, move.mLastBlock);
            }
            
            SlabAllocator( const SlabAllocator& ) = delete;
            SlabAllocator& operator = ( const SlabAllocator& ) = delete;
            
            
            T* allocate() {
                SlabBlock *block = mFirstBlock;
                
                while (block) {
                    if (T *ptr = block->allocate()) {
                        return ptr;
                    }
                    
                    block = block->next;
                }
                
                block = new (std::nothrow) SlabBlock;
                if (!block) return nullptr;
                
                if (mFirstBlock == nullptr) {
                    HALF_EDGE_MESH_ASSERT(mLastBlock == nullptr);
                    mFirstBlock = block;
                    mLastBlock = block;
                }
                else {
                    block->next = mFirstBlock;
                    mFirstBlock->prev = block;
                    
                    mFirstBlock = block;
                }
                
                return mFirstBlock->allocate();
            }
            
            void free( T *ptr ) {
                if (ptr == nullptr) return;
                
                SlabBlock *block = mFirstBlock;
                
                while (block) {
                    if (block->contains(ptr)) {
                        block->free(ptr);
                        return;
                    }
                    block = block->next;
                }
                HALF_EDGE_MESH_ASSERT (false && "Failed to free pointer!");
            }
            
#ifndef HALF_EDGE_MESH_NO_ITERATORS

            template<bool const_>
            class base_iterator;
            
            using iterator = base_iterator<false>;
            using const_iterator = base_iterator<true>;
            
            iterator begin() {
                if (mFirstBlock) {
                    iterator iter = iterator(mFirstBlock, 0);
                    if (mFirstBlock->isFree(0)) {
                        ++iter;
                    }
                    return iter;
                }
                else {
                    return end();
                }
            }
            iterator end() {
                return iterator(mLastBlock, BlockSize);
            }
            
            const_iterator begin() const {
                if (mFirstBlock) {
                    const_iterator iter = const_iterator(mFirstBlock, 0);
                    if (mFirstBlock->isFree(0)) {
                        ++iter;
                    }
                    return iter;
                }
                else {
                    return end();
                }
            }
            const_iterator end() const {
                return const_iterator(mLastBlock, BlockSize);
            }
#endif
        };
    }
}


#ifndef HALF_EDGE_MESH_NO_ITERATORS

#include <iterator>

namespace HalfEdgeMesh
{
    namespace internal
    {
        template< typename T, bool const_ >
        struct ConstType {
            using Type = T;
        };
        template<typename T>
        struct ConstType<T,true> {
            using Type = const T;
        };
            
        template< typename IteratorType, typename ValueType, typename PointerType=ValueType*, typename ReferenceType = ValueType& >
        struct IteratorAdopter
        {
            using value_type = ValueType;
            using pointer = PointerType;
            using reference = ValueType&;
            using difference_type = ptrdiff_t;
            using iterator_category = std::bidirectional_iterator_tag;
            
            IteratorType* iterator() {
                return static_cast<IteratorType*>(this);
            }
            const IteratorType* iterator() const {
                return static_cast<const IteratorType*>(this);
            }
            
            friend bool operator == ( const IteratorAdopter &lhs, const IteratorAdopter &rhs ) {
                return lhs.iterator()->equal(*rhs.iterator());
            }
            friend bool operator != ( const IteratorAdopter &lhs, const IteratorAdopter &rhs ) {
                return !lhs.iterator()->equal(*rhs.iterator());
            }
            
            IteratorType& operator ++ () {
                iterator()->increment();
                return *iterator();
            }
            
            IteratorType operator ++ (int) {
                IteratorType copy(*iterator());
                iterator()->increment();
                return *copy;
            }
            
            IteratorType& operator -- () {
                iterator()->decrement();
                return *iterator();
            }
            IteratorType operator -- (int) {
                IteratorType copy(*iterator());
                iterator()->decrement();
                return copy;
            }
        };
            
        template< typename T, size_t BlockSize >
        template< bool const_>
        class SlabAllocator<T,BlockSize>::base_iterator :
            public IteratorAdopter<base_iterator<const_>, typename ConstType<T,const_>::Type>
        {
            SlabBlock *mBlock = nullptr;
            size_t mSlab = 0;
            
            using Type = typename ConstType<T,const_>::Type;
            
        public:
            base_iterator() = default;
            base_iterator( SlabBlock *block, size_t slab ) :
                mBlock(block),
                mSlab(slab)
            {}
            base_iterator( const base_iterator& ) = default;
            base_iterator& operator = ( const base_iterator& ) = default;
            
            template< bool C, typename=typename std::enable_if<const_ || !C>::type>
            base_iterator( const base_iterator<C> &copy ) :
                mBlock(copy.mBlock),
                mSlab(copy.mSlab)
            {}
            
            template< bool C, typename=typename std::enable_if<const_ || !C>::type>
            base_iterator& operator = ( const base_iterator<C> &copy ) {
                mBlock = copy.mBlock;
                mSlab = mSlab;
            }
            
            bool equal( const base_iterator &other ) const {
                return mBlock == other.mBlock && 
                       mSlab  == other.mSlab;
            }
            
            Type& operator * () const {
                HALF_EDGE_MESH_ASSERT(mBlock && mSlab < BlockSize);
                return *mBlock->get(mSlab);
            }
            
            Type* operator -> () const {
                HALF_EDGE_MESH_ASSERT(mBlock && mSlab < BlockSize);
                return mBlock->get(mSlab);
            }
            
            void increment() {
                HALF_EDGE_MESH_ASSERT(mBlock);
                HALF_EDGE_MESH_ASSERT(mSlab < BlockSize);
                do {
                    mSlab++;
                    if (mBlock->next && mSlab == BlockSize) {
                        mSlab = 0;
                        mBlock = mBlock->next;
                    }
                    
                } while (mSlab != BlockSize && mBlock->isFree(mSlab));
                
            }
            
            void decrement() {
                HALF_EDGE_MESH_ASSERT(mBlock);
                
                do {
                    if (mSlab == 0 && mBlock->prev) {
                        mSlab = BlockSize;
                        mBlock = mBlock->prev;
                    }
                    if (mSlab > 0) mSlab--;
                    else break;
                } while (mBlock->isFree(mSlab));
            }
            
        };
    }
        
    template< typename FaceData, typename EdgeData, typename HalfEdgeData, typename VertexData >
    template< typename IteratorType, typename Entry, typename Handle, bool const_ >
    class HalfEdgeMesh<FaceData, EdgeData, HalfEdgeData, VertexData>::BaseHalfEdgeIterator :
        public internal::IteratorAdopter<IteratorType, typename internal::ConstType<Handle,const_>::Type*, typename internal::ConstType<Handle,const_>::Type*>
    {
        const HalfEdgeEntry *mRoot = nullptr,
                            *mHead = nullptr;
        int mLap = 0;
    protected:   
        using HandleType = typename internal::ConstType<Handle,const_>::Type*;
        using EntryType = Entry*;
    public:
        BaseHalfEdgeIterator() = default;
        BaseHalfEdgeIterator( const BaseHalfEdgeIterator& ) = default;
        BaseHalfEdgeIterator& operator = ( const BaseHalfEdgeIterator& ) = default;
        
        template< bool C, typename=typename std::enable_if<const_ || !C>::type>
        BaseHalfEdgeIterator( const BaseHalfEdgeIterator<IteratorType,Entry,Handle,C> &copy ) :
            mRoot(copy.mRoot),
            mHead(copy.mHead),
            mLap(copy.mLap)
        {}
        
        template< bool C, typename=typename std::enable_if<const_ || !C>::type>
        BaseHalfEdgeIterator& operator = ( const BaseHalfEdgeIterator<IteratorType,Entry,Handle,C> &copy )
        {
            mRoot = copy.mRoot;
            mHead = copy.mHead;
            mLap = copy.mLap;
        }
        
        BaseHalfEdgeIterator( EntryType root, int lap=0 ) :
            mRoot(root->hedge),
            mHead(root->hedge),
            mLap(lap)
        {
            if (mHead == nullptr) {
                mLap = 1;
                return;
            }
            while (mLap == 0 && IteratorType::get(mHead) == nullptr) {
                mHead = IteratorType::next(mHead);
                if (mHead == mRoot) {
                    mLap = 1;
                }
            }
        }
        
        bool equal( const BaseHalfEdgeIterator &other ) const {
            if (mLap != other.mLap) return false;
            if (mHead && other.mHead) {
                return mHead == other.mHead;
            }
            return true;
        }
        
        HandleType operator * () const {
            HALF_EDGE_MESH_ASSERT(mHead);
            return IteratorType::get(mHead);
        }
        
        HandleType operator -> () const {
            HALF_EDGE_MESH_ASSERT(mHead);
            return IteratorType::get(mHead);
        }
        
        void increment() {
            HALF_EDGE_MESH_ASSERT(mHead);
            do {
                mHead = IteratorType::next(mHead);
                if (mHead == mRoot) {
                    mLap++;
                    return;
                }
            } while (IteratorType::get(mHead) == nullptr);
        }
        void decrement() {
            HALF_EDGE_MESH_ASSERT(mHead);
            do {
                if (mHead == mRoot) {
                    mLap--;
                    // For supporting decrement end() iterator
                    mHead = IteratorType::prev(mHead);
                    return;
                }
                mHead = IteratorType::prev(mHead);
            } while  (IteratorType::get(mHead) == nullptr);
        }
    };
        
#define HALF_EDGE_MESH_CORE_ITER(Name, Handle, Member, iterator)                                          \
    template< typename FaceData, typename EdgeData, typename HalfEdgeData, typename VertexData >          \
    class HalfEdgeMesh<FaceData, EdgeData, HalfEdgeData, VertexData>::Name :                              \
        public internal::IteratorAdopter<Name, Handle>                                                    \
    {                                                                                                     \
        using Iterator = typename decltype(Member)::iterator;                                             \
        Iterator mIter;                                                                                   \
    public:                                                                                               \
        Name() = default;                                                                                 \
        Name( Iterator iter ) :                                                                           \
            mIter(iter)                                                                                   \
        {}                                                                                                \
        Name( const Name &copy ) :                                                                        \
            mIter(copy.mIter)                                                                             \
        {}                                                                                                \
        bool equal( const Name &other ) const {                                                           \
            return mIter == other.mIter;                                                                  \
        }                                                                                                 \
        void increment() {                                                                                \
            ++mIter;                                                                                      \
        }                                                                                                 \
        void decrement() {                                                                                \
            --mIter;                                                                                      \
        }                                                                                                 \
        Handle operator * () const {                                                                      \
            return &*mIter;                                                                               \
        }                                                                                                 \
        Handle operator -> () const {                                                                     \
            return &*mIter;                                                                               \
        }                                                                                                 \
    };
        
    HALF_EDGE_MESH_CORE_ITER(VertexIterator, VertexHandle, mVertexes, iterator);
    HALF_EDGE_MESH_CORE_ITER(EdgeIterator, EdgeHandle, mEdges, iterator);
    HALF_EDGE_MESH_CORE_ITER(HalfEdgeIterator, HalfEdgeHandle, mHalfEdges, iterator);
    HALF_EDGE_MESH_CORE_ITER(FaceIterator, FaceHandle, mFaces, iterator);
    
    HALF_EDGE_MESH_CORE_ITER(CVertexIterator, CVertexHandle, mVertexes, const_iterator);
    HALF_EDGE_MESH_CORE_ITER(CEdgeIterator, CEdgeHandle, mEdges, const_iterator);
    HALF_EDGE_MESH_CORE_ITER(CHalfEdgeIterator, CHalfEdgeHandle, mHalfEdges, const_iterator);
    HALF_EDGE_MESH_CORE_ITER(CFaceIterator, CFaceHandle, mFaces, const_iterator);
    
#undef HALF_EDGE_MESH_CORE_ITER

#define HALF_EDGE_MESH_HEDGE_ITER(Name, Type, Handle, const_, inc, dec, attr)                                       \
    template< typename FaceData, typename EdgeData, typename HalfEdgeData, typename VertexData >                    \
    class HalfEdgeMesh<FaceData, EdgeData, HalfEdgeData, VertexData>::Name :                                        \
        public BaseHalfEdgeIterator<Name, Type, Handle, const_>                                                     \
    {                                                                                                               \
        using HandleType = typename BaseHalfEdgeIterator<Name, Type, Handle, const_>::HandleType;                   \
    public:                                                                                                         \
        using BaseHalfEdgeIterator<Name, Type, Handle, const_>::BaseHalfEdgeIterator;                               \
                                                                                                                    \
        static HandleType get( const HalfEdgeEntry *head ) {                                                        \
            return head attr;                                                                                       \
        }                                                                                                           \
        static const HalfEdgeEntry* next( const HalfEdgeEntry *head ) {                                             \
            return head inc;                                                                                        \
        }                                                                                                           \
        static const HalfEdgeEntry* prev( const HalfEdgeEntry *head ) {                                             \
            return head dec;                                                                                        \
        }                                                                                                           \
    };
    
    
    HALF_EDGE_MESH_HEDGE_ITER(VertexVertexIterator, VertexEntry, VertexHandle_t, false, ->pair->next, ->prev->pair, ->vertex);
    HALF_EDGE_MESH_HEDGE_ITER(VertexOutHalfEdgeIterator, VertexEntry, HalfEdgeHandle_t, false, ->pair->next, ->prev->pair, ->asNonConst() );
    HALF_EDGE_MESH_HEDGE_ITER(VertexInHalfEdgeIterator, VertexEntry, HalfEdgeHandle_t, false, ->pair->next, ->prev->pair, ->pair );
    HALF_EDGE_MESH_HEDGE_ITER(VertexEdgeIterator, VertexEntry, EdgeHandle_t, false, ->pair->next, ->prev->pair, ->edge );
    HALF_EDGE_MESH_HEDGE_ITER(VertexFaceIterator, VertexEntry, FaceHandle_t, false, ->pair->next, ->prev->pair, ->face );
    
    HALF_EDGE_MESH_HEDGE_ITER(FaceVertexIterator, FaceEntry, VertexHandle_t, false, ->next, ->prev, ->vertex );
    HALF_EDGE_MESH_HEDGE_ITER(FaceHalfEdgeIterator, FaceEntry, HalfEdgeHandle_t, false, ->next, ->prev, ->asNonConst() );
    HALF_EDGE_MESH_HEDGE_ITER(FaceEdgeIterator, FaceEntry, EdgeHandle_t, false, ->next, ->prev, ->edge );
    HALF_EDGE_MESH_HEDGE_ITER(FaceFaceIterator, FaceEntry, FaceHandle_t, false, ->next, ->prev, ->pair->face );

    
    HALF_EDGE_MESH_HEDGE_ITER(CVertexVertexIterator, VertexEntry, VertexHandle_t, true, ->pair->next, ->prev->pair, ->vertex);
    HALF_EDGE_MESH_HEDGE_ITER(CVertexOutHalfEdgeIterator, VertexEntry, HalfEdgeHandle_t, true, ->pair->next, ->prev->pair, );
    HALF_EDGE_MESH_HEDGE_ITER(CVertexInHalfEdgeIterator, VertexEntry, HalfEdgeHandle_t, true, ->pair->next, ->prev->pair, ->pair );
    HALF_EDGE_MESH_HEDGE_ITER(CVertexEdgeIterator, VertexEntry, EdgeHandle_t, true, ->pair->next, ->prev->pair, ->edge );
    HALF_EDGE_MESH_HEDGE_ITER(CVertexFaceIterator, VertexEntry, FaceHandle_t, true, ->pair->next, ->prev->pair, ->face );
    
    HALF_EDGE_MESH_HEDGE_ITER(CFaceVertexIterator, FaceEntry, VertexHandle_t, true, ->next, ->prev, ->vertex );
    HALF_EDGE_MESH_HEDGE_ITER(CFaceHalfEdgeIterator, FaceEntry, HalfEdgeHandle_t, true, ->next, ->prev, );
    HALF_EDGE_MESH_HEDGE_ITER(CFaceEdgeIterator, FaceEntry, EdgeHandle_t, true, ->next, ->prev, ->edge );
    HALF_EDGE_MESH_HEDGE_ITER(CFaceFaceIterator, FaceEntry, FaceHandle_t, true, ->next, ->prev, ->pair->face );
 
#undef HALF_EDGE_MESH_HEDGE_ITER
    
#define HALF_EDGE_MESH_RANGE(Name)                                                                        \
    template< typename FaceData, typename EdgeData, typename HalfEdgeData, typename VertexData >          \
    class HalfEdgeMesh<FaceData, EdgeData, HalfEdgeData, VertexData>::Name##Range                         \
    {                                                                                                     \
        Name##Iterator mBegin, mEnd;                                                                      \
    public:                                                                                               \
        Name##Range() = default;                                                                          \
        Name##Range( const Name##Range& ) = default;                                                      \
        Name##Range& operator = ( const Name##Range& ) = default;                                         \
        Name##Range( Name##Iterator begin, Name##Iterator end ) :                                         \
            mBegin(begin),                                                                                \
            mEnd(end)                                                                                     \
        {}                                                                                                \
                                                                                                          \
        Name##Iterator begin() {                                                                          \
            return mBegin;                                                                                \
        }                                                                                                 \
        Name##Iterator end() {                                                                            \
            return mEnd;                                                                                  \
        }                                                                                                 \
    };
    
    HALF_EDGE_MESH_RANGE(Vertex);
    HALF_EDGE_MESH_RANGE(Edge);
    HALF_EDGE_MESH_RANGE(HalfEdge);
    HALF_EDGE_MESH_RANGE(Face);
    
    HALF_EDGE_MESH_RANGE(VertexVertex);
    HALF_EDGE_MESH_RANGE(VertexOutHalfEdge);
    HALF_EDGE_MESH_RANGE(VertexInHalfEdge);
    HALF_EDGE_MESH_RANGE(VertexEdge);
    HALF_EDGE_MESH_RANGE(VertexFace);
    
    HALF_EDGE_MESH_RANGE(FaceVertex);
    HALF_EDGE_MESH_RANGE(FaceHalfEdge);
    HALF_EDGE_MESH_RANGE(FaceEdge);
    HALF_EDGE_MESH_RANGE(FaceFace);
    
    
    HALF_EDGE_MESH_RANGE(CVertex);
    HALF_EDGE_MESH_RANGE(CEdge);
    HALF_EDGE_MESH_RANGE(CHalfEdge);
    HALF_EDGE_MESH_RANGE(CFace);
    
    HALF_EDGE_MESH_RANGE(CVertexVertex);
    HALF_EDGE_MESH_RANGE(CVertexOutHalfEdge);
    HALF_EDGE_MESH_RANGE(CVertexInHalfEdge);
    HALF_EDGE_MESH_RANGE(CVertexEdge);
    HALF_EDGE_MESH_RANGE(CVertexFace);
    
    HALF_EDGE_MESH_RANGE(CFaceVertex);
    HALF_EDGE_MESH_RANGE(CFaceHalfEdge);
    HALF_EDGE_MESH_RANGE(CFaceEdge);
    HALF_EDGE_MESH_RANGE(CFaceFace);
    
#undef HALF_EDGE_MESH_RANGE
}


#endif



#ifndef IN_IDE_PARSER
// include the implementation
#include "HalfEdgeMesh.impl.h"
#endif

