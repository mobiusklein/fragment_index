cimport numpy as np

ctypedef np.uint64_t uint64_t
ctypedef np.uint32_t uint32_t

ctypedef np.float32_t float32_t
ctypedef np.float64_t float64_t

ctypedef np.int32_t int32_t
ctypedef np.int64_t int64_t

ctypedef np.uint8_t uint8_t
ctypedef np.uint16_t uint16_t

ctypedef np.uint32_t int_id_t


cpdef enum SeriesEnum:
    parent = 0
    b = 1
    y = 2
    c = 3
    z = 4
    pY = 5


cpdef enum SortingEnum:
    unsorted = 0
    by_mass = 1
    by_parent = 2


cdef struct interval_t:
    size_t start
    size_t end


include "fragment_list.pxd"
include "parent_list.pxd"


cdef struct fragment_index_t:
    fragment_list_t* bins
    parent_list_t* parent_index
    size_t size
    uint32_t bins_per_dalton
    float64_t max_fragment_size
    uint8_t sort_type


cdef struct fragment_index_search_t:
    fragment_index_t* index
    double query
    double error_tolerance
    size_t low_bin
    size_t high_bin
    size_t current_bin
    interval_t position_range
    size_t position
    interval_t parent_id_interval


cdef struct fragment_index_traverse_t:
    fragment_index_t* index
    size_t current_bin
    size_t position
    interval_t parent_id_interval


# interval_t methods
cdef bint interval_contains(interval_t* self, size_t i) nogil
cdef bint interval_is_empty(interval_t* self) nogil
cdef bint interval_eq(interval_t* self, interval_t* other) nogil

# fragment_index_t methods
cdef int init_fragment_index(fragment_index_t* self, uint32_t bins_per_dalton=*, float64_t max_fragment_size=*) nogil
cdef int free_fragment_index(fragment_index_t* self) nogil
cpdef size_t total_bins_for_mass(int bins_per_dalton, double max_fragment_size) nogil
cdef size_t bin_for_mass(fragment_index_t* self, double mass) nogil
cdef void fragment_index_sort(fragment_index_t* self, SortingEnum sort_type) nogil
cdef int fragment_index_add_parent(fragment_index_t* self, double mass, int_id_t id, int_id_t parent_id=*, uint16_t start_position=*, uint16_t size=*) nogil
cdef int fragment_index_parents_for(fragment_index_t* self, double mass, double error_tolerance, interval_t* out) nogil
cdef int fragment_index_parents_for_range(fragment_index_t* self, double low, double high, double error_tolerance, interval_t* out) nogil

# fragment_index_search_t methods
cdef bint fragment_index_search_has_next(fragment_index_search_t* self) nogil
cdef int fragment_index_search_next(fragment_index_search_t* self, fragment_t* fragment) nogil
cdef int fragment_index_search(fragment_index_t* self, double mass, double error_tolerance,
                               fragment_index_search_t* iterator, interval_t parent_id_interval=*) nogil
cdef int fragment_index_search_set_parent_interval(fragment_index_search_t* self, interval_t parent_id_interval) nogil

# fragment_index_traverse_t methods
cdef int fragment_index_traverse(fragment_index_t* self, fragment_index_traverse_t* iterator, interval_t parent_id_interval=*) nogil
cdef bint fragment_index_traverse_has_next(fragment_index_traverse_t* self) nogil
cdef int fragment_index_traverse_next(fragment_index_traverse_t* self, fragment_t* fragment) nogil
cdef int fragment_index_traverse_seek(fragment_index_traverse_t* self, double query, double error_tolerance=*) nogil


cdef class FragmentIndex(object):
    cdef:
        fragment_index_t* index
        public list bins
        public bint owned
        public ParentList parent_index

    @staticmethod
    cdef FragmentIndex _create(fragment_index_t* pointer)

    cpdef _wrap_bins(self)
    cpdef clear(self, reinit=*)
    cpdef size_t bin_for(self, double mass)
    cpdef add(self, double mass, SeriesEnum series, int_id_t parent_id, uint16_t ordinal=*)
    cpdef add_parent(self, double mass, int_id_t id, int_id_t parent_id=*, uint16_t start_position=*, uint16_t size=*)
    cpdef sort(self, SortingEnum sort_type=?)
    cpdef size_t count(self)
    cpdef interval_t parents_for(self, double mass, double error_tolerance=*)
    cpdef interval_t parents_for_range(self, double low, double high, double error_tolerance=*)
    cpdef FragmentIndexSearchIterator search(self, double mass, double error_tolerance=*)
    cpdef FragmentIndexTraverseIterator traverse(self)


cdef class FragmentIndexSearchIterator(object):
    cdef:
        fragment_index_search_t* iterator
        public bint owned

    @staticmethod
    cdef FragmentIndexSearchIterator _create(fragment_index_search_t* iterator)
    cpdef bint set_parent_id_range(self, start, end)
    cpdef FragmentList all(self)

cdef class FragmentIndexTraverseIterator(object):
    cdef:
        fragment_index_traverse_t* iterator
        public bint owned

    @staticmethod
    cdef FragmentIndexTraverseIterator _create(fragment_index_traverse_t* iterator)

    cpdef int seek(self, double query, double error_tolerance=*)