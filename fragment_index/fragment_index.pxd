cimport numpy as np

ctypedef np.uint64_t uint64_t
ctypedef np.uint32_t uint32_t
ctypedef np.float32_t float32_t
ctypedef np.uint8_t uint8_t


cpdef enum SeriesEnum:
    parent = 0
    b = 1
    y = 2
    c = 3
    z = 4


cdef struct fragment_t:
    float32_t mass
    SeriesEnum series
    uint64_t parent_id


cpdef enum SortingEnum:
    unsorted = 0
    by_mass = 1
    by_parent = 2


cdef struct fragment_list_t:
    fragment_t* v
    size_t used
    size_t size
    SortingEnum sort_type
    float32_t min_mass
    float32_t max_mass


cdef struct interval_t:
    size_t start
    size_t end


cdef struct fragment_index_t:
    fragment_list_t* bins
    fragment_list_t* parent_index
    size_t size
    int bins_per_dalton
    double max_fragment_size
    SortingEnum sort_type


cdef struct fragment_index_search_t:
    fragment_index_t* index
    double query
    double error_tolerance
    size_t low_bin
    size_t high_bin
    size_t current_bin
    interval_t position_range
    size_t position


cdef struct fragment_index_traverse_t:
    fragment_index_t* index
    size_t current_bin
    size_t position


# fragment_list_t methods
cdef int init_fragment_list(fragment_list_t* self, size_t size) nogil
cdef int free_fragment_list(fragment_list_t* self) nogil
cdef int fragment_list_append(fragment_list_t* self, fragment_t fragment) nogil
cdef void fragment_list_sort(fragment_list_t* self, SortingEnum sort_type) nogil
cdef double fragment_list_lowest_mass(fragment_list_t* self) nogil
cdef double fragment_list_highest_mass(fragment_list_t* self) nogil
cdef int fragment_list_binary_search(fragment_list_t* self, double query, double error_tolerance, interval_t* out, size_t low_hint=*, size_t high_hint=*) nogil

# fragment_index_t methods
cdef int init_fragment_index(fragment_index_t* self, int bins_per_dalton=*, double max_fragment_size=*) nogil
cdef int free_fragment_index(fragment_index_t* self) nogil
cpdef size_t total_bins_for_mass(int bins_per_dalton, double max_fragment_size) nogil
cdef size_t bin_for_mass(fragment_index_t* self, double mass) nogil
cdef void fragment_index_sort(fragment_index_t* self, SortingEnum sort_type) nogil
cdef int fragment_index_add_parent(fragment_index_t* self, double mass, uint64_t parent_id) nogil
cdef int fragment_index_parents_for(fragment_index_t* self, double mass, double error_tolerance, interval_t* out) nogil
cdef int fragment_index_parents_for_range(fragment_index_t* self, double low, double high, double error_tolerance, interval_t* out) nogil

# fragment_index_search_t methods
cdef bint fragment_index_search_has_next(fragment_index_search_t* self) nogil
cdef int fragment_index_search_next(fragment_index_search_t* self, fragment_t* fragment) nogil
cdef int fragment_index_search(fragment_index_t* self, double mass, double error_tolerance, fragment_index_search_t* iterator) nogil

# fragment_index_traverse_t methods
cdef int fragment_index_traverse(fragment_index_t* self, fragment_index_traverse_t* iterator) nogil
cdef bint fragment_index_traverse_has_next(fragment_index_traverse_t* self) nogil
cdef int fragment_index_traverse_next(fragment_index_traverse_t* self, fragment_t* fragment) nogil
cdef int fragment_index_traverse_seek(fragment_index_traverse_t* self, double query, double error_tolerance=*) nogil


cdef class FragmentList(object):
    cdef:
        fragment_list_t* fragments
        public bint owned

    @staticmethod
    cdef FragmentList _create(fragment_list_t* pointer)

    cdef void _init_list(self)
    cpdef clear(self)
    cpdef append(self, float32_t mass, SeriesEnum series, uint64_t parent_id)
    cpdef sort(self, SortingEnum sort_type=?)
    cpdef interval_t search(self, double query, double error_tolerance=*)


cdef class FragmentIndex(object):
    cdef:
        fragment_index_t* index
        public list bins
        public bint owned
        public FragmentList parent_index

    @staticmethod
    cdef FragmentIndex _create(fragment_index_t* pointer)

    cpdef _wrap_bins(self)
    cpdef clear(self, reinit=*)
    cpdef size_t bin_for(self, double mass)
    cpdef add(self, double mass, SeriesEnum series, uint64_t parent_id)
    cpdef add_parent(self, double mass, uint64_t parent_id)
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


cdef class FragmentIndexTraverseIterator(object):
    cdef:
        fragment_index_traverse_t* iterator
        public bint owned

    @staticmethod
    cdef FragmentIndexTraverseIterator _create(fragment_index_traverse_t* iterator)

    cpdef int seek(self, double query, double error_tolerance=*)