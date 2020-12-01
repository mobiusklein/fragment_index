cimport numpy as np

ctypedef np.uint64_t uint64_t
ctypedef np.uint32_t uint32_t

ctypedef np.float32_t float32_t

ctypedef np.int16_t int16_t
ctypedef np.int32_t int32_t
ctypedef np.int64_t int64_t

ctypedef np.uint8_t uint8_t
ctypedef np.uint16_t uint16_t

ctypedef np.uint32_t int_id_t


cpdef enum SortingEnum:
    unsorted = 0
    by_mass = 1
    by_parent = 2


cdef struct interval_t:
    size_t start
    size_t end

include "peak_list.pxd"
include "parent_list.pxd"

cdef struct peak_index_t:
    peak_list_t* bins
    parent_list_t* parent_index
    size_t size
    int bins_per_dalton
    double max_peak_size
    uint8_t sort_type


cdef struct peak_index_search_t:
    peak_index_t* index
    double query
    double error_tolerance
    size_t low_bin
    size_t high_bin
    size_t current_bin
    interval_t position_range
    size_t position
    interval_t scan_id_interval


cdef struct peak_index_traverse_t:
    peak_index_t* index
    size_t current_bin
    size_t position
    interval_t scan_id_interval


# interval_t methods
cdef bint interval_contains(interval_t* self, size_t i) nogil
cdef bint interval_is_empty(interval_t* self) nogil
cdef bint interval_eq(interval_t* self, interval_t* other) nogil

# peak_index_t methods
cdef int init_peak_index(peak_index_t* self, int bins_per_dalton=*, double max_peak_size=*) nogil
cdef int free_peak_index(peak_index_t* self) nogil
cpdef size_t total_bins_for_mass(int bins_per_dalton, double max_peak_size) nogil
cdef size_t bin_for_mass(peak_index_t* self, double mass) nogil
cdef void peak_index_sort(peak_index_t* self, SortingEnum sort_type) nogil
cdef int peak_index_add_parent(peak_index_t* self, double mass, int_id_t id, int_id_t scan_id=*, uint16_t start_position=*, uint16_t size=*) nogil
cdef int peak_index_parents_for(peak_index_t* self, double mass, double error_tolerance, interval_t* out) nogil
cdef int peak_index_parents_for_range(peak_index_t* self, double low, double high, double error_tolerance, interval_t* out) nogil

# peak_index_search_t methods
cdef bint peak_index_search_has_next(peak_index_search_t* self) nogil
cdef int peak_index_search_next(peak_index_search_t* self, peak_t* fragment) nogil
cdef int peak_index_search(peak_index_t* self, double mass, double error_tolerance,
                               peak_index_search_t* iterator, interval_t scan_id_interval=*) nogil
cdef int peak_index_search_set_parent_interval(peak_index_search_t* self, interval_t scan_id_interval) nogil

# peak_index_traverse_t methods
cdef int peak_index_traverse(peak_index_t* self, peak_index_traverse_t* iterator, interval_t scan_id_interval=*) nogil
cdef bint peak_index_traverse_has_next(peak_index_traverse_t* self) nogil
cdef int peak_index_traverse_next(peak_index_traverse_t* self, peak_t* fragment) nogil
cdef int peak_index_traverse_seek(peak_index_traverse_t* self, double query, double error_tolerance=*) nogil



cdef class PeakIndex(object):
    cdef:
        peak_index_t* index
        public list bins
        public bint owned
        public ParentList parent_index

    @staticmethod
    cdef PeakIndex _create(peak_index_t* pointer)

    cpdef _wrap_bins(self)
    cpdef clear(self, reinit=*)
    cpdef size_t bin_for(self, double mass)
    cpdef add(self, float32_t mass, float32_t intensity, int16_t charge, int_id_t scan_id=*)
    cpdef add_parent(self, double mass, int_id_t id, int_id_t scan_id=*, uint16_t start_position=*, uint16_t size=*)
    cpdef sort(self, SortingEnum sort_type=?)
    cpdef size_t count(self)
    cpdef interval_t parents_for(self, double mass, double error_tolerance=*)
    cpdef interval_t parents_for_range(self, double low, double high, double error_tolerance=*)
    cpdef PeakIndexSearchIterator search(self, double mass, double error_tolerance=*)
    cpdef PeakIndexTraverseIterator traverse(self)


cdef class PeakIndexSearchIterator(object):
    cdef:
        peak_index_search_t* iterator
        public bint owned

    @staticmethod
    cdef PeakIndexSearchIterator _create(peak_index_search_t* iterator)
    cpdef bint set_scan_id_range(self, start, end)
    cpdef PeakList all(self)

cdef class PeakIndexTraverseIterator(object):
    cdef:
        peak_index_traverse_t* iterator
        public bint owned

    @staticmethod
    cdef PeakIndexTraverseIterator _create(peak_index_traverse_t* iterator)

    cpdef int seek(self, double query, double error_tolerance=*)