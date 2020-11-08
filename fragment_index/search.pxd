cimport numpy as np

ctypedef np.uint64_t uint64_t
ctypedef np.uint32_t uint32_t
ctypedef np.float32_t float32_t
ctypedef np.uint8_t uint8_t
ctypedef np.uint16_t uint16_t


from  fragment_index.fragment_index cimport (
    fragment_index_t,
    fragment_list_t,
    fragment_t,
    fragment_index_search_t,
    interval_t,
    FragmentIndex)


cdef struct peak_t:
    float32_t mass
    float32_t intensity
    int charge


cdef struct peak_list_t:
    peak_t* v
    size_t size
    size_t used


cdef struct match_t:
    uint32_t parent_id
    float32_t score
    uint32_t hit_count


cdef struct match_list_t:
    match_t* v
    size_t size
    size_t used


cdef struct fragment_search_t:
    fragment_index_t* index
    peak_list_t* peak_list
    match_list_t* match_list
    interval_t parent_interval


cdef int init_peak_list(peak_list_t* self, size_t size) nogil
cdef int free_peak_list(peak_list_t* self) nogil
cdef int peak_list_append(peak_list_t* self, peak_t peak) nogil


cdef int init_match_list(match_list_t* self, size_t size) nogil
cdef int free_match_list(match_list_t* self) nogil
cdef int match_list_append(match_list_t* self, match_t match) nogil


cdef int search_fragment_index(fragment_index_t* index, peak_list_t* peak_list, double precursor_mass, double parent_error_low,
                               double parent_error_high, double error_tolerance, fragment_search_t* result) nogil


cdef class PeakList(object):
    cdef:
        peak_list_t* peaks
        public bint owned

    @staticmethod
    cdef PeakList _create(peak_list_t* pointer)

    cdef void _init_list(self)
    cpdef clear(self)
    cpdef append(self, float32_t mass, float32_t intensity, int charge)


cdef class MatchList(object):
    cdef:
        match_list_t* matches
        public bint owned

    @staticmethod
    cdef MatchList _create(match_list_t* pointer)

    cdef void _init_list(self)
    cpdef clear(self)
    cpdef append(self, uint32_t parent_id, float32_t score, uint32_t hit_count)


cdef class SearchResult(object):
    cdef:
        public FragmentIndex index
        public MatchList matches
        public PeakList peak_list
        fragment_search_t* search
