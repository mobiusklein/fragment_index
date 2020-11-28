
cdef struct fragment_t:
    float32_t mass
    int_id_t parent_id
    uint8_t series
    uint16_t ordinal


cdef struct fragment_list_t:
    fragment_t* v
    size_t used
    size_t size
    uint8_t sort_type
    float32_t min_mass
    float32_t max_mass


# fragment_list_t methods
cdef int init_fragment_list(fragment_list_t* self, size_t size) nogil
cdef int free_fragment_list(fragment_list_t* self) nogil
cdef int fragment_list_append(fragment_list_t* self, fragment_t fragment) nogil
cdef void fragment_list_sort(fragment_list_t* self, SortingEnum sort_type) nogil
cdef double fragment_list_lowest_mass(fragment_list_t* self) nogil
cdef double fragment_list_highest_mass(fragment_list_t* self) nogil
cdef int fragment_list_binary_search(fragment_list_t* self, double query, double error_tolerance,
                                     interval_t* out, size_t low_hint=*, size_t high_hint=*) nogil
cdef int fragment_list_to_bytes(fragment_list_t* self, char** output_buffer, size_t* buffer_size) nogil
cdef int fragment_list_from_bytes(fragment_list_t* self, char*, size_t buffer_size) nogil


cdef class Fragment(object):
    cdef:
        fragment_t fragment

    @staticmethod
    cdef Fragment _create(fragment_t* fragment)

cdef class FragmentList(object):
    cdef:
        fragment_list_t* fragments
        public bint owned

    @staticmethod
    cdef FragmentList _create(fragment_list_t* pointer)

    cdef void _init_list(self)
    cpdef clear(self)
    cpdef append(self, float32_t mass, SeriesEnum series, int_id_t parent_id, uint16_t ordinal=*)
    cpdef sort(self, SortingEnum sort_type=?)
    cpdef interval_t search(self, double query, double error_tolerance=*)
    cpdef bytearray to_bytes(self)

