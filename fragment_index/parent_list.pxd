cdef struct parent_t:
    float32_t mass
    int_id_t id
    int_id_t parent_id
    uint16_t start_position
    uint16_t size


cdef struct parent_list_t:
    parent_t* v
    size_t used
    size_t size
    uint8_t sort_type


cdef int init_parent_list(parent_list_t* self, size_t size) nogil
cdef int free_parent_list(parent_list_t* self) nogil
cdef int parent_list_append(parent_list_t* self, parent_t parent) nogil
cdef void parent_list_sort(parent_list_t* self, SortingEnum sort_type) nogil
cdef int parent_list_binary_search(parent_list_t* self, double query, double error_tolerance,
                                     interval_t* out, size_t low_hint=*, size_t high_hint=*) nogil
cdef int parent_list_to_bytes(parent_list_t* self, char** output_buffer, size_t* buffer_size) nogil
cdef int parent_list_from_bytes(parent_list_t* self, char*, size_t buffer_size) nogil


cdef class ParentList(object):
    cdef:
        parent_list_t* members
        public bint owned

    @staticmethod
    cdef ParentList _create(parent_list_t* pointer)

    cdef void _init_list(self)
    cpdef clear(self)
    cpdef append(self, float32_t mass, int_id_t id, int_id_t parent_id, uint16_t start_position, uint16_t size)
    cpdef sort(self, SortingEnum sort_type=?)
    cpdef interval_t search(self, double query, double error_tolerance=*)
    cpdef bytearray to_bytes(self)