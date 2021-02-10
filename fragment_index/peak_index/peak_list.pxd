
cdef struct peak_t:
    float32_t mass
    float32_t intensity
    int16_t charge
    int_id_t scan_id


cdef struct peak_list_t:
    peak_t* v
    size_t used
    size_t size
    uint8_t sort_type
    float32_t min_mass
    float32_t max_mass


# peak_list_t methods
cdef int init_peak_list(peak_list_t* self, size_t size) nogil
cdef int free_peak_list(peak_list_t* self) nogil
cdef int peak_list_append(peak_list_t* self, peak_t peak) nogil
cdef void peak_list_sort(peak_list_t* self, SortingEnum sort_type) nogil
cdef double peak_list_lowest_mass(peak_list_t* self) nogil
cdef double peak_list_highest_mass(peak_list_t* self) nogil
cdef int peak_list_binary_search(peak_list_t* self, double query, double error_tolerance,
                                     interval_t* out, size_t low_hint=*, size_t high_hint=*) nogil
cdef int peak_list_to_bytes(peak_list_t* self, char** output_buffer, size_t* buffer_size) nogil
cdef int peak_list_from_bytes(peak_list_t* self, char*, size_t buffer_size) nogil


# cdef class Peak(object):
#     cdef:
#         peak_t peak

#     @staticmethod
#     cdef Peak _create(peak_t* peak)

cdef class PeakList(object):
    cdef:
        peak_list_t* peaks
        public bint owned

    @staticmethod
    cdef PeakList _create(peak_list_t* pointer)

    cdef void _init_list(self)
    cpdef clear(self)
    cpdef append(self, float32_t mass, float32_t intensity, int16_t charge, int_id_t parent_id=*)
    cpdef sort(self, SortingEnum sort_type=?)
    cpdef interval_t search(self, double query, double error_tolerance=*)
    cpdef bytearray to_bytes(self)

