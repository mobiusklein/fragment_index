from libc.stdio cimport (FILE, )

# The maximum length of a block of data to be encoded as a string.
cdef size_t MAX_SIZE = (<size_t>-1) - sizeof(size_t)

cdef struct serializer_block_t:
    char* data
    size_t size

cdef struct serializer_block_list_t:
    serializer_block_t* v
    size_t used
    size_t size


cdef int init_serializer_block_list(serializer_block_list_t* self, size_t size) nogil
cdef int free_serializer_block_list(serializer_block_list_t* self) nogil
cdef int serializer_block_list_append(serializer_block_list_t* self, serializer_block_t block) nogil

cdef int serializer_block_list_to_file(serializer_block_list_t* self, FILE* file_handle) nogil
