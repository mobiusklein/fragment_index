cimport cython
from libc.stdlib cimport malloc, realloc, calloc, free, qsort
from libc.string cimport memcpy
from libc.math cimport floor, fabs
from libc.stdio cimport FILE, fread, fwrite, fopen, fclose, fprintf

from cpython cimport Py_buffer
from cpython.string cimport PyString_AsString, PyString_AsStringAndSize, PyString_FromStringAndSize
from cpython.bytearray cimport PyByteArray_FromStringAndSize, PyByteArray_Size, PyByteArray_AsString

cdef extern from * nogil:
    int printf (const char *template, ...)


MAX_SIZE = (<size_t>-1) - sizeof(size_t)


cdef int init_serializer_block_list(serializer_block_list_t* self, size_t size) nogil:
    self.v = <serializer_block_t*>malloc(sizeof(serializer_block_t) * size)
    if self.v == NULL:
        return 1
    self.used = 0
    self.size = size
    return 0


cdef int release_serializer_block(serializer_block_t* self) nogil:
    if self.data != NULL:
        free(self.data)
        self.data = NULL
    return 0


cdef int free_serializer_block_list(serializer_block_list_t* self) nogil:
    cdef:
        size_t i, n

    n = self.used
    for i in range(n):
        release_serializer_block(&self.v[i])
    free(self.v)
    return 0


cdef int serializer_block_list_append(serializer_block_list_t* self, serializer_block_t block) nogil:
    if self.used >= self.size - 1:
        self.v = <serializer_block_t*>realloc(self.v, sizeof(serializer_block_t) * self.size * 2)
        if self.v == NULL:
            return 1
        self.size = self.size * 2
    self.v[self.used] = block
    self.used += 1
    return 0


cdef int serializer_block_list_to_file(serializer_block_list_t* self, FILE* file_handle) nogil:
    cdef:
        size_t i, n
        char* byte_buffer
    fwrite(<void*>&self.used, sizeof(size_t), 1, file_handle)
    for i in range(self.used):
        fwrite(<void*>&self.v[i].size, sizeof(size_t), 1, file_handle)
        fwrite(<void*>self.v[i].data, sizeof(char), self.v[i].size, file_handle)
    return 0


cdef Py_ssize_t[1] BLOCK_STRIDE = [sizeof(char)]

cdef serializer_block_t serializer_block_from_pybytes(bytes data):
    cdef:
        size_t n
        char* strdata
        char* copied
        serializer_block_t block

    n = len(data)
    if n > MAX_SIZE:
        raise ValueError("The data exceeds the size of a single block")
    copied = <char*>malloc(sizeof(char) * n)
    strdata = PyString_AsString(data)
    memcpy(copied, strdata, n)
    block.data = copied
    block.size = n
    return block


cdef class SerializerBlock(object):
    cdef:
        serializer_block_t* block
        public bint owned
        Py_ssize_t[1] shape

    @staticmethod
    cdef SerializerBlock _create(serializer_block_t* block):
        cdef SerializerBlock self = SerializerBlock.__new__(SerializerBlock)
        self.block = block
        self.owned = False
        self.shape[0] = block.size
        return self

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        itemsize = sizeof(char)

        buffer.buf = self.block.data
        buffer.len = self.block.size * itemsize # itemsize is 1
        buffer.format = 'c'
        buffer.itemsize = itemsize
        buffer.internal = NULL
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = False
        buffer.shape = self.shape
        buffer.strides = BLOCK_STRIDE
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

    def __init__(self, bytes data=None):
        if data is not None:
            self._block_from_bytes(data)
        else:
            self.block = NULL

    def __dealloc__(self):
        if self.owned:
            if self.block != NULL:
                release_serializer_block(self.block)
                free(self.block)
            self.block = NULL

    def __getitem__(self, i):
        if self.block != NULL and  self.block.data != NULL:
            if isinstance(i, slice):
                start = i.start or 0
                stop = min(i.stop or self.block.size, self.block.size)
                step = i.step or 1
                return [self[j] for j in range(start, stop, step)]
            elif i >= self.block.size:
                raise IndexError(i)
            return self.block.data[i]
        else:
            raise ValueError("Cannot access data block after releasing memory")

    def __setitem__(self, i, value):
        if not isinstance(value, (bytes, bytearray)) or len(value) > 1:
            raise ValueError("Can only set values of a single byte")
        if self.block != NULL and  self.block.data != NULL:
            if i >= self.block.size:
                raise IndexError(i)
            self.block.data[i] = value

    def __len__(self):
        return self.block.size

    def _block_from_bytes(self, bytes data):
        cdef:
            serializer_block_t* block
        block = <serializer_block_t*>malloc(sizeof(serializer_block_t))
        block[0] = serializer_block_from_pybytes(data)
        self.block = block
        self.owned = True
        self.shape[0] = block.size


cdef class SerializerBlockList(object):
    cdef:
        serializer_block_list_t* blocks
        public bint owned

    @staticmethod
    cdef SerializerBlockList _create(serializer_block_list_t* pointer):
        cdef SerializerBlockList self = SerializerBlockList.__new__(SerializerBlockList)
        self.blocks = pointer
        self.owned = False
        return self

    @property
    def allocated(self):
        return self.blocks.size

    def __init__(self, *args, **kwargs):
        self._init_list()

    cdef void _init_list(self):
        self.blocks = <serializer_block_list_t*>malloc(sizeof(serializer_block_list_t))
        self.owned = True
        init_serializer_block_list(self.blocks, 4)

    cpdef clear(self):
        free_serializer_block_list(self.blocks)
        free(self.blocks)
        self._init_list()

    def __dealloc__(self):
        if self.owned:
            free_serializer_block_list(self.blocks)
            free(self.blocks)

    def __len__(self):
        return self.blocks.used

    def __getitem__(self, i):
        if isinstance(i, slice):
            out = []
            for j in range(i.start or 0, min(i.stop or len(self), len(self)), i.step or 1):
                out.append(self[j])
            return out
        if i  >= self.blocks.used:
            raise IndexError(i)
        elif i < 0:
            j = len(self) + i
            if j < 0:
                raise IndexError(i)
            i = j
        return SerializerBlock._create(&self.blocks.v[i])

    def __iter__(self):
        for i in range(self.blocks.used):
            yield self.blocks.v[i]

    def __repr__(self):
        return "{self.__class__.__name__}({size})".format(self=self, size=len(self))

    cpdef append(self, bytes data):
        cdef:
            serializer_block_t block
        block = serializer_block_from_pybytes(data)
        out = serializer_block_list_append(self.blocks, block)
        if out == 1:
            raise MemoryError()

    def write(self, file_handle):
        cdef:
            size_t i, n, offset
            size_t block_size
            char* byte_buffer
            serializer_block_t* block
            bytes payload

        n = self.blocks.used
        byte_buffer = <char*>malloc(sizeof(size_t))
        memcpy(<void*>byte_buffer, <void*>&n, sizeof(size_t))
        payload = PyString_FromStringAndSize(byte_buffer, sizeof(size_t))
        file_handle.write(payload)

        for i in range(self.blocks.used):
            block = &self.blocks.v[i]
            block_size = block.size
            memcpy(<void*>byte_buffer, <void*>&(block_size), sizeof(size_t))
            payload = PyString_FromStringAndSize(byte_buffer, sizeof(size_t))
            file_handle.write(payload)
            payload = PyString_FromStringAndSize(block.data, block.size)
            file_handle.write(payload)
        free(byte_buffer)

    @classmethod
    def read(cls, file_handle):
        cdef:
            bytes data
            char* buff
            size_t i, n_blocks, block_size
            SerializerBlockList self
            serializer_block_t block
        data = file_handle.read(sizeof(size_t))
        buff = data
        n_blocks = (<size_t*>buff)[0]

        self = cls()
        for i in range(n_blocks):
            print("Reading %d bytes to get block size" % sizeof(size_t))
            data = file_handle.read(sizeof(size_t))
            buff = data
            block_size = (<size_t*>buff)[0]
            print("Block Size: %d" % block_size)
            data = file_handle.read(block_size)
            print("Block Read: %r" % data)
            block = serializer_block_from_pybytes(data)
            serializer_block_list_append(self.blocks, block)
        return self

    def to_file(self, path):
        cdef:
            FILE* stream
        if hasattr(path, 'write'):
            return self.write(path)
        stream = fopen(path, 'w')
        n = serializer_block_list_to_file(self.blocks, stream)
        fclose(stream)
        return n