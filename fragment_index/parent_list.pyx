

cdef int init_parent_list(parent_list_t* self, size_t size) nogil:
    self.v = <parent_t*>malloc(sizeof(parent_t) * size)
    self.used = 0
    self.size = size
    self.sort_type = SortingEnum.unsorted
    if self.v == NULL:
        return 1
    return 0


cdef int free_parent_list(parent_list_t* self) nogil:
    free(self.v)
    return 0


cdef int parent_list_append(parent_list_t* self, parent_t parent) nogil:
    if self.used >= self.size - 1:
        self.v = <parent_t*>realloc(self.v, sizeof(parent_t) * self.size * 2)
        if self.v == NULL:
            return 1
        self.size = self.size * 2
    self.v[self.used] = parent
    self.used += 1
    self.sort_type = SortingEnum.unsorted
    return 0


cdef int compare_by_mass_parent_t(const void * a, const void * b) nogil:
    if (<parent_t*>a).mass < (<parent_t*>b).mass:
        return -1
    elif (<parent_t*>a).mass == (<parent_t*>b).mass:
        return 0
    elif (<parent_t*>a).mass > (<parent_t*>b).mass:
        return 1


cdef int compare_by_parent_id_parent_t(const void * a, const void * b) nogil:
    if (<parent_t*>a).parent_id < (<parent_t*>b).parent_id:
        return -1
    elif (<parent_t*>a).parent_id == (<parent_t*>b).parent_id:
        return 0
    elif (<parent_t*>a).parent_id > (<parent_t*>b).parent_id:
        return 1


cdef void parent_list_sort(parent_list_t* self, SortingEnum sort_type) nogil:
    if sort_type == SortingEnum.by_mass:
        qsort(self.v, self.used, sizeof(parent_t), compare_by_mass_parent_t)
    elif sort_type == SortingEnum.by_parent:
        qsort(self.v, self.used, sizeof(parent_t), compare_by_parent_id_parent_t)
    self.sort_type = sort_type



@cython.cdivision(True)
cdef int parent_list_binary_search(parent_list_t* self, double query, double error_tolerance, interval_t* out, size_t low_hint=0, size_t high_hint=-1) nogil:
    cdef:
        size_t i, n, lo, hi, mid
        float x, err
    n = self.used
    lo = low_hint
    if high_hint == -1:
        hi = n
    else:
        hi = high_hint
    out.start = lo
    out.end = hi
    while hi != lo:
        mid = (hi + lo) // 2
        x = self.v[mid].mass
        err = (x - query) / query
        if lo == hi - 1 or fabs(err) < error_tolerance:
            i = mid
            while i > 0:
                x = self.v[i].mass
                if fabs(x - query) / query > error_tolerance:
                    break
                i -= 1
            out.start = i
            i = mid
            while i < n:
                x = self.v[i].mass
                if fabs(x - query) / query > error_tolerance:
                    break
                i += 1
            out.end = min(i + 1, self.used)
            return 0
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return 1


cdef int parent_list_to_bytes(parent_list_t* self, char** output_buffer, size_t* buffer_size) nogil:
    cdef:
        char* byte_buffer
        size_t i, n
        size_t offset
    offset = sizeof(size_t) + sizeof(uint8_t)
    n = sizeof(parent_t) * self.used + 1
    byte_buffer = <char*>malloc(offset + sizeof(char) * n + 1)
    if byte_buffer == NULL:
        return 1
    for i in range(offset + sizeof(char) * n + 1):
        byte_buffer[i] = '\0'
    memcpy(<void*>byte_buffer, <void*>&self.used, sizeof(size_t))
    memcpy(<void*>&byte_buffer[sizeof(size_t)], <void*>&self.sort_type, sizeof(uint8_t))
    memcpy(<void*>&byte_buffer[offset], <void*>self.v, n)
    output_buffer[0] = byte_buffer
    buffer_size[0] = n + offset
    return 0


cdef int parent_list_from_bytes(parent_list_t* self, char* buff, size_t buffer_size) nogil:
    cdef:
        int code
        size_t offset
        size_t i, n, list_size
        uint8_t sort_type
        parent_t f
    if self.v != NULL:
        free(self.v)
        self.v = NULL
        self.used = 0
        self.size = 0
    offset = 0
    if offset + sizeof(size_t) > buffer_size:
        return 3
    list_size = 0
    list_size = (<size_t*>buff)[0]
    offset += sizeof(size_t)
    code = init_parent_list(self, list_size)
    if code == 1:
        return 1
    sort_type = 0
    sort_type = (<uint8_t*>&buff[offset])[0]
    self.sort_type = sort_type
    offset += sizeof(uint8_t)
    for i in range(list_size):
        if offset + sizeof(parent_t) > buffer_size:
            return 3
        f = (<parent_t*>&buff[offset])[0]
        offset += sizeof(parent_t)
        code = parent_list_append(self, f)
        if code != 0:
            return 2
    return 0


cdef class ParentList(object):

    @staticmethod
    cdef ParentList _create(parent_list_t* pointer):
        cdef ParentList self = ParentList.__new__(ParentList)
        self.members = pointer
        self.owned = False
        return self

    @property
    def allocated(self):
        return self.members.size

    @property
    def sort_type(self):
        return self.members.sort_type

    def __init__(self, *args, **kwargs):
        self._init_list()

    cdef void _init_list(self):
        self.members = <parent_list_t*>malloc(sizeof(parent_list_t))
        self.owned = True
        init_parent_list(self.members, 32)

    cpdef clear(self):
        free_parent_list(self.members)
        free(self.members)
        self._init_list()

    def __dealloc__(self):
        if self.owned:
            free_parent_list(self.members)
            free(self.members)

    def __len__(self):
        return self.members.used

    def __getitem__(self, i):
        if isinstance(i, slice):
            out = []
            for j in range(i.start or 0, min(i.stop or len(self), len(self)), i.step or 1):
                out.append(self[j])
            return out
        if i  >= self.members.used:
            raise IndexError(i)
        elif i < 0:
            j = len(self) + i
            if j < 0:
                raise IndexError(i)
            i = j
        return self.members.v[i]

    def __iter__(self):
        for i in range(self.members.used):
            yield self.members.v[i]

    def __repr__(self):
        return "{self.__class__.__name__}({size})".format(self=self, size=len(self))

    cpdef append(self, float32_t mass, int_id_t id, int_id_t parent_id, uint16_t start_position, uint16_t size):
        cdef parent_t parent = parent_t(mass, id, parent_id, start_position, size)
        out = parent_list_append(self.members, parent)
        if out == 1:
            raise MemoryError()

    cpdef sort(self, SortingEnum sort_type=SortingEnum.by_mass):
        parent_list_sort(self.members, sort_type)

    cpdef interval_t search(self, double query, double error_tolerance=1e-5):
        cdef:
            int result
            interval_t interval_out
            double low, high
        result = parent_list_binary_search(self.members, query, error_tolerance, &interval_out)
        return interval_out

    cpdef bytearray to_bytes(self):
        cdef:
            char* buff
            size_t buffer_size
            int code
            bytearray out
        buffer_size = 0
        code = parent_list_to_bytes(self.members, &buff, &buffer_size)
        if code == 1:
            raise MemoryError()
        out = PyByteArray_FromStringAndSize(buff, buffer_size)
        free(buff)
        return out

    @classmethod
    def from_bytes(cls, bytearray bin_data):
        cdef:
            ParentList self
            int code
            size_t buffer_size
            char* buff
        self = ParentList()
        buffer_size = PyByteArray_Size(bin_data)
        buff = PyByteArray_AsString(bin_data)
        code = parent_list_from_bytes(self.members, buff, buffer_size)
        if code == 1:
            raise MemoryError()
        if code == 2:
            raise ValueError("Error adding members to list")
        if code == 3:
            raise ValueError("Malformed bytestring")
        return self

