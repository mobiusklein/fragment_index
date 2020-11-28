fragment_t_size = sizeof(fragment_t)


cdef int compare_by_mass(const void * a, const void * b) nogil:
    if (<fragment_t*>a).mass < (<fragment_t*>b).mass:
        return -1
    elif (<fragment_t*>a).mass == (<fragment_t*>b).mass:
        return 0
    elif (<fragment_t*>a).mass > (<fragment_t*>b).mass:
        return 1


cdef int compare_by_parent_id(const void * a, const void * b) nogil:
    if (<fragment_t*>a).parent_id < (<fragment_t*>b).parent_id:
        return -1
    elif (<fragment_t*>a).parent_id == (<fragment_t*>b).parent_id:
        return 0
    elif (<fragment_t*>a).parent_id > (<fragment_t*>b).parent_id:
        return 1


cdef int init_fragment_list(fragment_list_t* self, size_t size) nogil:
    self.v = <fragment_t*>malloc(sizeof(fragment_t) * size)
    self.used = 0
    self.size = size
    self.min_mass = self.max_mass = -1
    self.sort_type = SortingEnum.unsorted
    if self.v == NULL:
        return 1
    return 0


cdef int free_fragment_list(fragment_list_t* self) nogil:
    free(self.v)
    return 0


cdef int fragment_list_append(fragment_list_t* self, fragment_t fragment) nogil:
    if self.used >= self.size - 1:
        self.v = <fragment_t*>realloc(self.v, sizeof(fragment_t) * self.size * 2)
        if self.v == NULL:
            return 1
        self.size = self.size * 2
    self.v[self.used] = fragment
    self.used += 1
    # self.min_mass = self.max_mass = -1
    # self.sort_type = SortingEnum.unsorted
    return 0


cdef void fragment_list_sort(fragment_list_t* self, SortingEnum sort_type) nogil:
    if sort_type == SortingEnum.by_mass:
        qsort(self.v, self.used, sizeof(fragment_t), compare_by_mass)
    elif sort_type == SortingEnum.by_parent:
        qsort(self.v, self.used, sizeof(fragment_t), compare_by_parent_id)
    self.sort_type = sort_type
    # Update boundaries
    fragment_list_lowest_mass(self)
    fragment_list_highest_mass(self)

cdef double fragment_list_lowest_mass(fragment_list_t* self) nogil:
    cdef:
        size_t i, n
        double mass
    if self.used == 0:
        return 0
    if self.min_mass != -1:
        return self.min_mass
    if self.sort_type == SortingEnum.by_mass:
        self.min_mass = self.v[0].mass
    else:
        mass = INF
        for i in range(self.used):
            mass = min(self.v[i].mass, mass)
        self.min_mass = mass
    return self.min_mass


cdef double fragment_list_highest_mass(fragment_list_t* self) nogil:
    cdef:
        size_t i
        double mass
    if self.used == 0:
        return 0
    if self.max_mass != -1:
        return self.max_mass
    if self.sort_type == SortingEnum.by_mass:
        self.max_mass = self.v[self.used - 1].mass
    else:
        mass = 0
        for i in range(self.used):
            mass = max(self.v[i].mass, mass)
        self.max_mass = mass
    return self.max_mass


@cython.cdivision(True)
cdef int fragment_list_binary_search(fragment_list_t* self, double query, double error_tolerance, interval_t* out, size_t low_hint=0, size_t high_hint=-1) nogil:
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


cdef int fragment_list_to_bytes(fragment_list_t* self, char** output_buffer, size_t* buffer_size) nogil:
    cdef:
        char* byte_buffer
        size_t i, n
        size_t offset
    offset = sizeof(size_t) + sizeof(uint8_t)
    n = sizeof(fragment_t) * self.used + 1
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


cdef int fragment_list_from_bytes(fragment_list_t* self, char* buff, size_t buffer_size) nogil:
    cdef:
        int code
        size_t offset
        size_t i, n, list_size
        uint8_t sort_type
        fragment_t f
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
    code = init_fragment_list(self, list_size)
    if code == 1:
        return 1
    sort_type = 0
    sort_type = (<uint8_t*>&buff[offset])[0]
    self.sort_type = sort_type
    offset += sizeof(uint8_t)
    for i in range(list_size):
        if offset + sizeof(fragment_t) > buffer_size:
            return 3
        f = (<fragment_t*>&buff[offset])[0]
        offset += sizeof(fragment_t)
        code = fragment_list_append(self, f)
        if code != 0:
            return 2
    return 0


cdef class FragmentList(object):

    @staticmethod
    cdef FragmentList _create(fragment_list_t* pointer):
        cdef FragmentList self = FragmentList.__new__(FragmentList)
        self.fragments = pointer
        self.owned = False
        return self

    @property
    def allocated(self):
        return self.fragments.size

    @property
    def sort_type(self):
        return self.fragments.sort_type

    def __init__(self, *args, **kwargs):
        self._init_list()

    cdef void _init_list(self):
        self.fragments = <fragment_list_t*>malloc(sizeof(fragment_list_t))
        self.owned = True
        init_fragment_list(self.fragments, 32)

    cpdef clear(self):
        free_fragment_list(self.fragments)
        free(self.fragments)
        self._init_list()

    def __dealloc__(self):
        if self.owned:
            free_fragment_list(self.fragments)
            free(self.fragments)

    def __len__(self):
        return self.fragments.used

    def __getitem__(self, i):
        if isinstance(i, slice):
            out = []
            for j in range(i.start or 0, min(i.stop or len(self), len(self)), i.step or 1):
                out.append(self[j])
            return out
        if i  >= self.fragments.used:
            raise IndexError(i)
        elif i < 0:
            j = len(self) + i
            if j < 0:
                raise IndexError(i)
            i = j
        return self.fragments.v[i]

    def __iter__(self):
        for i in range(self.fragments.used):
            yield self.fragments.v[i]

    def __repr__(self):
        return "{self.__class__.__name__}({size})".format(self=self, size=len(self))

    cpdef append(self, float32_t mass, SeriesEnum series, int_id_t parent_id, uint16_t ordinal=0):
        cdef fragment_t fragment = fragment_t(mass, parent_id, series, ordinal)
        out = fragment_list_append(self.fragments, fragment)
        if out == 1:
            raise MemoryError()

    cpdef sort(self, SortingEnum sort_type=SortingEnum.by_mass):
        fragment_list_sort(self.fragments, sort_type)

    cpdef interval_t search(self, double query, double error_tolerance=1e-5):
        cdef:
            int result
            interval_t interval_out
            double low, high
        result = fragment_list_binary_search(self.fragments, query, error_tolerance, &interval_out)
        return interval_out

    @property
    def lowest_mass(self):
        return fragment_list_lowest_mass(self.fragments)

    @property
    def highest_mass(self):
        return fragment_list_highest_mass(self.fragments)

    cpdef bytearray to_bytes(self):
        cdef:
            char* buff
            size_t buffer_size
            int code
            bytearray out
        buffer_size = 0
        code = fragment_list_to_bytes(self.fragments, &buff, &buffer_size)
        if code == 1:
            raise MemoryError()
        out = PyByteArray_FromStringAndSize(buff, buffer_size)
        free(buff)
        return out

    @classmethod
    def from_bytes(cls, bytearray bin_data):
        cdef:
            FragmentList self
            int code
            size_t buffer_size
            char* buff
        self = FragmentList()
        buffer_size = PyByteArray_Size(bin_data)
        buff = PyByteArray_AsString(bin_data)
        code = fragment_list_from_bytes(self.fragments, buff, buffer_size)
        if code == 1:
            raise MemoryError()
        if code == 2:
            raise ValueError("Error adding fragments to list")
        if code == 3:
            raise ValueError("Malformed bytestring")
        return self


@cython.final
@cython.freelist(10000)
cdef class Fragment(object):
    @staticmethod
    cdef Fragment _create(fragment_t* fragment):
        cdef Fragment self = Fragment.__new__(Fragment)
        self.fragment = fragment[0]

    def __init__(self, *args, **kwargs):
        raise NotImplementedError()

    @property
    def mass(self):
        return self.fragment.mass

    @property
    def series(self):
        return  SeriesEnum[self.fragment.series]

    @property
    def ordinal(self):
        return self.fragment.ordinal

    @property
    def parent_id(self):
        return self.fragment.parent_id

    def __repr__(self):
        return "{self.__class__.__name__}({self.mass}, {self.series}, {self.parent_id})".format(self=self)

    def __getitem__(self, key):
        if key == 'mass':
            return self.fragment.mass
        if key == 'series':
            return self.fragment.series
        if key == 'parent_id':
            return self.fragment.parent_id
        if key == 'ordinal':
            return self.fragment.ordinal
        raise KeyError(key)