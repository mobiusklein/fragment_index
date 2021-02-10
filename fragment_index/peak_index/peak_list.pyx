peak_t_size = sizeof(peak_t)


cdef int compare_by_mass(const void * a, const void * b) nogil:
    if (<peak_t*>a).mass < (<peak_t*>b).mass:
        return -1
    elif (<peak_t*>a).mass == (<peak_t*>b).mass:
        return 0
    elif (<peak_t*>a).mass > (<peak_t*>b).mass:
        return 1


cdef int compare_by_scan_id(const void * a, const void * b) nogil:
    if (<peak_t*>a).scan_id < (<peak_t*>b).scan_id:
        return -1
    elif (<peak_t*>a).scan_id == (<peak_t*>b).scan_id:
        return 0
    elif (<peak_t*>a).scan_id > (<peak_t*>b).scan_id:
        return 1


cdef int init_peak_list(peak_list_t* self, size_t size) nogil:
    self.v = <peak_t*>malloc(sizeof(peak_t) * size)
    self.used = 0
    self.size = size
    self.min_mass = self.max_mass = -1
    self.sort_type = SortingEnum.unsorted
    if self.v == NULL:
        return 1
    return 0


cdef int free_peak_list(peak_list_t* self) nogil:
    free(self.v)
    return 0


cdef int peak_list_append(peak_list_t* self, peak_t fragment) nogil:
    if self.used >= self.size - 1:
        self.v = <peak_t*>realloc(self.v, sizeof(peak_t) * self.size * 2)
        if self.v == NULL:
            return 1
        self.size = self.size * 2
    self.v[self.used] = fragment
    self.used += 1
    # self.min_mass = self.max_mass = -1
    # self.sort_type = SortingEnum.unsorted
    return 0


cdef void peak_list_sort(peak_list_t* self, SortingEnum sort_type) nogil:
    if sort_type == SortingEnum.by_mass:
        qsort(self.v, self.used, sizeof(peak_t), compare_by_mass)
    elif sort_type == SortingEnum.by_parent:
        qsort(self.v, self.used, sizeof(peak_t), compare_by_scan_id)
    self.sort_type = sort_type
    # Update boundaries
    peak_list_lowest_mass(self)
    peak_list_highest_mass(self)

cdef double peak_list_lowest_mass(peak_list_t* self) nogil:
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


cdef double peak_list_highest_mass(peak_list_t* self) nogil:
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
cdef int peak_list_binary_search(peak_list_t* self, double query, double error_tolerance, interval_t* out, size_t low_hint=0, size_t high_hint=-1) nogil:
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


cdef int peak_list_to_bytes(peak_list_t* self, char** output_buffer, size_t* buffer_size) nogil:
    cdef:
        char* byte_buffer
        size_t i, n
        size_t offset
    offset = sizeof(size_t) + sizeof(uint8_t)
    n = sizeof(peak_t) * self.used + 1
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


cdef int peak_list_from_bytes(peak_list_t* self, char* buff, size_t buffer_size) nogil:
    cdef:
        int code
        size_t offset
        size_t i, n, list_size
        uint8_t sort_type
        peak_t f
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
    code = init_peak_list(self, list_size)
    if code == 1:
        return 1
    sort_type = 0
    sort_type = (<uint8_t*>&buff[offset])[0]
    self.sort_type = sort_type
    offset += sizeof(uint8_t)
    for i in range(list_size):
        if offset + sizeof(peak_t) > buffer_size:
            return 3
        f = (<peak_t*>&buff[offset])[0]
        offset += sizeof(peak_t)
        code = peak_list_append(self, f)
        if code != 0:
            return 2
    return 0


cdef class PeakList(object):

    @staticmethod
    cdef PeakList _create(peak_list_t* pointer):
        cdef PeakList self = PeakList.__new__(PeakList)
        self.peaks = pointer
        self.owned = False
        return self

    @property
    def allocated(self):
        return self.peaks.size

    @property
    def sort_type(self):
        return self.peaks.sort_type

    def __init__(self, *args, **kwargs):
        self._init_list()

    cdef void _init_list(self):
        self.peaks = <peak_list_t*>malloc(sizeof(peak_list_t))
        self.owned = True
        init_peak_list(self.peaks, 32)

    cpdef clear(self):
        free_peak_list(self.peaks)
        free(self.peaks)
        self._init_list()

    def __dealloc__(self):
        if self.owned:
            free_peak_list(self.peaks)
            free(self.peaks)

    def __len__(self):
        return self.peaks.used

    def __getitem__(self, i):
        if isinstance(i, slice):
            out = []
            for j in range(i.start or 0, min(i.stop or len(self), len(self)), i.step or 1):
                out.append(self[j])
            return out
        if i  >= self.peaks.used:
            raise IndexError(i)
        elif i < 0:
            j = len(self) + i
            if j < 0:
                raise IndexError(i)
            i = j
        return self.peaks.v[i]

    def __iter__(self):
        for i in range(self.peaks.used):
            yield self.peaks.v[i]

    def __repr__(self):
        return "{self.__class__.__name__}({size})".format(self=self, size=len(self))

    cpdef append(self, float32_t mass, float32_t intensity, int16_t charge, int_id_t parent_id=0):
        cdef peak_t peak = peak_t(mass, intensity, charge, scan_id)
        out = peak_list_append(self.peaks, peak)
        if out == 1:
            raise MemoryError()

    cpdef sort(self, SortingEnum sort_type=SortingEnum.by_mass):
        peak_list_sort(self.peaks, sort_type)

    cpdef interval_t search(self, double query, double error_tolerance=1e-5):
        cdef:
            int result
            interval_t interval_out
            double low, high
        result = peak_list_binary_search(self.peaks, query, error_tolerance, &interval_out)
        return interval_out

    @property
    def lowest_mass(self):
        return peak_list_lowest_mass(self.peaks)

    @property
    def highest_mass(self):
        return peak_list_highest_mass(self.peaks)

    cpdef bytearray to_bytes(self):
        cdef:
            char* buff
            size_t buffer_size
            int code
            bytearray out
        buffer_size = 0
        code = peak_list_to_bytes(self.peaks, &buff, &buffer_size)
        if code == 1:
            raise MemoryError()
        out = PyByteArray_FromStringAndSize(buff, buffer_size)
        free(buff)
        return out

    @classmethod
    def from_bytes(cls, bytearray bin_data):
        cdef:
            PeakList self
            int code
            size_t buffer_size
            char* buff
        self = PeakList()
        buffer_size = PyByteArray_Size(bin_data)
        buff = PyByteArray_AsString(bin_data)
        code = peak_list_from_bytes(self.peaks, buff, buffer_size)
        if code == 1:
            raise MemoryError()
        if code == 2:
            raise ValueError("Error adding fragments to list")
        if code == 3:
            raise ValueError("Malformed bytestring")
        return self

