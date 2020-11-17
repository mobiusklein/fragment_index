# cython: embedsignatures:True

cimport cython
from libc.stdlib cimport malloc, realloc, calloc, free, qsort
from libc.string cimport memcpy
from libc.math cimport floor, fabs

from cpython.bytearray cimport PyByteArray_FromStringAndSize, PyByteArray_Size, PyByteArray_AsString

cdef extern from * nogil:
    int printf (const char *template, ...)
    void qsort (void *base, unsigned short n, unsigned short w, int (*cmp_func)(void*, void*))


cdef double _round(double x) nogil:
    return floor(x + 0.5)


cdef double INF = float('inf')


cdef bint interval_contains(interval_t* self, size_t i) nogil:
    return self.start <= i < self.end


cdef bint interval_is_empty(interval_t* self) nogil:
    return self.start == self.end


cdef bint interval_eq(interval_t* self, interval_t* other) nogil:
    if other == NULL:
        return self == NULL
    return self.start == other.start and self.end == other.end


cdef interval_t OPEN_INTERVAL
OPEN_INTERVAL.start = 0
OPEN_INTERVAL.end = -1 # unsigned wrap-around to largest value here


cdef interval_t EMPTY_INTERVAL
EMPTY_INTERVAL.start = 0
EMPTY_INTERVAL.end = 0


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


cdef int init_fragment_index(fragment_index_t* self, int bins_per_dalton=1000, double max_fragment_size=3000) nogil:
    cdef:
        size_t total_bins, i
        int result

    self.bins_per_dalton = bins_per_dalton
    self.max_fragment_size = max_fragment_size
    total_bins = total_bins_for_mass(bins_per_dalton, max_fragment_size)
    self.size = total_bins
    self.sort_type = SortingEnum.unsorted
    self.parent_index = <fragment_list_t*>malloc(sizeof(fragment_list_t))
    init_fragment_list(self.parent_index, 64)
    self.bins = <fragment_list_t*>calloc(total_bins, sizeof(fragment_list_t))
    for i in range(self.size):
        result = init_fragment_list(&self.bins[i], 2)
        if result != 0:
            printf("Error when initializing fragment bin %d\n", i)
            for j in range(i):
                free_fragment_list(&self.bins[j])
            free(self.bins)
            return 1
    return 0


cdef int free_fragment_index(fragment_index_t* self) nogil:
    free_fragment_list(self.parent_index)
    free(self.parent_index)
    for i in range(self.size):
        free_fragment_list(&self.bins[i])
    free(self.bins)
    return 0


cpdef size_t total_bins_for_mass(int bins_per_dalton, double max_fragment_size) nogil:
    return bins_per_dalton * <size_t>_round(max_fragment_size)


cdef size_t bin_for_mass(fragment_index_t* self, double mass) nogil:
    cdef:
        size_t i
    i = <size_t>_round(mass * self.bins_per_dalton)
    if i > self.size:
        i = self.size - 1
    elif i < 0:
        i = 0
    return i


cdef void fragment_index_sort(fragment_index_t* self, SortingEnum sort_type) nogil:
    for i in range(self.size):
        fragment_list_sort(&self.bins[i], sort_type)
    self.sort_type = sort_type


cdef int fragment_index_add_parent(fragment_index_t* self, double mass, uint64_t parent_id) nogil:
    cdef:
        fragment_t f
    if self.parent_index.used > 0 and (self.parent_index.v[self.parent_index.used - 1].mass - mass) > 1e-3:
        return 2
    f.mass = mass
    f.parent_id = parent_id
    f.series = SeriesEnum.parent
    f.ordinal = 0
    return fragment_list_append(self.parent_index, f)


cdef int fragment_index_parents_for(fragment_index_t* self, double mass, double error_tolerance, interval_t* out) nogil:
    cdef:
        int code
    code = fragment_list_binary_search(self.parent_index, mass, error_tolerance, out)
    return code


cdef int fragment_index_parents_for_range(fragment_index_t* self, double low, double high, double error_tolerance, interval_t* out) nogil:
    cdef:
        interval_t tmp
    fragment_index_parents_for(self, low, error_tolerance, &tmp)
    out.start = tmp.start
    fragment_index_parents_for(self, high, error_tolerance, &tmp)
    out.end = tmp.end
    return 0


cdef bint fragment_index_search_has_next(fragment_index_search_t* self) nogil:
    if self.position < self.position_range.end:
        return True
    if self.current_bin < self.high_bin:
        return True
    return False


@cython.cdivision(True)
cdef int fragment_index_search_update_bin(fragment_index_search_t* self) nogil:
    cdef:
        size_t i

    while self.current_bin <= self.high_bin:
        self.current_bin += 1
        # If the index was sorted by mass alone, just use binary search to find the query start-stop interval
        if self.index.sort_type == SortingEnum.by_mass:
            fragment_list_binary_search(&self.index.bins[self.current_bin], self.query, self.error_tolerance, &self.position_range)
            self.position = self.position_range.start
            for i in range(self.position, min(self.position_range.end + 1, self.index.bins[self.current_bin].used)):
                if not interval_contains(&self.parent_id_interval, self.index.bins[self.current_bin].v[i].parent_id):
                    continue
                else:
                    self.position_range.start = i
                    break
            else:
                self.position_range.start = self.index.bins[self.current_bin].used
        else:
            # Otherwise we have to traverse the whole bin, but we can at least advance to the first good position
            self.position_range.start = 0
            self.position_range.end = self.index.bins[self.current_bin].used
            i = 0
            for i in range(self.index.bins[self.current_bin].used):
                if not interval_contains(&self.parent_id_interval, self.index.bins[self.current_bin].v[i].parent_id):
                    continue
                if fabs(self.index.bins[self.current_bin].v[i].mass - self.query) / self.query < self.error_tolerance:
                    self.position_range.start = i
                    break
            else:
                self.position_range.start = self.index.bins[self.current_bin].used

        self.position = self.position_range.start
        if self.position < self.position_range.end:
            break
    return 0


@cython.cdivision(True)
cdef int fragment_index_search_advance(fragment_index_search_t* self) nogil:
    cdef:
        size_t i
        fragment_t fragment

    # If the index was sorted by mass, we're guaranteed to be in the next valid position
    if self.index.sort_type == SortingEnum.by_mass:
        self.position += 1
        fragment_index_search_peek(self, &fragment)
        while not interval_contains(&self.parent_id_interval, fragment.parent_id) and self.position < self.index.bins[self.current_bin].used:
            self.position += 1
    else:
        # Otherwise we have to walk forward incrementally until we find the next
        # valid position or the end of the bin/interval
        self.position += 1
        i = self.position
        for i in range(self.position, self.index.bins[self.current_bin].used):
            if not interval_contains(&self.parent_id_interval, self.index.bins[self.current_bin].v[i].parent_id):
                continue
            if fabs(self.index.bins[self.current_bin].v[i].mass - self.query) / self.query < self.error_tolerance:
                self.position = i
                break
        else:
            self.position = self.index.bins[self.current_bin].used
    return 0


@cython.cdivision(True)
cdef int fragment_index_search_next(fragment_index_search_t* self, fragment_t* fragment) nogil:
    cdef:
        size_t i

    if self.position < self.position_range.end:
        fragment[0] = self.index.bins[self.current_bin].v[self.position]
        fragment_index_search_advance(self)

    # We reached the end of the interval, possibly the end of the bin
    if self.position > self.position_range.end or self.position >= self.index.bins[self.current_bin].used:
        fragment_index_search_update_bin(self)
    return 0


cdef int fragment_index_search_peek(fragment_index_search_t* self, fragment_t* fragment) nogil:
    fragment[0] = self.index.bins[self.current_bin].v[self.position]
    return 0


@cython.cdivision(True)
cdef int fragment_index_search_init_interval(fragment_index_search_t* self) nogil:
    cdef:
        fragment_list_t* fragment_bin
        size_t i
        fragment_t fragment

    fragment_bin = &self.index.bins[self.current_bin]

    self.position_range.start = 0
    self.position_range.end = fragment_bin.used
    self.position = 0
    # printf("Investigating current bin %d with size %d\n", self.current_bin, fragment_bin.used)
    if self.index.sort_type == SortingEnum.by_mass:
        result = fragment_list_binary_search(fragment_bin, self.query, self.error_tolerance, &self.position_range)
        self.position = self.position_range.start
        fragment_index_search_peek(self, &fragment)
        while not interval_contains(&self.parent_id_interval, fragment.parent_id) and self.position < self.index.bins[self.current_bin].used:
            self.position += 1
    else:
        i = 0
        for i in range(fragment_bin.used):
            if not interval_contains(&self.parent_id_interval, self.index.bins[self.current_bin].v[i].parent_id):
                # printf("Skipping %d, not in [%d, %d]\n",
                #        self.index.bins[self.current_bin].v[i].parent_id, self.parent_id_interval.start, self.parent_id_interval.end)
                continue
            if fabs(fragment_bin.v[i].mass - self.query) / self.query < self.error_tolerance:
                # printf("Found starting mass match at %d\n", i)
                self.position_range.start = i
                break
        else:
            self.position_range.start = self.position_range.end

    self.position = self.position_range.start

    # If the result set is empty, make the index think it is empty
    if self.position == self.position_range.end:
        if self.current_bin != self.high_bin:
            self.current_bin += 1
            return fragment_index_search_init_interval(self)
        return 1
    return 0


cdef int fragment_index_search(fragment_index_t* self, double mass, double error_tolerance,
                               fragment_index_search_t* iterator, interval_t parent_id_interval=OPEN_INTERVAL) nogil:
    cdef:
        int result
        size_t low_bin, high_bin, i
        double low, high
        fragment_list_t* fragment_bin
        interval_t bin_range
    iterator.query = mass
    iterator.error_tolerance = error_tolerance
    low = mass - (mass * error_tolerance)
    if low <= 0:
        low_bin = 0
    elif low >= self.max_fragment_size:
        low_bin = self.size - 1
    else:
        low_bin = bin_for_mass(self, low)
    if low_bin != 0:
        low_bin -= 1
        if fragment_list_highest_mass(&self.bins[low_bin]) < low:
            low_bin += 1
    elif low_bin >= self.size:
        low_bin = self.size - 1
    high = mass + (mass * error_tolerance)
    if high >= self.max_fragment_size:
        high_bin = self.size - 1
    else:
        high_bin = bin_for_mass(self, high)
    if high_bin < self.size - 1:
        high_bin += 1
        if fragment_list_lowest_mass(&self.bins[high_bin]) > high:
            high_bin -= 1

    iterator.index = self
    iterator.low_bin = low_bin
    iterator.current_bin = low_bin
    iterator.high_bin = high_bin
    iterator.position_range.start = 0
    iterator.position_range.end = 0
    iterator.parent_id_interval = parent_id_interval

    fragment_index_search_init_interval(iterator)
    return 0


cdef int fragment_index_search_set_parent_interval(fragment_index_search_t* self, interval_t parent_id_interval) nogil:
    cdef:
        fragment_t fragment
        int code
    self.parent_id_interval = parent_id_interval
    if fragment_index_search_has_next(self):
        code = fragment_index_search_peek(self, &fragment)
        if code != 0:
            return 0
        if not interval_contains(&self.parent_id_interval, fragment.parent_id):
            fragment_index_search_next(self, &fragment)
    return 0

cdef int fragment_index_traverse(fragment_index_t* self, fragment_index_traverse_t* iterator, interval_t parent_id_interval=OPEN_INTERVAL) nogil:
    iterator.index = self
    iterator.current_bin = 0
    iterator.position = 0
    iterator.parent_id_interval = parent_id_interval
    while self.bins[iterator.current_bin].used == 0 and iterator.current_bin < self.size:
        iterator.current_bin += 1
    return 0


cdef bint fragment_index_traverse_has_next(fragment_index_traverse_t* self) nogil:
    if self.current_bin < self.index.size - 1:
        return True
    if self.current_bin < self.index.size and self.position < self.index.bins[self.current_bin].used:
        return True
    return False


cdef int fragment_index_traverse_advance(fragment_index_traverse_t* self) nogil:
    self.position += 1


cdef int fragment_index_traverse_update_bin(fragment_index_traverse_t* self) nogil:
    while self.current_bin < self.index.size:
        self.current_bin += 1
        self.position = 0
        if self.position < self.index.bins[self.current_bin].used:
            break


cdef int fragment_index_traverse_next(fragment_index_traverse_t* self, fragment_t* fragment) nogil:
    if self.position < self.index.bins[self.current_bin].used:
        fragment[0] = self.index.bins[self.current_bin].v[self.position]
        fragment_index_traverse_advance(self)
    if self.position == self.index.bins[self.current_bin].used:
        fragment_index_traverse_update_bin(self)
    return 0


cdef int fragment_index_traverse_seek(fragment_index_traverse_t* self, double query, double error_tolerance=1e-5) nogil:
    cdef:
        size_t i
        interval_t q_range
        double lower_bound
    lower_bound = query - (query * error_tolerance)
    i = bin_for_mass(self.index, query)
    if i < 0:
        i = 0
    if i > 0:
        if fragment_list_highest_mass(&self.index.bins[i - 1]) > lower_bound:
            i -= 1
    self.current_bin = i
    if self.index.sort_type == SortingEnum.by_mass:
        fragment_list_binary_search(&self.index.bins[i], lower_bound, error_tolerance, &q_range)
        self.position = q_range.start
        if q_range.start == q_range.end:
            fragment_index_traverse_update_bin(self)
            self.position = 0
        return 0
    else:
        self.position = 0
        while self.position < self.index.bins[self.current_bin].used:
            if (self.index.bins[self.current_bin].v[self.position].mass < lower_bound):
                self.position += 1
            else:
                break
        if self.position == self.index.bins[self.current_bin].used:
            self.position = 0
            # If there's a next bin to go to, go to it
            if self.current_bin < self.index.size - 1:
                fragment_index_traverse_update_bin(self)
            else:
                # otherwise move the iterator position to the last position in the last bin
                self.position = self.index.bins[self.index.size - 1].used - 1
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
            for j in range(i.start, max(i.stop, len(self)), i.step):
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

    cpdef append(self, float32_t mass, SeriesEnum series, uint32_t parent_id, uint16_t ordinal=0):
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


cdef class FragmentIndex(object):

    @staticmethod
    cdef FragmentIndex _create(fragment_index_t* pointer):
        cdef FragmentIndex self = FragmentIndex.__new__(FragmentIndex)
        self.index = pointer
        self.owned = False
        self._wrap_bins()
        return self

    @property
    def bins_per_dalton(self):
        return self.index.bins_per_dalton

    @property
    def max_fragment_size(self):
        return self.index.max_fragment_size

    @property
    def sort_type(self):
        return self.index.sort_type

    def __init__(self, bins_per_dalton=10, max_fragment_size=3000):
        self._init_index(bins_per_dalton, max_fragment_size)

    def _init_index(self, bins_per_dalton, max_fragment_size):
        self.index = <fragment_index_t*>malloc(sizeof(fragment_index_t))
        result = init_fragment_index(self.index, bins_per_dalton, max_fragment_size)
        if result != 0:
            self.index = NULL
            self.owned = False
            raise MemoryError()
        self.owned = True
        self._wrap_bins()

    cpdef _wrap_bins(self):
        self.bins = list()
        for i in range(self.index.size):
            self.bins.append(FragmentList._create(&self.index.bins[i]))
        self.parent_index = FragmentList._create(self.index.parent_index)

    def __dealloc__(self):
        if self.owned and self.index != NULL:
            free_fragment_index(self.index)
            free(self.index)
            del self.bins[:]

    cpdef clear(self, reinit=True):
        bins_per_dalton = self.index.bins_per_dalton
        max_fragment_size = self.index.max_fragment_size
        self.bins = []
        free_fragment_index(self.index)
        free(self.index)
        self.index = NULL
        if reinit:
            self._init_index(bins_per_dalton, max_fragment_size)
            self._wrap_bins()
            self.owned = True
        else:
            self.owned = False
            self.index = NULL

    def __len__(self):
        return self.index.size

    def __getitem__(self, i):
        if i >= self.index.size:
            raise IndexError(i)
        return self.bins[i]

    cpdef size_t bin_for(self, double mass):
        return bin_for_mass(self.index, mass)

    cpdef add(self, double mass, SeriesEnum series, uint32_t parent_id, uint16_t ordinal=0):
        if mass > self.index.max_fragment_size:
            return
        if mass < 0:
            raise ValueError("Mass cannot be negative!")
        value = <uint32_t>_round(mass * self.bins_per_dalton)

        # At this point it must be within the maximum size
        if value >= self.index.size:
            value = self.index.size - 1
        fragment_list_append(&self.index.bins[value], fragment_t(mass, parent_id, series, ordinal))

    cpdef add_parent(self, double mass, uint32_t parent_id):
        cdef:
            int result
        result = fragment_index_add_parent(self.index, mass, parent_id)
        if result == 2:
            raise ValueError("Parents must be added in ascending mass order")
        elif result == 1:
            raise MemoryError()

    cpdef sort(self, SortingEnum sort_type=SortingEnum.by_mass):
        fragment_index_sort(self.index, sort_type)

    cpdef size_t count(self):
        total = 0
        for bin in self.bins:
            total += len(bin)
        return total

    cpdef interval_t parents_for(self, double mass, double error_tolerance=1e-5):
        cdef:
            interval_t out
        fragment_index_parents_for(self.index, mass, error_tolerance, &out)
        return out

    cpdef interval_t parents_for_range(self, double low, double high, double error_tolerance=1-5):
        cdef:
            interval_t out
        fragment_index_parents_for_range(self.index, low, high, error_tolerance, &out)
        return out

    cpdef FragmentIndexSearchIterator search(self, double mass, double error_tolerance=1e-5):
        cdef:
            fragment_index_search_t* iterator
            FragmentIndexSearchIterator iter_obj
        iterator = <fragment_index_search_t*>malloc(sizeof(fragment_index_search_t))
        fragment_index_search(self.index, mass, error_tolerance, iterator)
        iter_obj = FragmentIndexSearchIterator._create(iterator)
        iter_obj.owned = True
        return iter_obj

    cpdef FragmentIndexTraverseIterator traverse(self):
        cdef:
            fragment_index_traverse_t* iterator
            FragmentIndexTraverseIterator iter_obj
        iterator = <fragment_index_traverse_t*>malloc(sizeof(fragment_index_traverse_t))
        fragment_index_traverse(self.index, iterator)
        iter_obj = FragmentIndexTraverseIterator._create(iterator)
        iter_obj.owned = True
        return iter_obj


@cython.final
cdef class FragmentIndexSearchIterator(object):

    @staticmethod
    cdef FragmentIndexSearchIterator _create(fragment_index_search_t* iterator):
        cdef FragmentIndexSearchIterator self = FragmentIndexSearchIterator.__new__(FragmentIndexSearchIterator)
        self.iterator = iterator
        self.owned = False
        return self

    @property
    def current_bin(self):
        return self.iterator.current_bin

    @property
    def position(self):
        return self.iterator.position

    @property
    def position_range(self):
        return self.iterator.position_range

    @property
    def low_bin(self):
        return self.iterator.low_bin

    @property
    def high_bin(self):
        return self.iterator.high_bin

    def __init__(self, *args, **kwargs):
        raise NotImplementedError()

    def _allocate(self):
        self.iterator = <fragment_index_search_t*>malloc(sizeof(fragment_index_search_t))
        self.owned = True

    cpdef bint set_parent_id_range(self, start, end):
        cdef:
            interval_t interval
        interval.start = start
        interval.end = end
        fragment_index_search_set_parent_interval(self.iterator, interval)
        return 0

    cpdef FragmentList all(self):
        cdef:
            int code
            fragment_t f
            fragment_list_t* acc
            FragmentList result
        acc = <fragment_list_t*>malloc(sizeof(fragment_list_t))
        code = init_fragment_list(acc, 32)
        if code != 0:
            raise MemoryError("Cannot initialize fragment list")
        with nogil:
            while fragment_index_search_has_next(self.iterator):
                code = fragment_index_search_next(self.iterator, &f)
                if code != 0:
                    break
                code = fragment_list_append(acc, f)
                if code != 0:
                    with gil:
                        raise MemoryError("Cannot append to fragment list")
        result = FragmentList._create(acc)
        result.owned = True
        return result

    def __dealloc__(self):
        if self.owned:
            free(self.iterator)

    def __next__(self):
        cdef:
            fragment_t f
        if self.iterator != NULL and fragment_index_search_has_next(self.iterator):
            code = fragment_index_search_next(self.iterator, &f)
            if code != 0:
                raise StopIteration()
            return f
        else:
            raise StopIteration()

    def __iter__(self):
        return self


@cython.final
cdef class FragmentIndexTraverseIterator(object):

    @staticmethod
    cdef FragmentIndexTraverseIterator _create(fragment_index_traverse_t* iterator):
        cdef FragmentIndexTraverseIterator self = FragmentIndexTraverseIterator.__new__(FragmentIndexTraverseIterator)
        self.iterator = iterator
        self.owned = False
        return self

    @property
    def current_bin(self):
        return self.iterator.current_bin

    @property
    def position(self):
        return self.iterator.position

    def __init__(self, *args, **kwargs):
        raise NotImplementedError()

    def _allocate(self):
        self.iterator = <fragment_index_traverse_t*>malloc(sizeof(fragment_index_traverse_t))
        self.owned = True

    def __dealloc__(self):
        if self.owned:
            free(self.iterator)

    def __next__(self):
        cdef:
            fragment_t f
        if self.iterator != NULL and fragment_index_traverse_has_next(self.iterator):
            code = fragment_index_traverse_next(self.iterator, &f)
            if code != 0:
                raise StopIteration()
            return f
        else:
            raise StopIteration()

    def __iter__(self):
        return self

    cpdef int seek(self, double query, double error_tolerance=1e-5):
        return fragment_index_traverse_seek(self.iterator, query, error_tolerance)


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
    def parent_id(self):
        return self.fragment.parent_id

    def __repr__(self):
        return "{self.__class__.__name__}({self.mass}, {self.series}, {self.parent_id})".format(self=self)

    def __getitem__(self, key):
        if key == 'mass':
            return self.mass
        if key == 'series':
            return self.series
        if key == 'parent_id':
            return self.parent_id
        raise KeyError(key)