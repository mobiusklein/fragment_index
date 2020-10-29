# cython: embedsignatures:True

cimport cython
from libc.stdlib cimport malloc, realloc, calloc, free, qsort
from libc.math cimport floor, fabs

cdef extern from * nogil:
    int printf (const char *template, ...)
    void qsort (void *base, unsigned short n, unsigned short w, int (*cmp_func)(void*, void*))


cdef double _round(double x) nogil:
    return floor(x + 0.5)


cdef int compare_by_mass(const void * a, const void * b) nogil:
    if (<fragment_t*>a).mass < (<fragment_t*>b).mass:
        return -1
    elif (<fragment_t*>a).mass == (<fragment_t*>b).mass:
        return 0
    elif (<fragment_t*>a).mass > (<fragment_t*>b).mass:
        return 1


cdef int init_fragment_list(fragment_list_t* self, size_t size) nogil:
    self.v = <fragment_t*>malloc(sizeof(fragment_t) * size)
    self.used = 0
    self.size = size
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
    return 0


cdef void fragment_list_sort(fragment_list_t* self) nogil:
    qsort(self.v, self.used, sizeof(fragment_t), compare_by_mass)


cdef double fragment_list_lowest_mass(fragment_list_t* self) nogil:
    if self.used == 0:
        return 0
    return self.v[0].mass


cdef double fragment_list_highest_mass(fragment_list_t* self) nogil:
    if self.used == 0:
        return 0
    return self.v[self.used - 1].mass


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
        if lo == hi - 1 or err < error_tolerance:
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
            out.end = i
            return 0
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return 1


cdef int init_fragment_index(fragment_index_t* self, int bins_per_dalton=1000, double max_fragment_size=3000) nogil:
    cdef:
        size_t total_bins, i
        int result

    self.bins_per_dalton = bins_per_dalton
    self.max_fragment_size = max_fragment_size
    total_bins = total_bins_for_mass(bins_per_dalton, max_fragment_size)
    self.size = total_bins

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


cdef void fragment_index_sort(fragment_index_t* self) nogil:
    for i in range(self.size):
        fragment_list_sort(&self.bins[i])


cdef bint fragment_index_search_has_next(fragment_index_search_t* self) nogil:
    if self.position < self.position_range.end:
        return True
    if self.current_bin < self.high_bin:
        return True
    return False


cdef int fragment_index_search_next(fragment_index_search_t* self, fragment_t* fragment) nogil:
    if self.position < self.position_range.end:
        fragment[0] = self.index.bins[self.current_bin].v[self.position]
        self.position += 1
    if self.position == self.position_range.end:
        while self.current_bin <= self.high_bin:
            self.current_bin += 1
            fragment_list_binary_search(&self.index.bins[self.current_bin], self.query, self.error_tolerance, &self.position_range)
            self.position = self.position_range.start
            if self.position < self.position_range.end:
                break
    return 0


cdef int fragment_index_search(fragment_index_t* self, double mass, double error_tolerance, fragment_index_search_t* iterator) nogil:
    cdef:
        int result
        size_t low_bin, high_bin
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

    fragment_bin = &self.bins[low_bin]
    result = fragment_list_binary_search(fragment_bin, mass, error_tolerance, &iterator.position_range)
    iterator.position = iterator.position_range.start

    # If the result set is empty, make the iterator think it is empty
    if iterator.position == iterator.position_range.end:
        iterator.current_bin = iterator.high_bin
    return 0


cdef int fragment_index_traverse(fragment_index_t* self, fragment_index_traverse_t* iterator) nogil:
    iterator.index = self
    iterator.current_bin = 0
    iterator.position = 0
    while self.bins[iterator.current_bin].used == 0 and iterator.current_bin < self.size:
        iterator.current_bin += 1
    return 0


cdef bint fragment_index_traverse_has_next(fragment_index_traverse_t* self) nogil:
    if self.current_bin < self.index.size - 1:
        return True
    if self.current_bin < self.index.size and self.position < self.index.bins[self.current_bin].used:
        return True
    return False


cdef int fragment_index_traverse_next(fragment_index_traverse_t* self, fragment_t* fragment) nogil:
    if self.position < self.index.bins[self.current_bin].used:
        fragment[0] = self.index.bins[self.current_bin].v[self.position]
        self.position += 1
    if self.position == self.index.bins[self.current_bin].used:
        while self.current_bin < self.index.size:
            self.current_bin += 1
            self.position = 0
            if self.position < self.index.bins[self.current_bin].used:
                break
    return 0


cdef int fragment_index_traverse_seek(fragment_index_traverse_t* self, double query, double error_tolerance=1e-5) nogil:
    cdef:
        size_t i
        interval_t q_range
    query = query - (query * error_tolerance)
    i = bin_for_mass(self.index, query)
    if i < 0:
        i = 0
    if i > 0:
        if fragment_list_highest_mass(&self.index.bins[i - 1]) > query:
            i -= 1
    self.current_bin = i
    fragment_list_binary_search(&self.index.bins[i], query, error_tolerance, &q_range)
    self.position = q_range.start
    if q_range.start == q_range.end:
        return 1
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

    cpdef append(self, float32_t mass, SeriesEnum series, uint64_t parent_id):
        cdef fragment_t fragment = fragment_t(mass, series, parent_id)
        out = fragment_list_append(self.fragments, fragment)
        if out == 1:
            raise MemoryError()

    cpdef sort(self):
        fragment_list_sort(self.fragments)

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

    cpdef add(self, double mass, SeriesEnum series, uint64_t parent_id):
        if mass > self.index.max_fragment_size:
            return
        if mass < 0:
            raise ValueError("Mass cannot be negative!")
        value = <uint32_t>_round(mass * self.bins_per_dalton)

        # At this point it must be within the maximum size
        if value >= self.index.size:
            value = self.index.size - 1
        fragment_list_append(&self.index.bins[value], fragment_t(mass, series, parent_id))

    cpdef sort(self):
        fragment_index_sort(self.index)

    cpdef size_t count(self):
        total = 0
        for bin in self.bins:
            total += len(bin)
        return total

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
