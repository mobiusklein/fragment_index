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


include "peak_list.pyx"
include "parent_list.pyx"

# Fragment Index Methods

cdef int init_peak_index(peak_index_t* self, int bins_per_dalton=10, double max_peak_size=3000) nogil:
    cdef:
        size_t total_bins, i
        int result

    self.bins_per_dalton = bins_per_dalton
    self.max_peak_size = max_peak_size
    total_bins = total_bins_for_mass(bins_per_dalton, max_peak_size)
    self.size = total_bins
    self.sort_type = SortingEnum.unsorted
    self.parent_index = <parent_list_t*>malloc(sizeof(parent_list_t))
    init_parent_list(self.parent_index, 64)
    self.bins = <peak_list_t*>calloc(total_bins, sizeof(peak_list_t))
    for i in range(self.size):
        result = init_peak_list(&self.bins[i], 2)
        if result != 0:
            printf("Error when initializing peak bin %d\n", i)
            for j in range(i):
                free_peak_list(&self.bins[j])
            free(self.bins)
            return 1
    return 0


cdef int free_peak_index(peak_index_t* self) nogil:
    free_parent_list(self.parent_index)
    free(self.parent_index)
    for i in range(self.size):
        free_peak_list(&self.bins[i])
    free(self.bins)
    return 0


cpdef size_t total_bins_for_mass(int bins_per_dalton, double max_peak_size) nogil:
    return bins_per_dalton * <size_t>_round(max_peak_size)


cdef size_t bin_for_mass(peak_index_t* self, double mass) nogil:
    cdef:
        size_t i
    i = <size_t>_round(mass * self.bins_per_dalton)
    if i > self.size:
        i = self.size - 1
    elif i < 0:
        i = 0
    return i


cdef void peak_index_sort(peak_index_t* self, SortingEnum sort_type) nogil:
    for i in range(self.size):
        peak_list_sort(&self.bins[i], sort_type)
    self.sort_type = sort_type


cdef int peak_index_add_parent(peak_index_t* self, double mass, int_id_t id, int_id_t scan_id=0, uint16_t start_position=0, uint16_t size=0) nogil:
    cdef:
        parent_t f
    if self.parent_index.used > 0 and (self.parent_index.v[self.parent_index.used - 1].mass - mass) > 1e-3:
        return 2
    f.mass = mass
    f.id = id
    f.parent_id = scan_id
    f.start_position = start_position
    f.size = size
    return parent_list_append(self.parent_index, f)


cdef int peak_index_parents_for(peak_index_t* self, double mass, double error_tolerance, interval_t* out) nogil:
    cdef:
        int code
    code = parent_list_binary_search(self.parent_index, mass, error_tolerance, out)
    return code


cdef int peak_index_parents_for_range(peak_index_t* self, double low, double high, double error_tolerance, interval_t* out) nogil:
    cdef:
        interval_t tmp
    peak_index_parents_for(self, low, error_tolerance, &tmp)
    out.start = tmp.start
    peak_index_parents_for(self, high, error_tolerance, &tmp)
    out.end = tmp.end
    return 0


# Fragment Index Search Methods

cdef bint peak_index_search_has_next(peak_index_search_t* self) nogil:
    if self.position < self.position_range.end:
        return True
    if self.current_bin < self.high_bin:
        return True
    return False


@cython.cdivision(True)
cdef int peak_index_search_update_bin(peak_index_search_t* self) nogil:
    cdef:
        size_t i

    while self.current_bin <= self.high_bin:
        self.current_bin += 1
        # If the index was sorted by mass alone, just use binary search to find the query start-stop interval
        if self.index.sort_type == SortingEnum.by_mass:
            peak_list_binary_search(&self.index.bins[self.current_bin], self.query, self.error_tolerance, &self.position_range)
            self.position = self.position_range.start
            for i in range(self.position, min(self.position_range.end + 1, self.index.bins[self.current_bin].used)):
                if not interval_contains(&self.scan_id_interval, self.index.bins[self.current_bin].v[i].scan_id):
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
                if not interval_contains(&self.scan_id_interval, self.index.bins[self.current_bin].v[i].scan_id):
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
cdef int peak_index_search_advance(peak_index_search_t* self) nogil:
    cdef:
        size_t i
        peak_t peak

    # If the index was sorted by mass, we're guaranteed to be in the next valid position
    if self.index.sort_type == SortingEnum.by_mass:
        self.position += 1
        peak_index_search_peek(self, &peak)
        while not interval_contains(&self.scan_id_interval, peak.scan_id) and self.position < self.index.bins[self.current_bin].used:
            self.position += 1
    else:
        # Otherwise we have to walk forward incrementally until we find the next
        # valid position or the end of the bin/interval
        self.position += 1
        i = self.position
        for i in range(self.position, self.index.bins[self.current_bin].used):
            if not interval_contains(&self.scan_id_interval, self.index.bins[self.current_bin].v[i].scan_id):
                continue
            if fabs(self.index.bins[self.current_bin].v[i].mass - self.query) / self.query < self.error_tolerance:
                self.position = i
                break
        else:
            self.position = self.index.bins[self.current_bin].used
    return 0


@cython.cdivision(True)
cdef int peak_index_search_next(peak_index_search_t* self, peak_t* peak) nogil:
    cdef:
        size_t i

    if self.position < self.position_range.end:
        peak[0] = self.index.bins[self.current_bin].v[self.position]
        peak_index_search_advance(self)

    # We reached the end of the interval, possibly the end of the bin
    if self.position > self.position_range.end or self.position >= self.index.bins[self.current_bin].used:
        peak_index_search_update_bin(self)
    return 0


cdef int peak_index_search_peek(peak_index_search_t* self, peak_t* peak) nogil:
    peak[0] = self.index.bins[self.current_bin].v[self.position]
    return 0


@cython.cdivision(True)
cdef int peak_index_search_init_interval(peak_index_search_t* self) nogil:
    cdef:
        peak_list_t* peak_bin
        size_t i
        peak_t peak

    peak_bin = &self.index.bins[self.current_bin]

    self.position_range.start = 0
    self.position_range.end = peak_bin.used
    self.position = 0
    # printf("Investigating current bin %d with size %d\n", self.current_bin, peak_bin.used)
    if self.index.sort_type == SortingEnum.by_mass:
        result = peak_list_binary_search(peak_bin, self.query, self.error_tolerance, &self.position_range)
        self.position = self.position_range.start
        peak_index_search_peek(self, &peak)
        while not interval_contains(&self.scan_id_interval, peak.scan_id) and self.position < self.index.bins[self.current_bin].used:
            self.position += 1
    else:
        i = 0
        for i in range(peak_bin.used):
            if not interval_contains(&self.scan_id_interval, self.index.bins[self.current_bin].v[i].scan_id):
                # printf("Skipping %d, not in [%d, %d]\n",
                #        self.index.bins[self.current_bin].v[i].scan_id, self.scan_id_interval.start, self.scan_id_interval.end)
                continue
            if fabs(peak_bin.v[i].mass - self.query) / self.query < self.error_tolerance:
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
            return peak_index_search_init_interval(self)
        return 1
    return 0


cdef int peak_index_search(peak_index_t* self, double mass, double error_tolerance,
                           peak_index_search_t* iterator, interval_t scan_id_interval=OPEN_INTERVAL) nogil:
    cdef:
        int result
        size_t low_bin, high_bin, i
        double low, high
        peak_list_t* peak_bin
        interval_t bin_range
    iterator.query = mass
    iterator.error_tolerance = error_tolerance
    low = mass - (mass * error_tolerance)
    if low <= 0:
        low_bin = 0
    elif low >= self.max_peak_size:
        low_bin = self.size - 1
    else:
        low_bin = bin_for_mass(self, low)
    if low_bin != 0:
        low_bin -= 1
        if peak_list_highest_mass(&self.bins[low_bin]) < low:
            low_bin += 1
    elif low_bin >= self.size:
        low_bin = self.size - 1
    high = mass + (mass * error_tolerance)
    if high >= self.max_peak_size:
        high_bin = self.size - 1
    else:
        high_bin = bin_for_mass(self, high)
    if high_bin < self.size - 1:
        high_bin += 1
        if peak_list_lowest_mass(&self.bins[high_bin]) > high:
            high_bin -= 1

    iterator.index = self
    iterator.low_bin = low_bin
    iterator.current_bin = low_bin
    iterator.high_bin = high_bin
    iterator.position_range.start = 0
    iterator.position_range.end = 0
    iterator.scan_id_interval = scan_id_interval

    peak_index_search_init_interval(iterator)
    return 0


cdef int peak_index_search_set_parent_interval(peak_index_search_t* self, interval_t scan_id_interval) nogil:
    cdef:
        peak_t peak
        int code
    self.scan_id_interval = scan_id_interval
    if peak_index_search_has_next(self):
        code = peak_index_search_peek(self, &peak)
        if code != 0:
            return 0
        if not interval_contains(&self.scan_id_interval, peak.scan_id):
            peak_index_search_next(self, &peak)
    return 0

cdef int peak_index_traverse(peak_index_t* self, peak_index_traverse_t* iterator, interval_t scan_id_interval=OPEN_INTERVAL) nogil:
    iterator.index = self
    iterator.current_bin = 0
    iterator.position = 0
    iterator.scan_id_interval = scan_id_interval
    while self.bins[iterator.current_bin].used == 0 and iterator.current_bin < self.size:
        iterator.current_bin += 1
    return 0


cdef bint peak_index_traverse_has_next(peak_index_traverse_t* self) nogil:
    if self.current_bin < self.index.size - 1:
        return True
    if self.current_bin < self.index.size and self.position < self.index.bins[self.current_bin].used:
        return True
    return False


cdef int peak_index_traverse_advance(peak_index_traverse_t* self) nogil:
    self.position += 1


cdef int peak_index_traverse_update_bin(peak_index_traverse_t* self) nogil:
    while self.current_bin < self.index.size:
        self.current_bin += 1
        self.position = 0
        if self.position < self.index.bins[self.current_bin].used:
            break


cdef int peak_index_traverse_next(peak_index_traverse_t* self, peak_t* peak) nogil:
    if self.position < self.index.bins[self.current_bin].used:
        peak[0] = self.index.bins[self.current_bin].v[self.position]
        peak_index_traverse_advance(self)
    if self.position == self.index.bins[self.current_bin].used:
        peak_index_traverse_update_bin(self)
    return 0


cdef int peak_index_traverse_seek(peak_index_traverse_t* self, double query, double error_tolerance=1e-5) nogil:
    cdef:
        size_t i
        interval_t q_range
        double lower_bound
    lower_bound = query - (query * error_tolerance)
    i = bin_for_mass(self.index, query)
    if i < 0:
        i = 0
    if i > 0:
        if peak_list_highest_mass(&self.index.bins[i - 1]) > lower_bound:
            i -= 1
    self.current_bin = i
    if self.index.sort_type == SortingEnum.by_mass:
        peak_list_binary_search(&self.index.bins[i], lower_bound, error_tolerance, &q_range)
        self.position = q_range.start
        if q_range.start == q_range.end:
            peak_index_traverse_update_bin(self)
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
                peak_index_traverse_update_bin(self)
            else:
                # otherwise move the iterator position to the last position in the last bin
                self.position = self.index.bins[self.index.size - 1].used - 1
        return 0


cdef class PeakIndex(object):

    @staticmethod
    cdef PeakIndex _create(peak_index_t* pointer):
        cdef PeakIndex self = PeakIndex.__new__(PeakIndex)
        self.index = pointer
        self.owned = False
        self._wrap_bins()
        return self

    @property
    def bins_per_dalton(self):
        return self.index.bins_per_dalton

    @property
    def max_peak_size(self):
        return self.index.max_peak_size

    @property
    def sort_type(self):
        return self.index.sort_type

    def __init__(self, bins_per_dalton=10, max_peak_size=3000):
        self._init_index(bins_per_dalton, max_peak_size)

    def _init_index(self, bins_per_dalton, max_peak_size):
        self.index = <peak_index_t*>malloc(sizeof(peak_index_t))
        result = init_peak_index(self.index, bins_per_dalton, max_peak_size)
        if result != 0:
            self.index = NULL
            self.owned = False
            raise MemoryError()
        self.owned = True
        self._wrap_bins()

    cpdef _wrap_bins(self):
        self.bins = list()
        for i in range(self.index.size):
            self.bins.append(PeakList._create(&self.index.bins[i]))
        self.parent_index = ParentList._create(self.index.parent_index)

    def __dealloc__(self):
        if self.owned and self.index != NULL:
            free_peak_index(self.index)
            free(self.index)
            del self.bins[:]

    cpdef clear(self, reinit=True):
        bins_per_dalton = self.index.bins_per_dalton
        max_peak_size = self.index.max_peak_size
        self.bins = []
        free_peak_index(self.index)
        free(self.index)
        self.index = NULL
        if reinit:
            self._init_index(bins_per_dalton, max_peak_size)
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

    cpdef add(self, float32_t mass, float32_t intensity, int16_t charge, int_id_t scan_id=0):
        if mass > self.index.max_peak_size:
            return
        if mass < 0:
            raise ValueError("Mass cannot be negative!")
        value = <uint32_t>_round(mass * self.bins_per_dalton)

        # At this point it must be within the maximum size
        if value >= self.index.size:
            value = self.index.size - 1
        peak_list_append(&self.index.bins[value], peak_t(mass, intensity, charge, scan_id))

    cpdef add_parent(self, double mass, int_id_t id, int_id_t scan_id=0, uint16_t start_position=0, uint16_t size=0):
        cdef:
            int result
        result = peak_index_add_parent(self.index, mass, id, scan_id, start_position, size)
        if result == 2:
            raise ValueError("Parents must be added in ascending mass order")
        elif result == 1:
            raise MemoryError()

    cpdef sort(self, SortingEnum sort_type=SortingEnum.by_mass):
        peak_index_sort(self.index, sort_type)

    cpdef size_t count(self):
        total = 0
        for bin in self.bins:
            total += len(bin)
        return total

    cpdef interval_t parents_for(self, double mass, double error_tolerance=1e-5):
        cdef:
            interval_t out
        peak_index_parents_for(self.index, mass, error_tolerance, &out)
        return out

    cpdef interval_t parents_for_range(self, double low, double high, double error_tolerance=1-5):
        cdef:
            interval_t out
        peak_index_parents_for_range(self.index, low, high, error_tolerance, &out)
        return out

    cpdef PeakIndexSearchIterator search(self, double mass, double error_tolerance=1e-5):
        cdef:
            peak_index_search_t* iterator
            PeakIndexSearchIterator iter_obj
        iterator = <peak_index_search_t*>malloc(sizeof(peak_index_search_t))
        peak_index_search(self.index, mass, error_tolerance, iterator)
        iter_obj = PeakIndexSearchIterator._create(iterator)
        iter_obj.owned = True
        return iter_obj

    cpdef PeakIndexTraverseIterator traverse(self):
        cdef:
            peak_index_traverse_t* iterator
            PeakIndexTraverseIterator iter_obj
        iterator = <peak_index_traverse_t*>malloc(sizeof(peak_index_traverse_t))
        peak_index_traverse(self.index, iterator)
        iter_obj = PeakIndexTraverseIterator._create(iterator)
        iter_obj.owned = True
        return iter_obj


@cython.final
cdef class PeakIndexSearchIterator(object):

    @staticmethod
    cdef PeakIndexSearchIterator _create(peak_index_search_t* iterator):
        cdef PeakIndexSearchIterator self = PeakIndexSearchIterator.__new__(PeakIndexSearchIterator)
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
        self.iterator = <peak_index_search_t*>malloc(sizeof(peak_index_search_t))
        self.owned = True

    cpdef bint set_scan_id_range(self, start, end):
        cdef:
            interval_t interval
        interval.start = start
        interval.end = end
        peak_index_search_set_parent_interval(self.iterator, interval)
        return 0

    cpdef PeakList all(self):
        cdef:
            int code
            peak_t f
            peak_list_t* acc
            PeakList result
        acc = <peak_list_t*>malloc(sizeof(peak_list_t))
        code = init_peak_list(acc, 32)
        if code != 0:
            raise MemoryError("Cannot initialize peak list")
        with nogil:
            while peak_index_search_has_next(self.iterator):
                code = peak_index_search_next(self.iterator, &f)
                if code != 0:
                    break
                code = peak_list_append(acc, f)
                if code != 0:
                    with gil:
                        raise MemoryError("Cannot append to peak list")
        result = PeakList._create(acc)
        result.owned = True
        return result

    def __dealloc__(self):
        if self.owned:
            free(self.iterator)

    def __next__(self):
        cdef:
            peak_t f
        if self.iterator != NULL and peak_index_search_has_next(self.iterator):
            code = peak_index_search_next(self.iterator, &f)
            if code != 0:
                raise StopIteration()
            return f
        else:
            raise StopIteration()

    def __iter__(self):
        return self


@cython.final
cdef class PeakIndexTraverseIterator(object):

    @staticmethod
    cdef PeakIndexTraverseIterator _create(peak_index_traverse_t* iterator):
        cdef PeakIndexTraverseIterator self = PeakIndexTraverseIterator.__new__(PeakIndexTraverseIterator)
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
        self.iterator = <peak_index_traverse_t*>malloc(sizeof(peak_index_traverse_t))
        self.owned = True

    def __dealloc__(self):
        if self.owned:
            free(self.iterator)

    def __next__(self):
        cdef:
            peak_t f
        if self.iterator != NULL and peak_index_traverse_has_next(self.iterator):
            code = peak_index_traverse_next(self.iterator, &f)
            if code != 0:
                raise StopIteration()
            return f
        else:
            raise StopIteration()

    def __iter__(self):
        return self

    cpdef int seek(self, double query, double error_tolerance=1e-5):
        return peak_index_traverse_seek(self.iterator, query, error_tolerance)
