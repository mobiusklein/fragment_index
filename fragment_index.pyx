# cython: embedsignatures:True
cimport cython
from libc.stdlib cimport malloc, realloc, calloc, free, qsort
from libc.math cimport floor

cdef extern from * nogil:
    int printf (const char *template, ...)
    void qsort (void *base, unsigned short n, unsigned short w, int (*cmp_func)(void*, void*))

cimport numpy as np

ctypedef np.uint64_t uint64_t
ctypedef np.uint32_t uint32_t
ctypedef np.float32_t float32_t
ctypedef np.uint8_t uint8_t


cdef double _round(double x):
    return floor(x + 0.5)


cpdef enum SeriesEnum:
    b = 1
    y = 2
    c = 3
    z = 4


cdef struct fragment_t:
    float32_t mass
    SeriesEnum series
    uint64_t parent_id


cdef int compare_by_mass(const void * a, const void * b) nogil:
    if (<fragment_t*>a).mass < (<fragment_t*>b).mass:
        return -1
    elif (<fragment_t*>a).mass == (<fragment_t*>b).mass:
        return 0
    elif (<fragment_t*>a).mass > (<fragment_t*>b).mass:
        return 1



cdef struct fragment_list_t:
    fragment_t* v
    size_t used
    size_t size


cdef int init_fragment_list(fragment_list_t* self, size_t size):
    self.v = <fragment_t*>malloc(sizeof(fragment_t) * size)
    self.used = 0
    self.size = size
    if self.v == NULL:
        return 1
    return 0


cdef int free_fragment_list(fragment_list_t* self):
    free(self.v)
    return 0


cdef int fragment_list_append(fragment_list_t* self, fragment_t fragment):
    if self.used >= self.size - 1:
        self.v = <fragment_t*>realloc(self.v, sizeof(fragment_t) * self.size * 2)
        if self.v == NULL:
            return 1
        self.size = self.size * 2
    self.v[self.used] = fragment
    self.used += 1
    return 0


cdef void fragment_list_sort(fragment_list_t* self):
    qsort(self.v, self.used, sizeof(fragment_t), compare_by_mass)

cdef struct fragment_index_t:
    fragment_list_t* bins
    size_t size
    int bins_per_dalton
    double max_fragment_size


cdef int init_fragment_index(fragment_index_t* self, int bins_per_dalton=1000, double max_fragment_size=3000):
    cdef:
        size_t total_bins
        int result

    self.bins_per_dalton = bins_per_dalton
    self.max_fragment_size = max_fragment_size
    total_bins = total_bins_for_mass(bins_per_dalton, max_fragment_size)
    self.size = total_bins

    self.bins = <fragment_list_t*>calloc(total_bins, sizeof(fragment_list_t))
    for i in range(self.size):
        result = init_fragment_list(&self.bins[i], 2)
        if result != 0:
            print("Error when initializing fragment bin %d" % i)
            for j in range(i):
                free_fragment_list(&self.bins[j])
            free(self.bins)
            return 1
    return 0


cdef int free_fragment_index(fragment_index_t* self):
    for i in range(self.size):
        free_fragment_list(&self.bins[i])
    free(self.bins)
    return 0


cpdef size_t total_bins_for_mass(int bins_per_dalton, double max_fragment_size):
    return bins_per_dalton * <size_t>_round(max_fragment_size)


cdef size_t bin_for_mass(fragment_index_t* self, double mass):
    return <size_t>_round(mass * self.bins_per_dalton)


cdef void fragment_index_sort(fragment_index_t* self):
    for i in range(self.size):
        fragment_list_sort(&self.bins[i])


cdef class FragmentList(object):
    cdef:
        fragment_list_t* fragments
        public bint owned

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


cdef class FragmentIndex(object):
    cdef:
        fragment_index_t* index
        public list bins
        public bint owned

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