cimport cython
from libc.stdlib cimport malloc, realloc, calloc, free, qsort
from libc.string cimport memcpy
from libc.math cimport floor, fabs, log10

import numpy as np
cimport numpy as np

cdef extern from * nogil:
    int printf (const char *template, ...)
    void qsort (void *base, unsigned short n, unsigned short w, int (*cmp_func)(void*, void*))


from fragment_index.fragment_index cimport (
    fragment_index_t,
    fragment_list_t,
    fragment_t,
    interval_t,
    fragment_index_search_t,
    fragment_index_parents_for,
    fragment_index_parents_for_range,
    fragment_index_search,
    fragment_index_search_next,
    fragment_index_search_has_next)


np.import_array()


cdef int init_peak_list(peak_list_t* self, size_t size) nogil:
    self.v = <peak_t*>malloc(sizeof(peak_t) * size)
    self.used = 0
    self.size = size
    if self.v == NULL:
        return 1
    return 0


cdef int free_peak_list(peak_list_t* self) nogil:
    free(self.v)
    return 0


cdef int peak_list_append(peak_list_t* self, peak_t peak) nogil:
    if self.used >= self.size - 1:
        self.v = <peak_t*>realloc(self.v, sizeof(peak_t) * self.size * 2)
        if self.v == NULL:
            return 1
        self.size = self.size * 2
    self.v[self.used] = peak
    self.used += 1
    return 0


cdef int init_match_list(match_list_t* self, size_t size) nogil:
    self.v = <match_t*>malloc(sizeof(match_t) * size)
    self.used = 0
    self.size = size
    if self.v == NULL:
        return 1
    return 0


cdef int free_match_list(match_list_t* self) nogil:
    free(self.v)
    return 0


cdef int match_list_append(match_list_t* self, match_t match) nogil:
    if self.used >= self.size - 1:
        self.v = <match_t*>realloc(self.v, sizeof(match_t) * self.size * 2)
        if self.v == NULL:
            return 1
        self.size = self.size * 2
    self.v[self.used] = match
    self.used += 1
    return 0


cdef int compare_by_score_less_than(const void * a, const void * b) nogil:
    if (<match_t*>a).score < (<match_t*>b).score:
        return -1
    elif (<match_t*>a).score == (<match_t*>b).score:
        return 0
    elif (<match_t*>a).score > (<match_t*>b).score:
        return 1

cdef int compare_by_score_greater_than(const void * a, const void * b) nogil:
    return -compare_by_score_less_than(a, b)


cdef int match_list_sort(match_list_t* self) nogil:
    qsort(self.v, self.used, sizeof(match_t), compare_by_score_greater_than)
    return 0


cdef int score_matched_peak(peak_t* peak, fragment_t* fragment, match_t* match) nogil:
    match.score += log10(peak.intensity)
    match.hit_count += 1
    return 0


cdef int search_fragment_index(fragment_index_t* index, peak_list_t* peak_list, double precursor_mass, double parent_error_low,
                               double parent_error_high, double error_tolerance, fragment_search_t* result) nogil:
    cdef:
        interval_t parent_id_interval
        int code
        size_t n_parents, i, parent_offset
        match_list_t* matches
        match_t match
        peak_t* peak
        fragment_t fragment
        fragment_index_search_t iterator


    peak = NULL

    # Initialize match list and parent_id_interval
    parent_id_interval.start = 0
    parent_id_interval.end = -1
    fragment_index_parents_for(index, precursor_mass - parent_error_low, 1e-5, &parent_id_interval)
    fragment_index_parents_for(index, precursor_mass + parent_error_high, 1e-5, &parent_id_interval)
    fragment_index_parents_for_range(
        index,
        precursor_mass - parent_error_low,
        precursor_mass + parent_error_high,
        1e-5,
        &parent_id_interval)
    n_parents = parent_id_interval.end - parent_id_interval.start + 1
    matches = <match_list_t*>malloc(sizeof(match_list_t))
    if matches == NULL:
        return 1
    code = init_match_list(matches, n_parents)
    if code != 0:
        return 1
    for i in range(n_parents):
        match.parent_id = parent_id_interval.start + i
        match.score = 0
        match.hit_count = 0
        match_list_append(matches, match)

    # Search the index for each peak in the peak list
    for i in range(peak_list.used):
        peak = &peak_list.v[i]
        code = fragment_index_search(index, peak_list.v[i].mass, error_tolerance, &iterator, parent_id_interval)
        if code != 0:
            return 2
        while fragment_index_search_has_next(&iterator):
            code = fragment_index_search_next(&iterator, &fragment)
            if code != 0:
                break
            parent_offset = fragment.parent_id
            if parent_offset < parent_id_interval.start:
                printf("Parent ID %d outside of expected interval [%d, %d] for mass %f\n",
                       parent_offset, parent_id_interval.start, parent_id_interval.end, peak_list.v[i].mass)
                return 3
            parent_offset -= parent_id_interval.start
            score_matched_peak(peak, &fragment, &matches.v[parent_offset])
    match_list_sort(matches)
    result.index = index
    result.peak_list = peak_list
    result.match_list = matches
    result.parent_interval = parent_id_interval
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
        if self.owned and self.peaks != NULL:
            free_peak_list(self.peaks)
            free(self.peaks)
            self.peaks = NULL

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

    cpdef append(self, float32_t mass, float32_t intensity, int charge):
        cdef peak_t peak = peak_t(mass, intensity, charge)
        out = peak_list_append(self.peaks, peak)
        if out == 1:
            raise MemoryError()

    @classmethod
    def from_arrays(cls, numeric_collection1 mz_array, numeric_collection2 intensity_array, object charge_array=None):
        cdef:
            PeakList self
            Py_ssize_t i
            float32_t mz
            float32_t mass
            float32_t intensity
            int charge

        self = cls()
        n = len(mz_array)
        for i in range(n):
            mz = mz_array[i]
            intensity = intensity_array[i]
            if charge_array is not None:
                charge = charge_array[i]
            else:
                charge = 1
            if charge != 0:
                mass = neutral_mass(mz, charge)
            else:
                mass = mz
            self.append(mass, intensity, charge)
        return self


cdef double PROTON = 1.00727646677


cdef double neutral_mass(double mz,  int z, double charge_carrier=PROTON) nogil:
    return (mz * fabs(z)) - (z * charge_carrier)

@cython.cdivision
cdef double mass_charge_ratio(double neutral_mass, int z, double charge_carrier=PROTON) nogil:
    return (neutral_mass + (z * charge_carrier)) / fabs(z)


cdef class MatchList(object):

    @staticmethod
    cdef MatchList _create(match_list_t* pointer):
        cdef MatchList self = MatchList.__new__(MatchList)
        self.matches = pointer
        self.owned = False
        return self

    @property
    def allocated(self):
        return self.matches.size

    def __init__(self, *args, **kwargs):
        self._init_list()

    cdef void _init_list(self):
        self.matches = <match_list_t*>malloc(sizeof(match_list_t))
        self.owned = True
        init_match_list(self.matches, 32)

    cpdef clear(self):
        free_match_list(self.matches)
        free(self.matches)
        self._init_list()

    def __dealloc__(self):
        if self.owned and self.matches != NULL:
            free_match_list(self.matches)
            free(self.matches)
            self.matches = NULL

    def __len__(self):
        return self.matches.used

    def __getitem__(self, i):
        if isinstance(i, slice):
            out = []
            for j in range(i.start or 0, min(i.stop or len(self), len(self)), i.step or 1):
                out.append(self[j])
            return out
        if i  >= self.matches.used:
            raise IndexError(i)
        elif i < 0:
            j = len(self) + i
            if j < 0:
                raise IndexError(i)
            i = j
        return self.matches.v[i]

    def __iter__(self):
        for i in range(self.matches.used):
            yield self.matches.v[i]

    def __repr__(self):
        return "{self.__class__.__name__}({size})".format(self=self, size=len(self))

    cpdef append(self, uint32_t parent_id, float32_t score, uint32_t hit_count):
        cdef match_t match = match_t(parent_id, score, hit_count)
        out = match_list_append(self.matches, match)
        if out == 1:
            raise MemoryError()


def search_index(FragmentIndex index, PeakList peaks, double precursor_mass, double parent_error_low, double parent_error_high, double error_tolerance=2e-5):
    cdef:
        fragment_search_t* search_result
        int code
        MatchList matches
    search_result = <fragment_search_t*>malloc(sizeof(fragment_search_t))
    if search_result == NULL:
        raise MemoryError()
    search_result.index = NULL
    search_result.peak_list = NULL
    search_result.match_list = NULL
    with nogil:
        code = search_fragment_index(
            index.index,
            peaks.peaks,
            precursor_mass=precursor_mass,
            parent_error_low=parent_error_low,
            parent_error_high=parent_error_high,
            error_tolerance=error_tolerance,
            result=search_result)

    if code != 0:
        raise ValueError()

    matches = MatchList._create(search_result.match_list)
    matches.owned = True
    free(search_result)
    return matches