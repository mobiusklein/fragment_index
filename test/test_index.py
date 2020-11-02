import os
from pyteomics import fasta, parser, mass

data_path = os.path.join(os.path.dirname(__file__), "test_data")

from fragment_index import FragmentIndex, SeriesEnum, SortingEnum


def get_test_data(filename):
    return os.path.join(data_path, filename)


def test_index():
    path = get_test_data("yeast_glycoproteins.fa")
    peptides = []
    for prot in fasta.FASTA(path):
        peptides.extend(parser.cleave(prot.sequence, 'trypsin', 0))
    peptides.sort(key=mass.fast_mass)

    index = FragmentIndex()
    for j, peptide in enumerate(peptides):
        peptide_mass = mass.fast_mass(peptide)
        index.add_parent(peptide_mass, j)
        for i in range(1, len(peptide) - 1):
            m = mass.fast_mass(peptide[:i], ion_type='b')
            index.add(m, SeriesEnum.b, j)
            m = mass.fast_mass(peptide[i:], ion_type='y')
            index.add(m, SeriesEnum.y, j)

    index.sort(SortingEnum.by_parent)
    assert index.count() == 197755
    assert index.bin_for(113.084) == 1131
    assert index.parents_for_range(500, 1200) == {"end": 5796, "start": 2590}

    search = search = index.search(113.084)
    search.set_parent_id_range(**{"end": 5796, "start": 2590})
    found = list(search)
    assert len(found) == 569
    assert found[-1]['parent_id'] == 5785

    search = index.search(115.084)
    assert len(list(search)) == 0

    return index
