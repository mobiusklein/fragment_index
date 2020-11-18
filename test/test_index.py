import os
from pyteomics import fasta, parser, mass

data_path = os.path.join(os.path.dirname(__file__), "test_data")

from fragment_index import FragmentIndex, SeriesEnum, SortingEnum, PeakList, search_index


def get_test_data(filename):
    return os.path.join(data_path, filename)


def peptide_to_peaklist(peptide):
    pl = PeakList()
    for i in range(1, len(peptide)):
        pl.append(mass.fast_mass(
            peptide[i:], ion_type='y'), i * 100, 1)
        pl.append(mass.fast_mass(peptide[:i], ion_type='b'), i * 100, 1)
    return pl


def digest_proteins(fasta_path):
    peptides = []
    for prot in fasta.FASTA(fasta_path):
        peptides.extend(parser.cleave(prot.sequence, 'trypsin', 0))
    peptides.sort(key=mass.fast_mass)
    return peptides


def test_index():
    path = get_test_data("yeast_glycoproteins.fa")
    peptides = digest_proteins(path)

    index = FragmentIndex()
    for j, peptide in enumerate(peptides):
        peptide_mass = mass.fast_mass(peptide)
        index.add_parent(peptide_mass, j)
        for i in range(1, len(peptide) - 1):
            m = mass.fast_mass(peptide[:i], ion_type='b')
            index.add(m, SeriesEnum.b, j, i)
            m = mass.fast_mass(peptide[i:], ion_type='y')
            index.add(m, SeriesEnum.y, j, len(peptide) - i)

    index.sort(SortingEnum.by_parent)
    assert index.count() == 197755
    assert index.bin_for(113.084) == 1131
    assert index.parents_for_range(500, 1200) == {"end": 5797, "start": 2590}

    search = search = index.search(113.084)
    search.set_parent_id_range(**{"end": 5797, "start": 2590})
    found = list(search)
    assert len(found) == 569
    assert found[-1]['parent_id'] == 5785

    search = index.search(115.084)
    assert len(list(search)) == 0

    peptide = 'NINVLSDICFPLSNNAHDSLPTFNNGSDLFNPLYFAVLNAATPAR'
    pm = mass.fast_mass(peptide)
    pl = peptide_to_peaklist(peptide)
    assert len(pl) == 88
    matches = search_index(index, pl, pm, 10, 500)
    assert len(matches) == 97
    assert matches[0]['parent_id'] == 9791

    return index
