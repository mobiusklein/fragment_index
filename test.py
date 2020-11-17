import itertools
from pyteomics import mass, parser
import fragment_index
import faulthandler
from pyteomics import mass as masslib
faulthandler.enable()
# Ten bins per dalton, maximum fragment mass of 3000.0
index = fragment_index.FragmentIndex(10, 3000.0)

AGP2 = "MALSWVLTVLSLLPLLEAQIPLCANLVPVPITNATLDRITGKWFYIASAFRNEEYNKSVQEIQATFFYFTPNKTEDTIFLREYQTRQNQCFYNSSYLNVQRENGTVSRYEGGREHVAHLLFLRDTKTLMFGSYLDDEKNWGLSFYADKPETTKEQLGEFYEALDCLCIPRSDVMYTDWKKDKCEPLEKQHEKERKQEEGES"

ATG9 = "MFYQPAQNKKQYDDLADIEAQNNVPNTQEVLEAWQESLDSDEDESSPLEESNGFTISEHDDFVKSVPRKNNPTDLLYSGKLLDSDEPPSVHGNSSKVPSKHPSPSFPETTSLRNLQNGSKQKPALPNFNDPHFYNEDVTRSGHPNRSIYTQLPRNEFSNARVLWNRLSARDRVLWRWANVENLDSFLQQVYTYYTGKGLSCIIVHRLFQILTVSFVIGFTTFITSCIDWPAVTPHGSLAGVTKSQCIAQMSPITYLVLWLFLSFLLALWIYYLTDIPRLWQMREFYIHALKIATADMPTVSWQRVLYRLLKLKNVNALTAEDGRVVSLHNMKRLDAYAIANRIMRKDNYFIALINNGIINIELPLLHRRILTHTTEWNINWCIFNFVFDEQGQLRSAFRNPNSRKRLSEELRRRFIVAGFLNCLFAPIVAIYLVIHNFFRYFNEYHKNPGALSTRRYTPLALWTFREYNELQHFFDERINDSYAAASHYVSQFPDFNMIRLFKYISFILGSFTAILVIITVFDPELMVTFEITKDRSVLFYLGLFGSLIAVSRSIIPDETLVFAPEKALRRVITFTHYMPGWWSDNMHSKAVQQEFCSLYSYRIVNLLWEILGILLTPVLLFFTFPSCSQDIVDFFREHTINVEGVGYVCSYAVFQDNPPYESVASLVQSRKISPLIQNKPELSRISFYEQFNTEAPRRDLR"

AUR1 = "MSALSTLKKRLAACNRASQYKLETSLNPMPTFRLLRNTKWSWTHLQYVFLAGNLIFACIVIESPGFWGKFGIACLLAIALTVPLTRQIFFPAIVIITWAILFYSCRFIPERWRPPIWVRVLPTLENILYGSNLSSLLSKTTHSILDILAWVPYGVMHYSAPFIISFILFIFAPPGTLPVWARTFGYMNLFGVLIQMAFPCSPPWYENMYGLEPATYAVRGSPGGLARIDALFGTSIYTDGFSNSPVVFGAFPSLHAGWAMLEALFLSHVFPRYRFCFYGYVLWLCWCTMYLTHHYFVDLVGGMCLAIICFVFAQKLRLPQLQTGKILRWEYEFVIHGHGLSEKTSNSLARTGSPYLLGRDSFTQNPNAVAFMSGLNNMELANTDHEWSVGSSSPEPLPSPAADLIDRPASTTSSIFDASHLP"


peptide_iterator = itertools.chain.from_iterable(parser.cleave(
    protein, 'trypsin') for protein in (AGP2, ATG9, AUR1))


# Add each peptide's product ions to the index, keyed to the peptide's
# ascending mass order.
for j, peptide in enumerate(sorted(peptide_iterator, key=mass.fast_mass)):
    peptide_mass = mass.fast_mass(peptide)
    print(peptide, peptide_mass, j)
    index.add_parent(peptide_mass, j)
    for i in range(1, len(peptide) - 1):
        m = mass.fast_mass(peptide[:i], ion_type='b')
        index.add(m, fragment_index.SeriesEnum.b, j)
        m = mass.fast_mass(peptide[i:], ion_type='y')
        index.add(m, fragment_index.SeriesEnum.y, j)

index.sort(fragment_index.SortingEnum.by_parent)
print(index.count())

# A sequential iterator over fragments in the index that may quickly fast-forward
# towards a specific mass
trav = index.traverse()

f = next(trav)
print(f)

mass = f['mass']
i = index.bin_for(mass)
print("Bin for %f: %d" % (mass, i))
print(list(index.bins[i]))

print("Searching for first fragment")
iterator = index.search(mass, 1e-5)
for f in iterator:
    print(f)

print("Trying a fragment not in the index")
iterator = index.search(37.0320, 1e-5)
for f in iterator:
    print(f)

# fl = index[index.bin_for(57.021)]
# bin_dat = fl.to_bytes()
# print(bin_dat)
# print(fragment_index.FragmentList.from_bytes(bin_dat))

print("Building Peak List")
peptide = "SDVMYTDWK"
pl = fragment_index.PeakList()
for i in range(1, len(peptide) - 1):
    pl.append(masslib.fast_mass(peptide[i:], ion_type='y'), i * 100, 1)

precursor_mass = masslib.fast_mass(peptide)
print("Executing Index Search")
matches = fragment_index.search_index(index, pl, precursor_mass, 200, 700)
print(matches)

import IPython
IPython.embed()
