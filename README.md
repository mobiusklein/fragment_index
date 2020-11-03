# Fragment Indexing
A really simplistic product ion index.

```python
import itertools
from pyteomics import mass, parser
import fragment_index

# Ten bins per dalton, maximum fragment mass of 3000.0
index = fragment_index.FragmentIndex(10, 3000.0)

AGP2 = "MALSWVLTVLSLLPLLEAQIPLCANLVPVPITNATLDRITGKWFYIASAFRNEEYNKSVQEIQATFFYFTPNKTEDTIFLREYQTRQNQCFYNSSYLNVQRENGTVSRYEGGREHVAHLLFLRDTKTLMFGSYLDDEKNWGLSFYADKPETTKEQLGEFYEALDCLCIPRSDVMYTDWKKDKCEPLEKQHEKERKQEEGES"

ATG9 = "MFYQPAQNKKQYDDLADIEAQNNVPNTQEVLEAWQESLDSDEDESSPLEESNGFTISEHDDFVKSVPRKNNPTDLLYSGKLLDSDEPPSVHGNSSKVPSKHPSPSFPETTSLRNLQNGSKQKPALPNFNDPHFYNEDVTRSGHPNRSIYTQLPRNEFSNARVLWNRLSARDRVLWRWANVENLDSFLQQVYTYYTGKGLSCIIVHRLFQILTVSFVIGFTTFITSCIDWPAVTPHGSLAGVTKSQCIAQMSPITYLVLWLFLSFLLALWIYYLTDIPRLWQMREFYIHALKIATADMPTVSWQRVLYRLLKLKNVNALTAEDGRVVSLHNMKRLDAYAIANRIMRKDNYFIALINNGIINIELPLLHRRILTHTTEWNINWCIFNFVFDEQGQLRSAFRNPNSRKRLSEELRRRFIVAGFLNCLFAPIVAIYLVIHNFFRYFNEYHKNPGALSTRRYTPLALWTFREYNELQHFFDERINDSYAAASHYVSQFPDFNMIRLFKYISFILGSFTAILVIITVFDPELMVTFEITKDRSVLFYLGLFGSLIAVSRSIIPDETLVFAPEKALRRVITFTHYMPGWWSDNMHSKAVQQEFCSLYSYRIVNLLWEILGILLTPVLLFFTFPSCSQDIVDFFREHTINVEGVGYVCSYAVFQDNPPYESVASLVQSRKISPLIQNKPELSRISFYEQFNTEAPRRDLR"

AUR1 = "MSALSTLKKRLAACNRASQYKLETSLNPMPTFRLLRNTKWSWTHLQYVFLAGNLIFACIVIESPGFWGKFGIACLLAIALTVPLTRQIFFPAIVIITWAILFYSCRFIPERWRPPIWVRVLPTLENILYGSNLSSLLSKTTHSILDILAWVPYGVMHYSAPFIISFILFIFAPPGTLPVWARTFGYMNLFGVLIQMAFPCSPPWYENMYGLEPATYAVRGSPGGLARIDALFGTSIYTDGFSNSPVVFGAFPSLHAGWAMLEALFLSHVFPRYRFCFYGYVLWLCWCTMYLTHHYFVDLVGGMCLAIICFVFAQKLRLPQLQTGKILRWEYEFVIHGHGLSEKTSNSLARTGSPYLLGRDSFTQNPNAVAFMSGLNNMELANTDHEWSVGSSSPEPLPSPAADLIDRPASTTSSIFDASHLP"


peptide_iterator = itertools.chain.from_iterable(parser.cleave(protein, 'trypsin') for protein in (AGP2, ATG9, AUR1))


# Add each peptide's product ions to the index, keyed to the peptide's
# ascending mass order.
for j, peptide in enumerate(sorted(peptide_iterator, key=mass.fast_mass)):
    peptide_mass = mass.fast_mass(peptide)
    print(peptide, peptide_mass)
    index.add_parent(peptide_mass, j)
    for i in range(1, len(peptide) - 1):
        m = mass.fast_mass(peptide[:i], ion_type='b')
        index.add(m, fragment_index.SeriesEnum.b, j, i)
        m = mass.fast_mass(peptide[i:], ion_type='y')
        index.add(m, fragment_index.SeriesEnum.y, j, len(peptide) - i)

index.sort(fragment_index.SortingEnum.by_parent)

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

```

A more serious application using that index:
```python
from math import log
from collections import defaultdict

# Load some peaks to match to the neutral mass fragments
peak_list = load_peaks()

# A holder for some match statistic
scores = defaultdict(float)
interval = parent_range = index.parents_for_range(500, 1000)
parent_id_start = interval['start']
parent_id_end = interval['end']
for mass, intensity in peak_list:
    # Skip to the matches for mass
    for frag in index.search(mass, 1e-5):
        if parent_id <= frag['parent_id'] <= parent_id_end:
            scores[frag['parent_id']] += log(intensity)

best_matches_descending = sorted(scores.items(), key=lambda x: x[1], reverse=True)
```