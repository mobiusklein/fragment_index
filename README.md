# Fragment Indexing
A really simplistic product ion index.

```python
from pyteomics import mass, parser
import fragment_index

# Ten bins per dalton, maximum fragment mass of 3000.0
index = fragment_index.FragmentIndex(10, 3000.0)

AGP2 = "MALSWVLTVLSLLPLLEAQIPLCANLVPVPITNATLDRITGKWFYIASAFRNEEYNKSVQEIQATFFYFTPNKTEDTIFLREYQTRQNQCFYNSSYLNVQRENGTVSRYEGGREHVAHLLFLRDTKTLMFGSYLDDEKNWGLSFYADKPETTKEQLGEFYEALDCLCIPRSDVMYTDWKKDKCEPLEKQHEKERKQEEGES"

# Add each peptide's product ions to the index, keyed to the peptide's
# yield order.
precursor_index = []
for j, peptide in enumerate(parser.cleave(AGP2, 'trypsin')):
    print(peptide)
    precursor_index.append(j)
    for i in range(1, len(peptide) - 1):
        m = mass.fast_mass(peptide[:i], ion_type='b')
        index.add(m, fragment_index.SeriesEnum.b, j)
        m = mass.fast_mass(peptide[i:], ion_type='y')
        index.add(m, fragment_index.SeriesEnum.y, j)
index.sort()
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
iterator = index.search(57.0320, 1e-5)
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

for mass, intensity in peak_list:
    # Skip to the matches for mass
    for frag in index.search(mass, 1e-5):
        scores[frag['parent_id']] += log(intensity)

best_matches_descending = sorted(scores.items(), key=lambda x: x[1], reverse=True)
```