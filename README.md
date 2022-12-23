# Residue Interaction Network Plotter for MD Simulations

The following notebook is composed of functions from the MDAnalysis and ProLif documentations that are altered to take inputs for your protien of interest. 
**NOTE** Be careful with your selections! When selecting for a specific segid or residue type, make sure you have an intimate knowledge of your own topology and trjaectory files.

Citations:
Bouysset, C., Fiorucci, S. ProLIF: a library to encode molecular interactions as fingerprints.
J Cheminform 13, 72 (2021). https://doi.org/10.1186/s13321-021-00548-6
N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319–2327. doi:10.1002/jcc.21787
R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler, D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein. MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations. In S. Benthall and S. Rostrup, editors, Proceedings of the 15th Python in Science Conference, pages 98-105, Austin, TX, 2016. SciPy. doi:10.25080/Majora-629e541a-00e

### Remember to run in a folder with your desired PSF and DCD file.


```python
import glob
import MDAnalysis as mda
import prolif as plf
import pandas as pd
import numpy as np
import networkx as nx
from pyvis.network import Network
from matplotlib import cm, colors
from IPython.display import IFrame
from MDAnalysis.topology.guessers import guess_types
```

### Making a Universe

MDAnalysis requires a topology file (such as a psf) and a trajectory file (like a dcd). The cells below will display the current working directory and the files inside. Please make sure your desired topology and trajectory files are inside. 


```python
%pwd

```




    '/mnt/c/users/rosha/desktop/MDContact'




```python
%ls
```

    [0m[01;32mUntitled.ipynb[0m*  [01;32mfinal_eq.dcd[0m*  [34;42mprodigy[0m/  [01;32msystem_autopsf.psf[0m*



```python
top = input('Please input the filename of your topology file with its specified file tag. For example, a psf called "file" should be exactly written as "file.psf"')
```

    Please input the filename of your topology file with its specified file tag. For example, a psf called "file" should be exactly written as "file.psf"system_autopsf.psf



```python
dcd = input('Please input the filename of your trajectory file with its specified file tag. For example, a dcd called "traj" should be exactly written as "traj.dcd"')
```

    Please input the filename of your trajectory file with its specified file tag. For example, a dcd called "traj" should be exactly written as "traj.dcd"final_eq.dcd



```python
# load trajectory
u = mda.Universe(top , dcd)
guessed_elements = guess_types(u.atoms.names)
    
u.add_TopologyAttr('elements', guessed_elements)
print(u.atoms.elements)  # returns an array of guessed elements
```

    ['N' 'H' 'H' ... 'O' 'H' 'H']


### Selecting center of interaction

Enter your selection (ex. 'segid A and resid 1:100') and step size (in frames) for reading your trajectory. For information on how to format your selection, use: https://userguide.mdanalysis.org/1.0.0/selections.html.


```python
sele = input('Please input your selection')

```

    Please input your selectionsegid XP1 and resid 1:100



```python
itr = int(input('Choose the step size for reading your trajectory in frames. For example, inputting "1" will read your trajectory every 1 frame.'))

```

    Choose the step size for reading your trajectory in frames. For example, inputting "1" will read your trajectory every 1 frame.2



```python
selection = u.select_atoms(""+sele+"")
print(selection)
prot = u.select_atoms("protein and not group selection", selection=selection)
fp = plf.Fingerprint()
### The "::1" means that the dcd trajectory will be read every 2 frames.
fp.run(u.trajectory[::itr], selection, prot)
df = fp.to_dataframe()
df.head()
```

    <AtomGroup [<Atom 1: N of type NH3 of resname MET, resid 17 and segid XP1>, <Atom 2: HT1 of type HC of resname MET, resid 17 and segid XP1>, <Atom 3: HT2 of type HC of resname MET, resid 17 and segid XP1>, ..., <Atom 1252: HZ3 of type HC of resname LYS, resid 100 and segid XP1>, <Atom 1253: C of type C of resname LYS, resid 100 and segid XP1>, <Atom 1254: O of type O of resname LYS, resid 100 and segid XP1>]>



      0%|          | 0/100 [00:00<?, ?it/s]





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead tr th {
        text-align: left;
    }

    .dataframe thead tr:last-of-type th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th>ligand</th>
      <th colspan="4" halign="left">MET17</th>
      <th colspan="5" halign="left">ILE18</th>
      <th>TRP19</th>
      <th>...</th>
      <th colspan="2" halign="left">ASN79</th>
      <th colspan="3" halign="left">ASP84</th>
      <th>LEU95</th>
      <th>GLY99</th>
      <th colspan="3" halign="left">LYS100</th>
    </tr>
    <tr>
      <th>protein</th>
      <th colspan="2" halign="left">LYS109</th>
      <th colspan="2" halign="left">LYS110</th>
      <th>LYS109</th>
      <th colspan="4" halign="left">LYS110</th>
      <th>ARG107</th>
      <th>...</th>
      <th colspan="2" halign="left">LYS103</th>
      <th colspan="3" halign="left">ARG107</th>
      <th>ILE104</th>
      <th>ASN101</th>
      <th colspan="3" halign="left">ASN101</th>
    </tr>
    <tr>
      <th>interaction</th>
      <th>Hydrophobic</th>
      <th>VdWContact</th>
      <th>Cationic</th>
      <th>VdWContact</th>
      <th>VdWContact</th>
      <th>Hydrophobic</th>
      <th>HBDonor</th>
      <th>HBAcceptor</th>
      <th>VdWContact</th>
      <th>Hydrophobic</th>
      <th>...</th>
      <th>HBAcceptor</th>
      <th>VdWContact</th>
      <th>HBAcceptor</th>
      <th>Anionic</th>
      <th>VdWContact</th>
      <th>Hydrophobic</th>
      <th>VdWContact</th>
      <th>Hydrophobic</th>
      <th>HBAcceptor</th>
      <th>VdWContact</th>
    </tr>
    <tr>
      <th>Frame</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>...</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>...</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>4</th>
      <td>True</td>
      <td>False</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>...</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>6</th>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>True</td>
      <td>True</td>
      <td>...</td>
      <td>False</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>8</th>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>...</td>
      <td>False</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 154 columns</p>
</div>



### Produce Residue Interaction Map.

The script below from PyVis will now make a residue interaction map. The red-tone circles represent your selected residues, while the blue circles are the other interacting residues. Connecting line widths are calculated based on the amount of occurence of the interaction between two residues.  


```python
### Convert pandas dataframe to a NetworkX object
def make_graph(
    values,
    df=None,
    node_color=["#FFB2AC", "#ACD0FF"],
    node_shape="dot",
    edge_color="#a9a9a9",
    width_multiplier=1,):
    lig_res = values.index.get_level_values("ligand").unique().tolist()
    prot_res = values.index.get_level_values("protein").unique().tolist()

    G = nx.Graph()
    # add nodes
    # https://pyvis.readthedocs.io/en/latest/documentation.html#pyvis.network.Network.add_node
    for res in lig_res:
        G.add_node(
            res, title=res, shape=node_shape, color=node_color[0], dtype="ligand"
        )
    for res in prot_res:
        G.add_node(
            res, title=res, shape=node_shape, color=node_color[1], dtype="protein"
        )

    for resids, value in values.items():
        label = "{} - {}<br>{}".format(
            *resids,
            "<br>".join(
                [
                    f"{k}: {v}"
                    for k, v in (
                        df.xs(resids, level=["ligand", "protein"], axis=1)
                        .sum()
                        .to_dict()
                        .items()
                    )
                ]
            ),
        )
        # https://pyvis.readthedocs.io/en/latest/documentation.html#pyvis.network.Network.add_edge
        G.add_edge(
            *resids,
            title=label,
            color=edge_color,
            weight=value,
            width=value * width_multiplier,
        )

    return G
```


```python
data = (
    df.groupby(level=["ligand", "protein"], axis=1, sort=False)
    .sum()
    .astype(bool)
    .mean()
)

G = make_graph(data, df, width_multiplier=8)

# color each node based on its degree
max_nbr = len(max(G.adj.values(), key=lambda x: len(x)))
blues = cm.get_cmap("Blues", max_nbr)
reds = cm.get_cmap("Reds", max_nbr)
for n, d in G.nodes(data=True):
    n_neighbors = len(G.adj[n])
    # show TM3 in red and the rest of the protein in blue
    palette = reds if d["dtype"] == "ligand" else blues
    d["color"] = colors.to_hex(palette(n_neighbors / max_nbr))

# convert to pyvis network
net = Network(width=640, height=500, notebook=True, heading="")
net.from_nx(G)
net.write_html("prot-prot_graph.html")
IFrame("prot-prot_graph.html", width=650, height=510)
```

    Local cdn resources have problems on chrome/safari when used in jupyter-notebook. 






<iframe
    width="650"
    height="510"
    src="prot-prot_graph.html"
    frameborder="0"
    allowfullscreen

></iframe>





```python
**Note** do the blue circles include residues from the other protein? 
```
