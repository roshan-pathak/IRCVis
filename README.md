# Residue Interaction Network Plotter for MD Simulations

The following notebook is composed of functions from MDAnalysis, ProLif, and Pyvis documentations that are altered to take inputs for your protein of interest. It uses much of the documentation from both modules, but the goal is to give it a more interactive format that can be used on a variety of selection types and files, without having to change it in the actual code. 

**NOTE** Be careful with your selections! When selecting for a specific segid or residue type, make sure you have an intimate knowledge of your own topology and trjaectory files.

A Google Colaboratory version of this notebook can be found here: https://colab.research.google.com/drive/1TuLKJOJS0L2p_WfO-gi6C2ljk6hpZE6u?usp=sharing


### Steps
This notebook works with the following steps:
1. Importing the required modules
2. Making an MDAnalysis Universe (a class that describes everything in your system)
3. Selecting the center of a protein-protein interaction (this can be a whole protein, a portion of it, or any selection of residues/atoms permitted by MDAnalysis).
4. Making a 'fingerprint' of the AtomGroup from Step 3, reading the trajectory for interactions every X frames.
5. Producing a graph with Pyvis.

### Example of a Residue Interaction Map
(https://i.imgur.com/tiL7FmR.png)

*This is an interaction network of all Alanine residues in a sample trajectory/topology pair for reference. The selected residues (center of interaction) are colored in red as nodes. Other interacting residues are colored blue. The degree of interaction over the course of the trajectory for each red node darkens the color. The thicknes off the edges between nodes indicates the number of interactions observed over the course of the trajectory.*

*Citations*:
Bouysset, C., Fiorucci, S. ProLIF: a library to encode molecular interactions as fingerprints.
J Cheminform 13, 72 (2021). https://doi.org/10.1186/s13321-021-00548-6

N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319–2327. doi:10.1002/jcc.21787

R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler, D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein. MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations. In S. Benthall and S. Rostrup, editors, Proceedings of the 15th Python in Science Conference, pages 98-105, Austin, TX, 2016. SciPy. doi:10.25080/Majora-629e541a-00e
