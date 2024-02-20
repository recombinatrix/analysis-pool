# Generic analysis scripts

Commonly used analysis scripts

## Contact frequency and distance over time

### `contacts_mda/`

`contacts_mda/contacts.py` uses mda to calculate the frequency of contact between two groups of molecules, eg between a ligand the residues of a protein.  Can also calculate pairwise distances between the two groups, such as to determine the distance between gating residues, but this becomes very slow when used with large groups of residues

## Clustering

### `cluster/`

`cluster/cluster.sh` calculates protein and ligand clusters using gromacs tools.  **This is old and needs to be reworked.**  Requires stripped down trajectories fit to the protein backbone. `clust_trj` is a function to generate those, but in future projects I'll just create them in the initial trajectory processing

## RMSD

### `rmsd_gmx`

`rmsd_gmx/rmsd_gmx_simple.sh` is a simple example of how to calculate protein backbone RMSD with gromacs tools

`rmsd_gmx/rmsd_gmx.sh` is a more complex example using gromacs tools to calculate RMSD, and python to process data.  It runs multiple different calculations including protein backbone RMSD, protein subdomain backbone RMSD, and ligand all atom RMSD, then processes the data into a single csv using python scripts called `concat.py`

### `rmsd_py`

`rmsd_py/rmsd.py` is an overly complex, poorly structured example of calculating RMSD with MD Analysis.  **It needs to be reorked** and doesn't reflect my current practices.

## RMSF

### `rmsf_gmx`

`rmsf_gmx/rmsf_gmx_simple.sh` is a simple example of how to calculate RMSF using gromacs tools, also producing a pdb with the RMSF values saved as a β-factor

`rmsf_gmx/rmsf_gmx.sh` is a more complex example using gromacs tools to calculate RMSF, producing data files and pdb files as above.  Uses python to process data.  It runs multiple different calculations including protein RMSF and protein subdomain RMSF, then processes the data into a small number of csv using python scripts called `concat.py`
