# ATLAS-toolkit
Simulations toolkit for building molecular simulations systems, applying force field parameters, Generating lamps input files, and trajectory analysis

NOTE:  
Many of the scripts in the `scripts` folder rely on Biograf style atom selection strings. This allows for manipulation and operations across entire structure files. For any script where biograf field/atom selection strings can be used, you can select atoms with the following fields:  

"resid", "resname", "atmname", "index", "x", "y", "z", "fftype", "charge"

the input flag for atom selection string varies, but here are some simple examples of how this would be used

Convert a `psf` file into `bgf` format for use with other ATLAS-toolkit scripts. This only works if you also have a matching `crd` file in the working directory with the same prefix. Otherwise, you should specify the name and path to the desired crd file.
`~/ATLAS-toolkit/scripts/getBGFAtoms.pl -b myStruct.psf -a "index > 0" -s all_my_atoms.bgf`

Remove spurious hydrogen-hydrogen bonds from a BGF generated from CHARMM-GUI outputs using the prior script

`~/ATLAS-toolkit/scripts/removeBond.pl -b step5_assembly.bgf -i "fftype eq 'HT'" -j "fftype eq 'HT'" -s test`

note that strings, like those used to identify an fftype or resid, should be enclosed in single quotes within the selection string, which should go in double quotes.

You can also use regular expressions for this (details coming soon)

