pdbTool: An object-oriented Julia tool to parse PDB files and work with them
=============================================================================

Overview
--------

This is an object-oriented Julia tool to parse PDB files (http://www.rcsb.org/)
and work with them. 

In this moment it can parse coordinate information of amino acid chains, RNA
chains and secondary structure information. All information is stored in
hierarchical Julia types and easily accessible. Some functions are implemented
to work with the data, e.g. calculate residue distances with different distance
measures.

A backmapper is implemented that allows to identify residues in the PDB with
residues in a multiple sequence alignment if the Hidden-Markov model is
available. This identification is stored in the PDB object and can be used to
directly work in the alignment representation, e.g. to calculate residue
distances between alignment positions.

Dependencies
------------

The parser should not have any dependencies, but the backmapper expects to find `hmmalign`
and (for the function mapChainToHmmLegacy if you want to use it) `hmmsearch`
on the path. Furthermore, if you want to backmap RNA chains you need `cmsearch`
as well. 

Everything has been tested with HMMER 3.1b1.

Installation
------------

To install the package, use the command

```
	julia>Pkg.clone("https://github.com/christophfeinauer/pdbTool")
```

Alternatively, clone/download the repository and do a
	
```
	julia>include("REPO_DIR/src/pdbTool.jl")
```

whith REPO_DIR replaced with the direcrory you download the repository to.


Test the installation
---------------------

To test some major functions like the parser, the backmapper and the accessability of external functions, first run 

```
julia>using pdbTool
```

and then

```
julia>pdbTool.testall()
```

Documentation
-------------

Not yet available.



