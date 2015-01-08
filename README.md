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

A real documentation is not available yet, but here are some usage examples to get you started:

Load the module:

```
julia>using pdbTool
```

Now let us parse the PDB file `5PTI.pdb`. You can find it in the `./testall`
directory of the repository. Exchange `REPO_DIR` with the directory of the
repository in the following command and run it:

```
julia>pdb=pdbTool.parsePdb("REPO_DIR/testall/5PTI.pdb");
```

The object `pdb` contains collections of `Chain` objects, which contain
collections of `Residue` objects, which contain collections of `Atom` objects.
Notice the `;` at the end of the parsing command - if you forget it every object
will be printed (may take minutes) because `pdbTool` does not contain good object printing
yet.

To get the `Atom` object labeled with "N" in residue "10" of chain "A" type

```
julia>pdb.chain["A"].residue["10"].atom["N"]
pdbTool.Atom(154,(30.651,4.022,2.2))
```

Objects of type `Atom` have two values. The first, here 154, is the numerical identifier of the atom in the PDB file. The second are its coordinates. To get them as a tuple write
```
julia>pdb.chain["A"].resdiue["10"].atom["N"].coordinates
(30.651,4.022,2.2)
```

Notice that the key for the residue is a `String` while the identifier for the
atom is a `Int64`. This is because in PDB files consecutive residues are often
labeled e.g. with `10` and `10A`.

To calculate the distance between the atom "N" and the atom "CA" in the same residue, type

```
julia>pdbTool.atomDist(pdb.chain["A"].residue["10"].atom["N"],pdb.chain["A"].residue["10"].atom["CA"])
1.465626487206069
```

To get the distance between residue "10" and "20", type
```
julia>pdbTool.residueDist(pdb.chain["A"].residue["10"],pdb.chain["A"].residue["20"])
6.790367515827108
```
where as default the distance between the two closest heavy atoms in the residues is taken.

To map the chain to the Hidden-Markov model `Kunitz_BPTI.hmm` you find in the `testall` directory type

```
julia>pdbTool.mapChainToHmm(pdb.chain["A"],"REPO_DIR/testall/Kunitz_BPTI.hmm")
"REPO_DIR/testall/Kunitz_BPTI.hmm"
```

The residue "10" in chain "A" now has been identified with a position in the Hidden-Markov model:

```
julia>pdb.chain["A"].residue["10"].alignmentPos
7
```

Residue "1" on the other hand could not be mapped and has retained the default `alignmentPos` value of `-1`
```
julia>pdb.chain["A"].residue["1"].alignmentPos
-1
```

After the mapping, chain "A" has a non-empty collection `align` of length 53,
containing `Residue` objects.  This is the inverse information to the
`alignmentPos` values of the residues: The residue corresponding to position 7
in the Hidden-Markov model (or the 7th column of the corresponding multiple
sequence alignment) should be the residue labeled "10". 

The `identifier` field of a `Residue` object gives the label:

```
julia>pdb.chain["A"].align[7].identifier
"10"
```

Notice that residues with `alignmentPos=-1` do not show up in the `align` collection.

To get a distance matrix for the columns of a multiple sequence alignment corresponding to the Hidden-Markov model you now just have to type

```
julia>ind=sort([k for k in keys(pdb.chain["A"].align)]) # Collect indices
julia>align=pdb.chain["A"].align; # Alias "align" 
julia>[pdbTool.residueDist(align[k1],align[k2]) for k1 in ind, k2 in ind] 
53x53 Array{Float64,2}:
 0.0      1.33243  …   3.68234   8.16824
 1.33243  0.0          2.04305   5.68643
 2.60381  1.29087      5.97459   7.68815
 1.99134  2.83832      6.94746  10.9973 
 5.73631  6.12322     10.1668   13.4135 
 ⋮                 ⋱                    
 7.65084  7.58533      3.44798   3.88685
 2.14929  3.52895  …   1.33761   3.46812
 3.68234  2.04305      0.0       1.27819
 8.16824  5.68643      1.27819   0.0  
```





