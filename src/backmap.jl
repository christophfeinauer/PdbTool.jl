using PdbTool

function print_help()
	@printf("\nScript to write out mappings/distances between a HMM and\na PDB file. Note that this is inefficient for\nbatch jobs.\n\n")
	@printf("%sUSAGE%s:  julia backmap.jl [pdb_file] [chain] [hmm_file] [output_file] [options...]\n\n",PdbTool.KRED,PdbTool.KRES)
	@printf("\tThis will parse the PDB file [pdb_file], then map [chain] to\n\tthe Hidden-Markov model in [hmm_file] and write the mapping\n\tto [output_file].\n")
	@printf("\tFormat:\tColumn 1 = HMM position\n\t\tColumn 2 = PDB Residue Name\n\t\tColumn 3 = Amino Acid\n\n")
	@printf("%sOPTIONS%s:\n\n",PdbTool.KGRN,PdbTool.KRES)	
	@printf("\t-d [distance_output_file]\n")
	@printf("\tCalculate distances between mapped residues and write\n\tto [distance_output_file].\n\tFormat:\tColumn 1 = HMM Pos 1\n\t\tColumn 2 = HMM Pos 2\n\t\tColumn 3 = PDB Residue Name 1\n\t\tColumn 4 = PDB Residue Name 2\n\t\tColumn 5 = Distance in Angstrom\n\n")
@printf("\t-dt [distance_type]\n")
	@printf("\tDistance type used with -d option. Should be 'heavyMin' (distance\n\tbetween closest heavy atoms in residues, default) or\n\t'ca' (distance between CA atoms).\n")
end


# Parse Input

if length(ARGS)==1 && ARGS[1]=="-h"
	print_help()
	exit(0)
end

if length(ARGS)<4
	@printf("%sERROR:%s expected at least 4 arguments, %d given\n",PdbTool.KRED,PdbTool.KRES,length(ARGS))
	@printf("Usage: julia backmap.jl [pdb_file] [chain] [hmm_file] [output_file] [options...]\n")
	@printf("For a list of possible options: julia backmap.jl -h\n")
	exit(1)
end


pdb_file=ARGS[1]
chain=ARGS[2]
hmm_file=ARGS[3]
output_file=ARGS[4]

MAKE_DIST=false
dist_type="heavyMin"
dist_output_file=nothing

options=("-d","-dt")
if length(ARGS)>4
	l=5
	while l<=length(ARGS)
		if !(ARGS[l] in options)
			@printf("%sERROR:%s argument \"%s\" not recognized\n",PdbTool.KRED,PdbTool.KRES,ARGS[l])
			exit(1)
		end
		if ARGS[l]=="-d"
			MAKE_DIST=true
			if length(ARGS)==l
				@printf("%sERROR:%s please provide a output file for distances after -d option\n",PdbTool.KRED,PdbTool.KRES)
				exit(1)
			end
			if (ARGS[l+1] in options)
				@printf("%sERROR:%s please provide a output file for distances after -d option\n",PdbTool.KRED,PdbTool.KRES)
				exit(1)
			end
			dist_output_file=ARGS[l+1]
			l+=2
			continue
		end
		if ARGS[l]=="-dt"
			if length(ARGS)==l
				@printf("%sERROR:%s please provide a valid distance measure (ca or heavyMin) after -dt option\n",PdbTool.KRED,PdbTool.KRES)
				exit(1)
			end
			if (ARGS[l+1] in options)
				@printf("%sERROR:%s please provide a valid distance measure (ca or heavyMin) after -dt option\n",PdbTool.KRED,PdbTool.KRES)
				exit(1)
			end
			if !(ARGS[l+1] in ("ca","heaveMin"))
				@printf("%sERROR:%s please provide a valid distance measure (ca or heavyMin) after -dt option\n",PdbTool.KRED,PdbTool.KRES)
				exit(1)
			end
			dist_type=ARGS[l+1]
			l+=2
			continue
		end
		error("Parsing arguments failed")
	end
end

if !isfile(pdb_file)
	@printf("%sERROR: %sPDB file not found\n",PdbTool.KRED,PdbTool.KRES)

	exit(1)
end

if !isfile(hmm_file)
	@printf("%sERROR: %sHMM file not found\n",PdbTool.KRED,PdbTool.KRES)
	exit(1)
end

@printf("Parsing %s ...",basename(pdb_file))
pdb=PdbTool.parsePdb(pdb_file);
@printf("\%sSUCCESS\%s\n",PdbTool.KGRN,PdbTool.KRES)

if !haskey(pdb.chain,chain)
	@printf("%sERROR: %sPDB object has no chain %s\n",PdbTool.KRED,PdbTool.KRES,chain)
	exit(1)
end

@printf("Mapping %s to chain %s ...",basename(hmm_file),chain)
PdbTool.mapChainToHmm(pdb.chain[chain],hmm_file)
if length(pdb.chain[chain].align)==0 
	@printf("\%sERROR\%s: no residue could be mapped\n",PdbTool.KRED,PdbTool.KRES)
	exit(1)
end
@printf("%sSUCCESS%s (mapped %d residues)\n",PdbTool.KGRN,PdbTool.KRES,length(pdb.chain[chain].align))

@printf("Writing mapping to \"%s\" ... ",output_file)
oFid=open(output_file,"w")
for l in sort(collect(keys(pdb.chain[chain].align)))
	id=pdb.chain[chain].align[l].identifier
	aa=pdb.chain[chain].align[l].aminoAcid
	@printf(oFid,"%d\t%s\t%s\n",l,id,aa)
end
close(oFid)
@printf("done.\n")

@printf("Writing distances to \"%s\" ... ",dist_output_file)
oFid=open(dist_output_file,"w")
if MAKE_DIST
	ind=sort(collect(keys(pdb.chain[chain].align)))
	for i=1:length(ind)-1
		for j=(i+1):length(ind)
			ri=pdb.chain[chain].align[ind[i]]
			rj=pdb.chain[chain].align[ind[j]]
			d=PdbTool.residueDist(ri,rj;distType=dist_type)
			@printf(oFid,"%d\t%d\t%s\t%s\t%f\n",ind[i],ind[j],ri.identifier,rj.identifier,d)
		end
	end
end
close(oFid)
@printf("done.\n")

