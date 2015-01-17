import Base.show

function show(io::IO,pdb::PdbTool.Pdb)
	write(io,"Pdb with $(length(pdb.chain)) chain(s)")
end

function show(io::IO,chain::PdbTool.Chain)
	write(io,"Chain with $(length(chain.residue)) residue(s), $(length(chain.align)) mapped")
end

function show(io::IO,res::PdbTool.Residue)
	write(io,"Residue with $(length(res.atom)) atom(s)")
end
