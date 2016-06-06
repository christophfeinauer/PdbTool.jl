# MIT License (MIT)
# 
# Copyright (c) 2015 Christoph Feinauer
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# For questions and suggestions please use the Github page https://github.com/christophfeinauer/PdbTool

module PdbTool
using Compat

	macro spath()
	        return dirname(Base.source_path()) * "/"
	end

	KRED="\x1B[31m"
	KGRN="\x1B[32m"
	KRES="\033[0m"

	#include(@spath()*"../test/"*"runtests.jl")
	include("EXT_TEST.jl")
	include("aminoAcidDict.jl")
	using .EXT_TEST

	######################################################################	
	# TYPE DEFINITIONS
	######################################################################	
	type Atom
	        identifier::Int64
		@compat coordinates::Tuple{Float64,Float64,Float64}
	end
	type Residue
		aminoAcid::AbstractString
	        atom::Dict{AbstractString,Atom}
		pdbPos::Int64
		alignmentPos::Int64
		identifier::AbstractString
		naccess::Float64
		naccess_rel::Float64
		Residue(aa,pp,id)=new(aa,Dict{AbstractString,Atom}(),pp,-1,id,-1.0,-1.0)
	end
	type Strand
		identifier::Int64
		startRes::AbstractString
		endRes::AbstractString
		sense::Int64
		bondThis::AbstractString
		bondPrev::AbstractString
		Strand(identifier::Int64,startRes::AbstractString,endRes::AbstractString,sense::Int64,bondThis::AbstractString,bondPrev::AbstractString)=new(identifier::Int64,startRes::AbstractString,endRes::AbstractString,sense::Int64,bondThis::AbstractString,bondPrev::AbstractString)
	end
	type Sheet
		identifier::AbstractString
		strand::Dict{Int64,Strand}
		Sheet(identifier::AbstractString)=new(identifier,Dict{Int64,Strand}())
	end
	type Helix
		startRes::AbstractString
		endRes::AbstractString
		identifier::AbstractString
		Helix(sR,eR,id)=new(sR,eR,id)
	end
	type Chain
	        residue::Dict{AbstractString,Residue}
	        length::Int64
		mappedTo::AbstractString
		identifier::AbstractString
		align::Dict{Int64,Residue}
		helix::Dict{AbstractString,Helix}
		sheet::Dict{AbstractString,Sheet}
		isRNA::Bool
	        Chain()=new(Dict{AbstractString,Residue}(),0,"","",Dict{Int64,Residue}(),Dict{AbstractString,Helix}(),Dict{AbstractString,Sheet}(),false)
	end
	type Pdb
	        chain::Dict{AbstractString,Chain}
	        pdbName::AbstractString
	        fileName::AbstractString
	        Pdb()=new(Dict{AbstractString,Chain}(),"","")
	end
	include("show_methods.jl")

	######################################################################	
	# FUNCTION:		 parsePdb             	
	######################################################################	
	function parsePdb(pdbFile::AbstractString)
		!isfile(pdbFile) && error("File not found")
		pdb=Pdb()
		for l in eachline(open(pdbFile))
			if l[1:4]=="ATOM" && (l[17]==' ' || l[17]=='A')
				ch=strip(l[22:22])
				res=strip(l[23:27])
				if !haskey(pdb.chain,ch)
					pdb.chain[ch]=Chain()
					pdb.chain[ch].identifier=ch
					if length(lstrip(l[18:20]))<3
						pdb.chain[ch].isRNA=true
						println("Chain $ch treated as RNA")
					end
				end
				if !haskey(pdb.chain[ch].residue,res)
					pdb.chain[ch].length += 1
					pdb.chain[ch].residue[res]=Residue(l[18:20],pdb.chain[ch].length,res)
				end
				pdb.chain[ch].residue[res].atom[strip(l[13:16])]=Atom(parse(Int,l[7:11]),(parse(Float64,l[31:38]),parse(Float64,l[39:46]),parse(Float64,l[47:54])))
				if pdb.chain[ch].isRNA
					for res in values(pdb.chain[ch].residue)
						res.aminoAcid=lstrip(res.aminoAcid)
					end
				end

				
			end
			if l[1:5]=="HELIX"
				ch=strip(l[20:20])
				hel=strip(l[12:14])
				if !haskey(pdb.chain,ch)
					pdb.chain[ch]=Chain()
					pdb.chain[ch].identifier=ch
					if length(lstrip(l[18:20]))<3
						pdb.chain[ch].isRNA=true
						println("Chain $ch treated as RNA")
					end
				end
				if haskey(pdb.chain[ch].helix,hel)
					error("Found the same helix twice and panicked")
				end
				pdb.chain[ch].helix[hel]=Helix(strip(l[22:26]),strip(l[34:38]),hel)
			end
			if l[1:5]=="SHEET"
				if l[22:22]!=l[33:33]
					warn("PDB contains sheets across several chains - ignoring them for now!")
					continue
				end
				ch=strip(l[22:22])
				sh=strip(l[12:14])
				str=parse(Int64,strip(l[8:10]))
				if !haskey(pdb.chain,ch)
					pdb.chain[ch]=Chain()
					pdb.chain[ch].identifier=ch
					if length(lstrip(l[18:20]))<3
						pdb.chain[ch].isRNA
						println("Chain $ch treated as RNA")
					end
				end
				if !haskey(pdb.chain[ch].sheet,sh)
					pdb.chain[ch].sheet[sh]=Sheet(sh)
				end
				if haskey(pdb.chain[ch].sheet[sh].strand,str)
					error("Found the same strand twice and panicked")
				end
				sense=parse(Int64,l[39:40])
				if sense!=0
					pdb.chain[ch].sheet[sh].strand[str]=Strand(str,strip(l[23:27]),strip(l[34:38]),sense,strip(l[51:55]),strip(l[66:70]))
				else
					pdb.chain[ch].sheet[sh].strand[str]=Strand(str,strip(l[23:27]),strip(l[34:38]),sense,"","")
				end
			end
					
		end
		return pdb;
	end

	######################################################################	
	# FUNCTION:		 atomDist             	
	######################################################################	
	function atomDist(atom1::Atom, atom2::Atom)
			return sqrt((atom1.coordinates[1]-atom2.coordinates[1])^2 + (atom1.coordinates[2] - atom2.coordinates[2])^2 + (atom1.coordinates[3] - atom2.coordinates[3])^2)
	end

	######################################################################	
	# FUNCTION:		 residueDist
	######################################################################	
	function residueDist(res1::Residue, res2::Residue; distType="heavyMin")
		if distType=="heavyMin"
			d=Inf
			for atom1 in values(res1.atom)
				for atom2 in values(res2.atom)	
					d=min(d,atomDist(atom1,atom2))
				end
			end
			return d
		end
		if distType=="ca"
			if haskey(res1.atom,"CA") && haskey(res2.atom,"CA")
				return atomDist(res1.atom["CA"],res2.atom["CA"])
			else
			 	error("error calculating distances: Residues do not have CA entries")
			end
		end
		error("Error calculating distances: supported types are heavyMin or ca")
	end

	######################################################################	
	# FUNCTION:		 chainDist
	######################################################################	
	function chainDist(ch1::Chain, ch2::Chain,out="undef")
		if out=="undef"
			for r1=sort([k for k in keys(ch1.residue)])
				for r2=sort([k for k in keys(ch2.residue)])
					@printf("%s %s %f\n", r1, r2, residueDist(ch1.residue[r1],ch2.residue[r2]))
				end
			end
		else 
			fid=open(out,"w")
			for r1=sort([k for k in keys(ch1.residue)])
				for r2=sort([k for k in keys(ch2.residue)])
					@printf(fid,"%s %s %f\n", r1, r2, residueDist(ch1.residue[r1],ch2.residue[r2]))
				end
			end
		end
	end

	######################################################################	
	# FUNCTION:		 interAlignDist
	######################################################################	
	function interAlignDist(ch1::Chain, ch2::Chain;out="undef")
		if ch1.mappedTo==""
			error("chain 1 not mapped")
		end
		if ch2.mappedTo==""
			error("chain 2 not mapped")
		end
		LENG1=getHmmLength(ch1.mappedTo)
		LENG2=getHmmLength(ch2.mappedTo)
		completeAlign=Dict{Int64,Residue}()
		for i in keys(ch1.align)
			completeAlign[i]=ch1.align[i]
		end
		for i in keys(ch2.align)
			completeAlign[i+LENG1]=ch2.align[i]
		end
		ind=sort([x for x in keys(completeAlign)])
		if out=="undef"
			for i=1:length(ind)
				for j=(i+1):length(ind)
					@printf("%s %s %f\n", ind[i], ind[j], residueDist(completeAlign[ind[i]],completeAlign[ind[j]]))
				end
			end
		else 
			fid=open(out,"w")
			for i=1:length(ind)
				for j=(i+1):length(ind)
					@printf(fid,"%s %s %f\n", ind[i], ind[j], residueDist(completeAlign[ind[i]],completeAlign[ind[j]]))
				end
			end
			close(fid)
		end
	end

	######################################################################	
	# FUNCTION:		 chainSeq
	######################################################################	
	function chainSeq(chain::Chain)
		# Unordered pdbSeq
		if !chain.isRNA
			pdbSeq=[aminoAcidDict[chain.residue[k].aminoAcid] for k in keys(chain.residue)]
		else
			pdbSeq=[chain.residue[k].aminoAcid for k in keys(chain.residue)]
		end

		ind=sortperm([chain.residue[k].pdbPos for k in keys(chain.residue)])
		return join(pdbSeq[ind])
	end

	
	######################################################################	
	# FUNCTION:		 mapChainToHmm
	######################################################################	
	function mapChainToHmm(chain::Chain,hmmFile::AbstractString)
		# Check if the chain already has a mapping - and delete it if yes
		if chain.mappedTo!=""
			println("chain $(chain.identifier) already has a mapping.")
		end
		
		# Actual mapping
		pdbSeq=chainSeq(chain)
		tempFile=tempname()
		fid=open(tempFile,"w")
		@printf(fid,">temp\n%s",pdbSeq)
		close(fid)
		if !chain.isRNA
			if !EXT_TEST_hmmalign()
				error("cannot run hmmalign - please check that it is on the path")
			end
                    @compat run(pipeline(`hmmalign $hmmFile $tempFile`,"$tempFile.out"))
		else
			if !EXT_TEST_cmsearch()
				error("cannot run cmsearch - please check that it is on the path")
			end
                    
                    @compat run(pipeline(`cmsearch -A $tempFile.out $hmmFile $tempFile`,DevNull))
		end

		st2fa("$tempFile.out";oFile=tempFile)
	        @compat align=collect(split(readall(tempFile),'\n'))

		rm("$tempFile")
		rm("$tempFile.out")		
		pdbIndices=find([align[2][x]!='-' for x=1:length(align[2])]) 
		cleanIndices=find(![islower(align[2][x]) for x=1:length(align[2])])

		fakeAlign2pdb=-ones(Int64,length(align[2]))
		@compat fakeAlign2pdb[pdbIndices]=collect(1:length(pdbSeq))
		align2pdb=fakeAlign2pdb[cleanIndices]

		for k in keys(chain.residue)
			if chain.residue[k].pdbPos > 0
				x=find(align2pdb.==chain.residue[k].pdbPos)
				if length(x)>1
					error("found several pdb positions")
				end
				if length(x)==1
					chain.residue[k].alignmentPos=x[1]
					chain.align[x[1]]=chain.residue[k]
				end
			end
		end	
		chain.mappedTo=hmmFile
	end

	######################################################################	
	# FUNCTION:		 alignSeqToHmm
	######################################################################	
	function alignSeqToHmm(seq::AbstractString,hmmFile::AbstractString)
		
		# Actual mapping
		tempFile=tempname()
		fid=open(tempFile,"w")
		@printf(fid,">temp\n%s",seq)
		close(fid)
		if !EXT_TEST_hmmalign()
			error("cannot run hmmalign - please check that it is on the path")
		end
                @compat run(pipeline(`hmmalign $hmmFile $tempFile`,"$tempFile.out"))

		st2fa("$tempFile.out";oFile=tempFile)
	        @compat align=collect(split(readall(tempFile),'\n'))

		rm("$tempFile")
		rm("$tempFile.out")		
		return filter(x->!islower(x),align[2])
	end

	######################################################################	
	# FUNCTION:		 mapSeqToHmm
	######################################################################	
	function mapSeqToHmm(seq::AbstractString,hmmFile::AbstractString)
		
		# Actual mapping
		tempFile=tempname()
		fid=open(tempFile,"w")
		@printf(fid,">temp\n%s",seq)
		close(fid)
		if !EXT_TEST_hmmalign()
			error("cannot run hmmalign - please check that it is on the path")
		end
                @compat run(pipeline(`hmmalign $hmmFile $tempFile`,"$tempFile.out"))

		st2fa("$tempFile.out";oFile=tempFile)
	        @compat align=collect(split(readall(tempFile),'\n'))

		rm("$tempFile")
		rm("$tempFile.out")		
		seqIndices=find([align[2][x]!='-' for x=1:length(align[2])]) 
		cleanIndices=find(![islower(align[2][x]) for x=1:length(align[2])])

		fakeAlign2seq=-ones(Int64,length(align[2]))
		@compat fakeAlign2seq[seqIndices]=collect(1:length(seq))
		align2seq=fakeAlign2seq[cleanIndices]
		
		hmm2seq = Dict{Int64,Int64}()
		for (hmm_ind,seq_ind) in enumerate(align2seq)
			seq_ind == -1 && continue
			hmm2seq[hmm_ind] = seq_ind
		end
		return hmm2seq
	end

	######################################################################	
	# FUNCTION:		 mapChainToHmmLegacy
	######################################################################	
	function mapChainToHmmLegacy(chain::Chain,hmmFile::AbstractString)
		# Check if the chain already has a mapping - and delete it if yes
		if chain.mappedTo!=""
			println("chain $(chain.identifier) already has a mapping.")
		end
		
		# Actual mapping
		pdbSeq=chainSeq(chain)
		tempFile=tempname()
		fid=open(tempFile,"w")
		@printf(fid,">temp\n%s",pdbSeq)
		close(fid)
		if !chain.isRNA
			if !EXT_TEST_hmmsearch()
				error("cannot run hmmsearch - please check that it is on the path")
			end
                    @compat run(pipeline(`hmmsearch -A $tempFile.out $hmmFile $tempFile`, DevNull))                    
		else
			if !EXT_TEST_cmsearch()
				error("cannot run cmsearch - please check that it is on the path")
			end
		    @compat run(pipeline(`cmsearch -A $tempFile.out $hmmFile $tempFile`, DevNull))
                    
		end
		st2fa("$tempFile.out";oFile=tempFile)
		align=[split(readall(tempFile),'\n')]	
		rm("$tempFile")
		rm("$tempFile.out")
		(pdbStart,pdbStop)=parse(Int64,matchall(r"\d+",align[1]))
		pdbIndices=find([align[2][x]!='-' for x=1:length(align[2])]) 
		cleanIndices=find(![islower(align[2][x]) for x=1:length(align[2])])
		fakeAlign2pdb=-ones(Int64,length(align[2]))
		fakeAlign2pdb[pdbIndices]=[pdbStart:pdbStop]
		align2pdb=fakeAlign2pdb[cleanIndices]
		for k in keys(chain.residue)
			if chain.residue[k].pdbPos > 0
				x=find(align2pdb.==chain.residue[k].pdbPos)
				if length(x)>1
					error("found several pdb positions")
				end
				if length(x)==1
					chain.residue[k].alignmentPos=x[1]
					chain.align[x[1]]=chain.residue[k]
				end
			end
		end	
		chain.mappedTo=hmmFile
	end

	######################################################################	
	# FUNCTION:		 intraAlignDist
	######################################################################	
	function intraAlignDist(chain::Chain;out="distMat")
		if out=="distMat"
			if chain.mappedTo=="" 
				error("chain has no mapping")
			end
			LENG=getHmmLength(chain.mappedTo)	
			distMat=-ones(LENG,LENG)
			for k1 in keys(chain.residue)
				for k2 in keys(chain.residue)
					if k1==k2
						continue
					end
					r1=chain.residue[k1]
					r2=chain.residue[k2]
					if r1.alignmentPos > 0 && r2.alignmentPos > 0
						distMat[r1.alignmentPos,r2.alignmentPos] = residueDist(r1,r2)
					end
				end
			end
		end
		return distMat
	end

	######################################################################	
	# FUNCTION:		 makeIntraRoc
	######################################################################	
	@compat function makeIntraRoc(score::Array{Tuple{Int64,Int64,Float64},1},chain::Chain;sz=200,cutoff::Float64=8.0,out::AbstractString="return",pymolMode::Bool=false,minSeparation::Int64=4)
		if chain.mappedTo==""
			error("chain has no mapping")
		end
		if out!="return"
			fid=open(out)
		end
		if !pymolMode
			roc=Array(Tuple{AbstractString,AbstractString,Float64,Float64},0)
		else
			roc=Array(Tuple{AbstractString,AbstractString,Int64,Float64},0)
		end
			s::Int64=0
			i::Int64=0
			hits::Int64=0
			while s<sz && i<length(score)
				i+=1
				if abs(score[i][1]-score[i][2]) <= minSeparation 
					continue
				end
				if haskey(chain.align,score[i][1]) && haskey(chain.align,score[i][2])
					s+=1
					if residueDist(chain.align[score[i][1]],chain.align[score[i][2]])<cutoff
						hits+=1
						if pymolMode
							x=1
						else 
							x=hits/s	
						end
					else
						x=hits/s
						if pymolMode
							x=0
						end
					end
					push!(roc,(chain.align[score[i][1]].identifier,chain.align[score[i][2]].identifier,x,score[i][3]))
				end
			end
		return roc
	end

	######################################################################	
	# FUNCTION:		 filterInterScore
	######################################################################	
	@compat function filterInterScore(score::Array{Tuple{Int64,Int64,Float64},1},chain1::Chain,chain2::Chain;sz=200,cutoff::Float64=8.0,out::AbstractString="return")

		# Check if mapping is existent
		if chain1.mappedTo == ""
				error("chain 1 has no mapping")
		elseif chain2.mappedTo == ""
				error("chain 2 has no mapping")
		end
		LENG1=getHmmLength(chain1.mappedTo)
		LENG2=getHmmLength(chain2.mappedTo)

		if out=="return"
			newScore=Array((AbstractString,AbstractString,Float64),sz)
			s::Int64=0
			i::Int64=0
			hits::Int64=0
			positives::Int64=0
			while s<sz && i<size(score,1)
				i+=1
				if score[i][1] <= LENG1 && score[i][2] > LENG1
					ind1=score[i][1]
					ind2=score[i][2]-LENG1
				elseif score[i][2] <= LENG1 && score[i][1] > LENG1
					ind1=score[i][2]
					ind2=score[i][1]-LENG1
				else
					continue
				end
				if haskey(chain1.align,ind1) && haskey(chain2.align,ind2)
					s+=1
					id1=chain1.align[ind1].identifier
					id2=chain2.align[ind2].identifier
					newScore[s]=(id1,id2,score[i][3])
				end
			end
			return newScore
		end
	end
	@compat function filterInterScore(score::Array{Tuple{Int64,Int64,Float64},1},hmm1::AbstractString,hmm2::AbstractString;sz=200,cutoff::Float64=8.0,out::AbstractString="return")

		# Check if mapping is existent
		LENG1=getHmmLength(hmm1)	
		LENG2=getHmmLength(hmm2)

		if out=="return"
			interScore=Array(Tuple{Int64,Int64,Float64},0)
			score1=Array(Tuple{Int64,Int64,Float64},0)
			score2=Array(Tuple{Int64,Int64,Float64},0)

			s::Int64=0
			i::Int64=0
			hits::Int64=0
			positives::Int64=0
			for i=1:length(score)
				
				if score[i][1] <= LENG1 && score[i][2] > LENG1
					ind1=score[i][1]
					ind2=score[i][2]-LENG1
					push!(interScore,(ind1,ind2,score[i][3]))

				elseif score[i][2] <= LENG1 && score[i][1] > LENG1
					ind1=score[i][2]
					ind2=score[i][1]-LENG1
					push!(interScore,(ind1,ind2,score[i][3]))

				elseif score[i][1] <= LENG1 && score[i][2] <= LENG1
					ind1=score[i][1]
					ind2=score[i][2]
					push!(score1,(ind1,ind2,score[i][3]))

				elseif score[i][1] > LENG1 && score[i][2] > LENG1
					ind1=score[i][1]-LENG1
					ind2=score[i][2]-LENG1
					push!(score2,(ind1,ind2,score[i][3]))
				end
			end
			return interScore,score1,score2
		end
	end

	######################################################################	
	# FUNCTION:		 makeInterRoc
	######################################################################	
	@compat function makeInterRoc(score::Array{Tuple{Int64,Int64,Float64},1},chain1::Chain,chain2::Chain;sz=200,cutoff::Float64=8.0,out::AbstractString="return",pymolMode::Bool=false,naccessRatio::Float64=1.0)

		# Check if mapping is existent
		if chain1.mappedTo == ""
				error("chain 1 has no mapping")
		elseif chain2.mappedTo == ""
				error("chain 2 has no mapping")
		end
		LENG1=getHmmLength(chain1.mappedTo)
		LENG2=getHmmLength(chain2.mappedTo)

		# Get naccess cutoff if necessary
		if naccessRatio<1.0
			naList1=zeros(LENG1);
			naList2=zeros(LENG2);
			i=1;
			for r1 in values(chain1.align)
				naList1[i]=r1.naccess;
				i+=1;
			end
			i=1;
			for r2 in values(chain2.align)
				naList2[i]=r2.naccess;
				i+=1;
			end
			naList1=sort(naList1,rev=true);
			na1Cutoff=naList1[parse(Int64,round(LENG1*naccessRatio))];
			println(na1Cutoff)
			naList2=sort(naList2,rev=true);
			na2Cutoff=naList2[parse(Int64,round(LENG2*naccessRatio))];
			println(na2Cutoff)
		end

		if out=="return"
			roc=Array(Tuple{AbstractString,AbstractString,Float64,Float64},0)
			s::Int64=0
			i::Int64=0
			hits::Int64=0
			positives::Int64=0
			while s<sz && i<size(score,1)
				i+=1
				if score[i][1] <= LENG1 && score[i][2] > LENG1
					ind1=score[i][1]
					ind2=score[i][2]-LENG1
				elseif score[i][2] <= LENG1 && score[i][1] > LENG1
					ind1=score[i][2]
					ind2=score[i][1]-LENG1
				else
					continue
				end
				if haskey(chain1.align,ind1) && haskey(chain2.align,ind2)
					if naccessRatio<1.0
						if chain1.align[ind1].naccess<na1Cutoff || chain2.align[ind2].naccess<na2Cutoff
							continue;
						end
					end
					s+=1
					id1=chain1.align[ind1].identifier
					id2=chain2.align[ind2].identifier
					if residueDist(chain1.align[ind1],chain2.align[ind2])<cutoff
						hits+=1
						if pymolMode
							hits=s
						end
						push!(roc,(id1,id2,hits/s,score[i][3]))
					else
						if pymolMode
							hits=0
						end
						push!(roc,(id1,id2,hits/s,score[i][3]))
					end
				end
			end
			return roc
		end
	end

	######################################################################	
	# FUNCTION:		 getHmmLength
	######################################################################	
	function getHmmLength(hmmFile::AbstractString)
		!isfile(hmmFile) && error("File not readable")
		for l in eachline(open(hmmFile))
			if l[1:4] == "LENG" || l[1:4]=="CLEN"
				LENG=parse(Int64,match(r"\d+",l).match)
				return LENG
			end
		end
		error("unable to get hmm length in $hmmFile")
	end

	######################################################################	
	# FUNCTION:		 makeMarriedContactMap
	######################################################################	
	## IDIOCY: THIS DOES THE SAME THING AS THE FUNCTION "interAlignDist"
	function makeMarriedContactMap(chain1::Chain,chain2::Chain;output::AbstractString="default")
		chain1map=chain1.mappedTo
		chain2map=chain2.mappedTo
		chain1map=="" && error("chain $(chain1.identifier) has no mapping")
		chain2map=="" && error("chain $(chain2.identifier) has no mapping")
		LENG1=getHmmLength(chain1map)
		LENG2=getHmmLength(chain2map)
		contactMap=-ones(Float64,LENG1+LENG2,LENG1+LENG2)
		# Protein 1
		for r1 in values(chain1.align)
			for r2 in values(chain1.align)
				r1.alignmentPos==-1 && error("incoherent alignment information")
				r2.alignmentPos==-1 && error("incoherent alignment information")
				contactMap[r1.alignmentPos,r2.alignmentPos]=residueDist(r1,r2)
			end
		end
		# Protein 2
		for r1 in values(chain2.align)
			for r2 in values(chain2.align)
				r1.alignmentPos==-1 && error("incoherent alignment information")
				r2.alignmentPos==-1 && error("incoherent alignment information")
				contactMap[r1.alignmentPos+LENG1,r2.alignmentPos+LENG1]=residueDist(r1,r2)
			end
		end
		# Protein 1 - Protein 2
		for r1 in values(chain1.align)
			for r2 in values(chain2.align)
				r1.alignmentPos==-1 && error("incoherent alignment information")
				r2.alignmentPos==-1 && error("incoherent alignment information")
				d=residueDist(r1,r2)
				contactMap[r1.alignmentPos,r2.alignmentPos+LENG1]=d
				contactMap[r2.alignmentPos+LENG1,r1.alignmentPos]=d
			end
		end
		if output!="default"
			fid=open(output,"w")	
			for i=1:LENG1+LENG2
				for j=(i+1):LENG1+LENG2
					contactMap[i,j]<0.0 && continue
					@printf(fid,"%d\t%d\t%f\n",i,j,contactMap[i,j])
				end
			end
			close(fid)
		end
		return contactMap
		
		
	end

	######################################################################	
	# FUNCTION:		 countContacts
	######################################################################	
	function countContacts(chain::Chain;min_separation=5,cutoff=8.0)
		nums=sort([n for n in keys(chain.align)])
		contacts=0
		for i=1:length(nums)
			for j=(i+1):length(nums)	
				id1=nums[i]
				id2=nums[j]
				if abs(id1-id2)<min_separation
					continue
				end
				res1=chain.align[id1]
				res2=chain.align[id2]
				if res1.alignmentPos<1 || res2.alignmentPos<1
					error("Something went horribly wrong here")
				end
				d=residueDist(res1,res2)
				if d<cutoff
					contacts=contacts+1
				end
			end
		end
		return contacts
	end

	######################################################################	
	# FUNCTION:		 st2fa
	######################################################################	
	function st2fa(iFile;oFile="default") 
		!isfile(iFile) && error("File not found")
		if oFile=="default"
			oFile="$iFile.st2fa"
		end
		iFid=open(iFile,"r")
		oFid=open(oFile,"w")
		seqDict=Dict{AbstractString,AbstractString}()
		for line in eachline(iFid)
			line[1]=='#' && continue
			s=split(line)
			length(s)!=2 && continue
			if haskey(seqDict,s[1])
				seqDict[s[1]]*=s[2]
			else
				seqDict[s[1]]=s[2]
			end
		end
		close(iFid)
		for k in keys(seqDict)
			@printf(oFid,">%s \n%s\n",k,seqDict[k])
		end
		close(oFid)
		return oFile
	end

	######################################################################	
	# FUNCTION:		 fasta2seqDict
	######################################################################	
	function fasta2seqDict(iFile::AbstractString)
		!isfile(iFile) && error("File not found")
		iFid=open(iFile,"r")
		seqDict=Dict{AbstractString,AbstractString}()
		for head in eachline(iFid)
			( head[1]!='>' || eof(iFid) ) && error("Fasta not readable")
			seq=readline(iFid)
			haskey(seqDict,head) && error("Fasta contains several entries for $head")
			seqDict[head]=seq
		end
		close(iFid)
		return seqDict
	end

	######################################################################	
	# FUNCTION:		 interactionSurface(chain1,chain2)
	######################################################################	
	function interactionSurface(chain1,chain2;numbersOnly::Bool=true,alignedOnly::Bool=true,cutoff::Float64=8.0)
		if chain1.mappedTo=="" || chain2.mappedTo==""
			return interactionSurface_unmapped(chain1,chain2;numbersOnly=numbersOnly,cutoff=cutoff)
		end
		pdbPairs=0;
		alignmentPairs=0;
		if !numbersOnly
			iS=Array((AbstractString,AbstractString),0)
		end
		for r1 in values(chain1.residue)
			for r2 in values(chain2.residue)
				if residueDist(r1,r2)<cutoff
					pdbPairs+=1;
					if !numbersOnly && !alignedOnly
						push!(iS,(r1.identifier,r2.identifier))
					end
					if r1.alignmentPos>0 && r2.alignmentPos>0
						alignmentPairs+=1;
						if !numbersOnly && alignedOnly
							push!(iS,(r1.identifier,r2.identifier))
						end
					end

				end
			end
		end
		if numbersOnly
			return pdbPairs,alignmentPairs
		else
			return iS
		end
						
	end
	function interactionSurface_unmapped(chain1,chain2;numbersOnly::Bool=true,cutoff::Float64=8.0)
		pdbPairs=0;
		if !numbersOnly
			iS=Array((String,String),0)
		end
		for r1 in values(chain1.residue)
			for r2 in values(chain2.residue)
				if residueDist(r1,r2)<cutoff
					pdbPairs+=1;
					if !numbersOnly
						push!(iS,(r1.identifier,r2.identifier))
					end
				end
			end
		end
		if numbersOnly
			return pdbPairs
		else
			return iS
		end
						
	end



#</module>
end

