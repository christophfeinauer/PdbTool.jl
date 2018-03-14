KRED="\x1B[31m"
KGRN="\x1B[32m"
KRES="\033[0m"

using PdbTool

macro spath()
	return dirname(Base.source_path()) * "/"
end

function TEST_parsePdb()
	@printf("Testing parsePdb ...")
	try
		pdb=PdbTool.parsePdb(@spath() * "5PTI.pdb");
		@printf("%s SUCCES%s\n",KGRN,KRES)
	catch
		@printf("%s FAIL (function exited with error)%s\n", KRED,KRES)
		return false
	end
	return true
end

function TEST_mapChainToHmm()
	@printf("Testing mapChainToHmm ...")
	try
		pdb=PdbTool.parsePdb(@spath() * "5PTI.pdb");
		PdbTool.mapChainToHmm(pdb.chain["A"],@spath() * "Kunitz_BPTI.hmm")
		@printf("%s SUCCESS %s\n",KGRN,KRES)
	catch
		@printf("%s FAIL (function exited with error) %s\n", KRED, KRES);
		return false
	end
	return true
end

function TEST_mapChainToHmmLegacy()
	@printf("Testing mapChainToHmmLegacy ...")
	try
		pdb=PdbTool.parsePdb(@spath() * "5PTI.pdb");
		PdbTool.mapChainToHmmLegacy(pdb.chain["A"],@spath() * "Kunitz_BPTI.hmm")
		@printf("%s SUCCESS\n%s",KGRN,KRES)
	catch
		@printf("%s FAIL (function exited with error) %s\n", KRED, KRES);
		return false
	end
	return true
end

# This function returns always true - the user is just informed that things
# don't work. This is because some tools (cmsearch e.g.) are not necessary for
# most functions and IF necessary the functions will check the external
# functions explicitly
function TEST_externals()

	@printf("Testing access to external functions ... \n")
	@printf("\thmmalign ... ")
	if !PdbTool.EXT_TEST.EXT_TEST_hmmalign()
		@printf("%s NOT AVAILABLE (mapChainToHmm will not work)%s\n", KRED, KRES);
	else 
		@printf("%s AVAILABLE\n%s",KGRN,KRES)
	end
        # FIXME!!!
	# @printf("\thmmsearch ... ")
	# if !PdbTool.EXT_TEST.EXT_TEST_hmmsearch()
	# 	@printf("%s NOT AVAILABLE (mapChainToHmmLegacy will not work)%s\n", KRED, KRES);
	# else 
	# 	@printf("%s AVAILABLE\n%s",KGRN,KRES)
	# end

	@printf("\tcmsearch ... ")
	if !PdbTool.EXT_TEST.EXT_TEST_hmmsearch()
		@printf("%s NOT AVAILABLE (mapping RNAs will not work)%s\n", KRED, KRES);
	else 
		@printf("%s AVAILABLE\n%s",KGRN,KRES)
	end

	return true
end



function TEST_all()

	if !TEST_externals()
		return false
	end
	
	if !TEST_parsePdb()
		return false
	end

	if !TEST_mapChainToHmm()
		return false
	end

	# if !TEST_mapChainToHmmLegacy() #FIXME!
	# 	return false
	# end

	return true
end

function testall()
	if !TEST_all()
		@printf("%s Stopped - if you cannot fix the problem please consider reporting to https://github.com/christophfeinauer/PdbTool/issues %s\n", KRED,KRES);
	else
		@printf("%sAll tests finished with success. %s", KGRN,KRES)
	end
end

testall()
