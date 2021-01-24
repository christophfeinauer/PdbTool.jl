module EXT_TEST

export EXT_TEST_hmmalign
export EXT_TEST_hmmsearch
export EXT_TEST_cmsearch

######################################################################
# FUNCTION:              EXT_TEST_hmmalign
######################################################################
function EXT_TEST_hmmalign()
	try
	    run(pipeline(`hmmalign -h`, devnull))
	catch	
	    return false
	end
    return true
end

######################################################################
# FUNCTION:              EXT_TEST_hmmsearch
######################################################################
function EXT_TEST_hmmsearch()
	try
	    run(pipeline(`hmmsearch -h`, devnull))
	catch	
	    return false
	end
    return true
end

######################################################################
# FUNCTION:              EXT_TEST_cmsearch
######################################################################
function EXT_TEST_cmsearch()
	try
	    run(pipeline(`cmsearch -h`, devnull))
	catch	
	    return false
	end
    return true
end

#/module
end
