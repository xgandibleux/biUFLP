# =============================================================================
# parser.jl
# August 2022
# =============================================================================


# =============================================================================
"""
    getfnames(fdirectory::String)

    collect the un-hidden filenames available in a given directory
"""
function getfnames(fdirectory::String)
    # fdirectory : string := path + name of the directory where are located the datasets

    # save the path to the current directory
    pwdSaved = pwd()

    # setup the current directory to the targeted directory
    cd(joinpath(pwd(),fdirectory)) 

    # the current directory
    verboseDev ? println("pwd = ", pwd()) : nothing

    # get all the files stored in the directory
    allfiles = readdir()

    # vector of booleans indicating the valid files
    flag = trues(size(allfiles))

    k=1  
    for f in allfiles
        # check all the files available in the directory
        if f[1] != '.'
            # not a hidden file => get it
            verboseDev ? println("fname = ", f) : nothing
        else
            # a hidden file => ignore it
            flag[k] = false
        end
        k = k+1
    end

    # filter on the valid files (not hidden files)
    fnames = allfiles[flag]

    # restore the original pwd
    cd(pwdSaved)

    # return the corresponding vector of filenames
    return fnames
end


# =============================================================================
"""
    load2UFLP(fdirectory::String, fname::String)

Load an instance of bi-objective UFLP
-  ↓  fdirectory :  path to the data file
-  ↓  fname : file name of the data
-  ↑ :  data of an instance
"""

function load2UFLP(fdirectory::String, fname::String)

    f=open(fdirectory*"/"*fname)
 
    # read the number of users (nI)
    nI = parse(Int, readline(f) )
    # read the number of services (nJ) 
    nJ = parse(Int, readline(f) )
    # read the following line (separator -> line without information)
    useless = readline(f)

    # assignment costs users-services (matrix nI x nJ) ------------------------
    c1 = Array{Int64,2}(undef,nI,nJ)
    c2 = Array{Int64,2}(undef,nI,nJ)

    # objective 1
    for i=1:nI
        c1[i,:] = parse.(Int64, split(readline(f)) )
    end
    # read the following line (separator -> line without information)
    useless = readline(f)

    # objective 2
    for i=1:nI
        c2[i,:] = parse.(Int64, split(readline(f)) )
    end
    # read the following line (separator -> line without information)
    useless = readline(f)

    # running costs of services (vector nJ) ------------------------------
    r1 = Array{Int64,1}(undef,nJ)
    r2 = Array{Int64,1}(undef,nJ)

    # objective 1
    r1[:] = parse.(Int64, split(readline(f)) )
    # read the following line (separator -> line without information)
    useless = readline(f)

    # objective 2
    r2[:] = parse.(Int64, split(readline(f)) )
   
    close(f)

    return Instance(fname,nI,nJ,c1,c2,r1,r2)

end