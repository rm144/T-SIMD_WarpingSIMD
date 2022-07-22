#!/usr/bin/tcsh
# ===========================================================================
# 
# linkit --
# set links to the code in the PROG path / copy / remove / list / revision
# 
# This file is part of the following software:
# 
#    - the low-level C++ template SIMD library
#    - the SIMD implementation of the MinWarping and the 2D-Warping methods 
#      for local visual homing.
# 
# The software is provided based on the accompanying license agreement
# in the file LICENSE or LICENSE.doc. The software is provided "as is"
# without any warranty by the licensor and without any liability of the
# licensor, and the software may not be distributed by the licensee; see
# the license agreement for details.
# 
# (C) Ralf Möller
#     Computer Engineering
#     Faculty of Technology
#     Bielefeld University
#     www.ti.uni-bielefeld.de
# 
# ===========================================================================

# create tarfile.list
cat simdimage*.files tsimd*.files tiltsearch*.files standalone.files \
    >! tarfile.list

# only tsimd
# create tsimd_tarfile.list
cat tsimd*.files standalone_tsimd.files >! tsimd_tarfile.list

# external file lists
set extlist = ()
set extlist_prefixed = ()
foreach dirfn (`ls -1 *.dir`)
    set extdir = $HOME/`cat $dirfn`
    set filesfn = $dirfn:t:r.files
    # echo "processing $dirfn $filesfn"
    foreach extfn (`cat $filesfn`)
	set extlist = ($extlist $extfn)
	set extlist_prefixed = ($extlist_prefixed $extdir/$extfn)
    end
end
# echo $extlist
# echo $extlist_prefixed

if ($# != 1) then
   echo "linkit remove|link|copy|list|revision"
   exit
endif

if ($1 == "remove") then
   foreach extfile ($extlist)
	echo "removing $extfile"
	/bin/rm -f $extfile
   end
else if ($1 == "link") then
   foreach extfile ($extlist_prefixed)
	echo "linking $extfile"
	/bin/ln -s $extfile .
   end
else if ($1 == "copy") then
   foreach extfile ($extlist_prefixed)
   	echo "copying $extfile"
	/bin/cp $extfile .
   end   
else if ($1 == "list") then
   foreach extfile ($extlist)
   	/bin/ls -lisa $extfile
   end
   echo "---------------------------------------------------------------"
   foreach extfile ($extlist_prefixed)
        /bin/ls -lisa $extfile
   end
else if ($1 == "revision") then
   /bin/rm -f revision
   foreach extfile ($extlist_prefixed)
     svn info $extfile >>! revision
   end
else
  echo "invalid option: $1"
  exit
endif
