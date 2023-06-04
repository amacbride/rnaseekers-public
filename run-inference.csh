#!/bin/csh -f

if ($#argv != 1) then
	echo "usage: $0 <file of h5ad files>"
	exit -1
endif

set samplefiles=$1

touch out-$$

foreach file (`cat $samplefiles`)
    python infer-2class.py $file | grep -v WARN >> out-$$
end
