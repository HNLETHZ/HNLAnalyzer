#!/bin/bash

# This is the name of the original header file from the MakeClass() function. Change this accordingly..."
	ClassHeaderName=AnalysisTreeClass.h 
    
# Make Backup copy of this header (just in case)
	cp ${ClassHeaderName} ${ClassHeaderName}~
    
# Make temp copy of this header, which will be modified below
	cp ${ClassHeaderName} temp.h
    
# This line in the "ClassHeaderName" will seed the initial iteration of replace statements ("sed" commands)
	initialSeed="<TFile.h>"
	commentBegin="// Begin list of arrays sizes"
	commentEnd="// End list of arrays sizes"
	
	sed 's|'"${initialSeed}"'|'"${initialSeed}"'\n'"${commentBegin}"'\n'"${commentEnd}"'|' < temp.h > temp1.h
	mv temp1.h temp.h

# Make file which contains arrays sizes
	list_of_arrays_sizes=list_of_arrays_sizes.txt
	# Search for all lines in header files with pattern used for arrays size definition
	find ../../ -name "*.h" -exec grep -h -G 'const *int *n' {} > ${list_of_arrays_sizes}.tmp1 \;
	# Remove all commented lines
	grep -h -G -v '// *const *int *n' ${list_of_arrays_sizes}.tmp1 > ${list_of_arrays_sizes}
	rm -rf ${list_of_arrays_sizes}.tmp1

# Read File which contains Max variables
	while read const datatype var equals max ; do
	
		# Insert Max variables into header file "  
		INSERT="${const} ${datatype} ${var} ${equals} ${max}"
		echo "Insert:" $INSERT
		sed 's|'"${commentEnd}"'|'"${INSERT}"'\n'"${commentEnd}"'|' < temp.h > temp1.h
	
		# Remove "Max" from end of variable (e.g. nMuMax --> nMu).  
		var_length=${#var};
		n=`expr ${var_length} - 3`;
		m="${var:0:${n}}";
		# "m" will be piped into "grep" command in order to stream in all lines that have "//[nMu]" appended to them in the .h file
    	
		cp temp1.h temp.h;	
		# Replace array sizes of all arrays that are bound by this variable
		for i in `grep "${m}\]" temp.h` ; do 
			echo ${i} | grep -q ";"; 
			if [ $? -eq 0 ] ; then
				i=$(echo "${i}" | sed "s|\([]]\)|\\\]|g")
				i=$(echo "${i}" | sed "s|\([[]\)|\\\[|g")
				j=$(echo "${i}" | sed "s|\[.*\]|\\[${m}Max\\\]|")
#				echo "Replacing   ${i}   with    ${j}"
				sed "s|${i}|${j}|" temp1.h > temp2.h
				mv temp2.h temp1.h;
			fi;
		done;
		mv temp1.h temp.h 	
	done < ${list_of_arrays_sizes}
	
	mv temp.h ${ClassHeaderName}
