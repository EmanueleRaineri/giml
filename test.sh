TESTFILE=C001UYA3bs_cpg.chr1.counts.slice.txt
time -p ./giml ${TESTFILE} 1:10:100:1000 2> ${TESTFILE}.log > ${TESTFILE}.giml 
