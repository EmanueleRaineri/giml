zcat /home/emanuele/Desktop/meth_data/C000S5A1bs_cpg.txt.gz | awk '$1=="chr1" && $2<10000000' | awk '{if ($2>=10^6 && $2<=10^6+5000) {print $1,$2,0.5,$6+$7} else{print $1,$2,0.9,$6+$7}}' | python sim_meth.py | gzip -c > /home/emanuele/Desktop/meth_data/fake2_cpg.stripped.txt.gz
zcat /home/emanuele/Desktop/meth_data/C001UYA3bs_cpg.txt.gz | awk '$1=="chr1" && $2<10000000' | awk '{if ($2>=10^6 && $2<=10^6+5000) {print $1,$2,0.75,$6+$7} else{print $1,$2,0.9,$6+$7}}' | python sim_meth.py | gzip -c > /home/emanuele/Desktop/meth_data/fake1_cpg.stripped.txt.gz

