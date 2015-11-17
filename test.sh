rm  -f ~/Desktop/meth_data/089B_cpg.gimli.log 
rm  -f ~/Desktop/meth_data/089B_cpg.gimli
time -p ./gimli <(head ~/Desktop/meth_data/089B_cpg.stripped.txt) 1:10:100:1000 2> ~/Desktop/meth_data/089B_cpg.gimli.log > ~/Desktop/meth_data/089B_cpg.gimli 
