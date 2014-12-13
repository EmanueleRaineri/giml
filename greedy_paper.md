#abstract

Due to the growing numbers of "tracks" of numerical values attached to the human genome,
it is very useful to design algorithms which summarize the data
and extract patterns to help automatizing at least part of the
biological analysis. In this respect, I am looking for ways to compute 
segmentations of  measured methylation values. A segmentation joins adjacent position
which have similar properties, while boundary between segments point to more or less abrupt changes
in the signal which might relate to biological mechanisms. For example
a drop in methylation in the promoter area is is some cases associated with 
transcriptional activity of the corresponding gene; 
segments with intermediate values might be associated to imprinted regions; and so on and so forth.

When it comes to segmenting methylation values there are two salient aspects which cannot be ignored : first, this epigenetic mark seems to have effect at various genomic scales (ranging from hundred of bases to megabases), hence a multiscale method is called for. Second, measurements through whole genome bisulfite sequencing are affected by sampling error which becomes more important at low coverage.

I present a fast probabilistic model which, coupled with a known greedy algorithm for copy number detection (of which I improved the implementation) provides a sampling-error-aware method for producing automatic segmentations at different scales.

I am using this software in my day-to-day work analyzing data for the Blueprint project and will present some practical examples of its utility. 
I am currently extending it so that I can use it to define differentially methylated regions (DMRs).

#introduction

GIMLI is  a software that computes a multiscale segmentation of DNA methylation data acquired through whole genome bisulfite sequencng (EGBS).
A segmentation is a statistical model of a dataset. When building a statistical model one treis to find a reasonable trade off between two 
features of it: its goodness of fit and its complexity.
In our case these two aspects have a very simple interpretation : the goodness of fit increases with the number of segments $N_S$ used to describe the 
methylation data wherease the complexity of the model decreases with $N_S$.  
At one extreme, representing the complete dataset with one big segment would give a very simple model with the worst possible fit and 
at the other using many segments of length $1$ would give us a perfectly fitting complex summary of the data.

The complexity of the model is linked to the typical scale of the segments used 


#the algorithm

#examples

#implementation

#bibliography


