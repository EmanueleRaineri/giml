
---
output: pdf_document
---

\newcommand{\lik}{\ensuremath{\mathcal{L}}}

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

The complexity of the model is linked to the typical scale of the segments used in the sense that if the typical block is very long one needs few of them to cover the genomic region under consideration.



#the algorithm

The algorithm consists of three parts : a formula to evaluate quickly the goodness of fit of any given segmentation
and a protocol to decide (greedily) wich close by segments to join depending on a coarseness parameter $\lambda$. This second
part is an implementation of the algorithm described in the VEGA paper where it is used for copy number variation detection (which is in turn an adaptation of a 2 dimensional image analysis algorithm due to Mumford and Shah, cite xxx)


One starts with a collection of $N$ methylated positions (CpGs or not CpGs) indexed by $j=1,\dots,N$.. For human chromosome $1$ and at a 
decent coverage one has $N \approx 2E6$. One also has a collection of $N_S$ segments (indexed by $i=1,\dots,N_S$), 
which at the beginning
coincide with the methylated positions (i.e. is a collection of $N_S$ segments of length $1$ with $N_S=N$.).
Last, one needs to assign a value to a parameter, $\lambda$ which controls the complexity of the final segmentation.

the initial segmentation has a certain likelihood, which is obtained by summing the log likelihood of each point. In turn, 
the log likelihood of a single CpG $i$ is evaluated as 

\[\lik_i={\overline{n}_i+n \choose \overline{n}_i}\hat{\theta}_i^{\overline{n}_i}(1-\hat{\theta}_i)^{n_i+\overline{n}_i}\]
\label{loglik}
where

$n_i$ and $\overline{n}_i$ are the converted and unconverted reads respectively and

\[\hat{\theta}_i=\frac{\overline{n}_i}{\overline{n}_i+n_i}\]

is the maximum likelihood estimation for the parameter of the binomial process which generates $\overline{n}_i$.

More in general, the likelihdd of a segment containing of more than a single CpG  (for example a segment which includes all the positions from $j_1$ to $j_2$) is computed as in equation (\ref{loglik}) except for $\hat{\theta}$ which is now:

\[
\hat{\theta}=\sum_{j=j_1}^{j_2} \frac{\overline{n}_j}{n_j+\overline{n}_j}
\]


the total likelihood of the segmentation is given by $\mathcal{L}=\sum_i\mathcal{L}_i$.

Now for each pair of adjacent segments $(i,i+1)$ we can compute the loss in total likelihood that we get if we merge them
in one single segment

\[\lik^{\prime}=\lik-L_i-L_{i+1}+L_{i,i+1}\]

notice that $L \leq l_1+l_2$ hence this always results in a loss in total likelihood. This loss is compensated, though, by the fact the the number of segments decreases by $1$. To take into account this decrease in complexity we consider an adjustet likelihood
$\tilde{\lik}=\lik+\lambda$.

We look for the maximum $\tilde{\lik}$ and if it is positive, we merge the corresponding segments into one segment.

We then repeat this operation until we can't merge any more pair of adjacent segments; in this case we can decide to increment 
$\lambda$ and try again, or to give up and exit the algorithm.

##sparsity of CpGs

methylated sites are irregularly spaced along the genome and the variability introduced by the sequencing process might further increase the distance between adjacent sites.
To avoid building segments which span long regions where no CpGs exist, we can multiply $\lambda$ by a penalty term of the
form $\exp(-\frac{D_{i,i+1}}{D_0})$ where $D_{i,i+1}$ is the distance between the the rigth end of segment $i$ and the left end of segment $i+1$ and $D_0$ is aconstant which can be reasonably set to $1000$. 


#examples

#implementation

#bibliography


