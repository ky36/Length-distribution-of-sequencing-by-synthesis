# Length-distribution-of-sequencing-by-synthesis

On this site you can find the source code of the simulation program 
`sbs_ffcm_simul.c` to simulate the distribution of sequence length for 
sequencing by synthesis (SBS), and compare the simulation results 
with the analytical results.  The program is written in C programming 
language. To compile the program, use something like this 
(different platforms/compilers differ slightly):

`gcc -Wall -o sbs_ffcm_simul sbs_ffcm_simul.c`

I also include two compiled binary (executable) files, one for 32-bit computers
(`sbs_ffcm_simul_32`), one for 64-bit computers (`sbs_ffcm_simul_64`).

To use the program:

`./sbs_ffcm_simul -f <flow cycles> -r <repeats> -a <prob of a> -b <prob of b> -c <prob of c> -d <prob of g> -A a0,a1,... -B b0,b1,...  -C c0,c1,... -D d0,d1,...`

The parameters are:

- -f: number of flow cycles
- -r: number of simulations
- -a: nucleotide context probability of "a"
- -b: nucleotide context probability of "b"
- -c: nucleotide context probability of "c"
- -d: nucleotide context probability of "d"
- -A: a list of nucleotide incorporation probabilities of nucleotide "a"
- -B: a list of nucleotide incorporation probabilities of nucleotide "b"
- -C: a list of nucleotide incorporation probabilities of nucleotide "c"
- -D: a list of nucleotide incorporation probabilities of nucleotide "d"

For example,

`./sbs_ffcm_simul -f 100 -r 200000 -a 0.3333 -b 0.0909 -c 0.4329 -d 0.1429 -A .1090909091,.5,.3,.09090909091 -B .3166666667,.25,.3333333333,.1 -C .646031
7460,.1428571429,.1,.1111111111 -D .425,.2,.25,.125`

See the references for the background and theoretical results.

## References
1. Kong Y. Statistical distributions of pyrosequencing. Journal of computational biology: a journal of computational molecular cell biology. 2009;16(1):31-42. Epub 2008/12/17. [doi: 10.1089/cmb.2008.0106](https://www.liebertpub.com/doi/10.1089/cmb.2008.0106) PubMed PMID: 19072582.
2. Kong Y. Statistical distributions of sequencing by synthesis with probabilistic nucleotide incorporation. Journal of computational biology: a journal of computational molecular cell biology. 2009;16(6):817-27. Epub 2009/06/16. [doi: 10.1089/cmb.2008.0215](https://www.liebertpub.com/doi/10.1089/cmb.2008.0215). PubMed PMID: 19522665.
3. Kong Y. Length distribution of sequencing by synthesis: fixed flow cycle model. J Math Biol. 2013;67(2):389-410. Epub 2012/06/13. [doi: 10.1007/s00285-012-0556-3](https://link.springer.com/article/10.1007/s00285-012-0556-3). PubMed PMID: 22689207.
4. Kong Y. Distributions of positive signals in pyrosequencing. J Math Biol. 2014;69(1):39-54. Epub 2013/06/01. [doi: 10.1007/s00285-013-0691-5](https://link.springer.com/article/10.1007/s00285-013-0691-5). PubMed PMID: 23722629; PMCID: PMC3795870.
