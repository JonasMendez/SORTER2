# SORTER2-Toolkit
Sorter of Orthologous Regions for Target Enrichment Reads v2 Toolkit

For detailed documentation on the pipeline and user tutorials please visit www.endemicbio.info/sorter2

We are in the process of publishing SORTER2, but you may cite the following publications if you use SORTER2 prior to publication:

1.Mendez-Reneau, J., Gordon Burleigh, J. and Sigel, E.M., 2023. Target capture methods offer insight into the evolution of rapidly diverged taxa and resolve allopolyploid homeologs in the fern genus Polypodium ss. Systematic Botany, 48(1), pp.96-109.
2.Mendez‐Reneau, J.I., Richards, J.L., Hobbie, J., Bollich, E., Kooyers, N.J. and Sigel, E.M., Lineage diversification and rampant hybridization among subspecies explain taxonomic confusion in the endemic Hawaiian fern Polypodium pellucidum. American Journal of Botany, p.e16379.

SORTER2-Toolkit comprises several python scripts to process paired-end illumina reads for Target-Capture datasets (e.g. Reduced Representation, Probe or Bait-Capture, Target-Enrichment, Hyb-Seq, Amplicon Sequencing, Exome Sequencing, etc..). The pipeline outputs orthologous and phased fasta alignments for non-hybrid organisms as well as alignments including allopolyploid or homoploid hybrid sequences phased into putative progenitor sequence sets, relative to non-hybrid (i.e. 'progenitor') samples. With SORTER2-Processor, users can generate filtered references and get a VCF file amenable for STRUCTURE/ADMIXTURE analyses to help guide the identification of putative progenitor and hybrid samples. In combination with output summary statistics for ortholog copy number (i.e. allopolyploids will have a higher number of sequence copies per ortholog) and other system specific information (i.e. morphology, geography, ecology), SORTER2 provides a robust framework for incorporating hybrid organisms into phylogenetic analyses to identify progenitors for numerous hybrid samples in a single analysis.

SORTER2-Toolkit is an updated version of the previous SORTER pipeline to provide an easier user experience, improve compatibility, as well as provide additional summary statistics to inform hybrid vs non-hybrid hypotheses. Additionally, we have included a new tool (SORTER2_Processor.py) that can: A) Filter output alignments for all stages based on sample representation and number of ortholog clusters per targeted loci (i.e. if a targeted locus resulted in more than one paralogous copies, you can choose how many or these orthologs to retain for downstream analyses), and B) convert the latter filtered orthologous data into consensus references and call SNPs to output a VCF (Variant Call Format) file that can be used for ADMIXTURE or STRUCTURE analyses and help guide the identification of hybrids.

The three overarching goals of the pipeline are to:

1. Build multiple sequence alignments of orthologous sequences for loci generated from paired-end reads. The initial purpose of this pipeline was to reduce paralogous multi-copy sequences and handle heterozygosity by using identity clustering within and among samples, respectively, for the same reference locus. The pipeline generates consensus alleles from heterozygous variants by clustering highly similar contigs (i.e. 99% sequence similarity) associated with the same reference locus and sample to generate consensus sequences presumably representing allelic variation at IUPAC ambiguity sites. The pipeline then maps consensus alleles for all samples to the target references and does a second round of clustering among samples to separate potential paralogs. Given appropriate clustering identity, this can filter and separate paralogs derived from the same locus into separate orthologous sets across all samples, effectively generating additional loci for analysis when paralogs are present. We are in the process of developing a clustering threshold optimizer to assess the best clustering threshold relative to sample representation across clusters.

2. Consensus allele sequences are phased into respective bi-allelic haplotypes for each presumably diploid sample. Phased bi-allelic haplotypes derived from previously determined orthologs are used to build alignments where each sample is represented by two sequences representing heterozygous or homozygous alleles.

​3. Infer hybrid haplotypes based on similarity to potential progenitor samples. This generates a multiple sequence alignment where hybrid samples can have several (phased) tips, corresponding to alternative progenitors, depending on the number of hybrid haplotypes present in the sample. This usually works better for allopolyploid hybrids where we expect low levels of inter-homeologous recombination, but the pipeline has been successful in phasing haplotypes from homoploid hybrids (i.e. Parental 'pseudo-haploblocks' which may have variable levels of diploid recombination between progenitors). The level of recombination in homoploid hybrids and the divergence between their progenitors will influence the relative resolution of the hybrids' phylogenetic relationships. Despite recombination, the pipeline usually generates resolved progenitor hybrid hypotheses at shallow phylogenetic scales (i.e. subspecies) even if they have lower support due to recombination based discordance (Mendez-Reneau et al., 2024)
