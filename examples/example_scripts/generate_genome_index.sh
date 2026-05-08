# Adjust genomeSAindexNbases and genomeChrBinNbits to the size of your genome. 
# Quoting from the STAR manual:
# --genomeSAindexNbases
#    default: 14
#    int: length (bases) of the SA pre-indexing string. Typically between 10 and 15.
#    Longer strings will use much more memory, but allow faster searches. For small
#    genomes, the parameter –genomeSAindexNbases must be scaled down to
#    min(14, log2(GenomeLength)/2 - 1)

# --genomeChrBinNbits
#     default: 18
#     int: =log2(chrBin), where chrBin is the size of the bins for genome storage:
#     each chromosome will occupy an integer number of bins. For a genome with
#     large number of contigs, it is recommended to scale this parameter as min(18,
#     log2[max(GenomeLength/NumberOfReferences,ReadLength)]).

STAR \  
--runMode genomeGenerate \  
--genomeDir path/to/output/dir \  
--genomeFastaFiles path/to/genome.fa \  
--sjdbGTFfile path/to/genome_annotations.gtf \  
--runThreadN 1 \  
--genomeSAindexNbases 3 \  
--genomeChrBinNbits 3