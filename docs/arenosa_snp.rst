.. _ArenosaSNP:

Variant Calling in *Arabidopsis arenosa*
========================================



Obtaining the Reference Genome
------------------------------

We got the reference genome (v1.0) from Phytozome, along with v2.1 of the gene annotations.

.. code-block:: bash

  # Unzip the file
  gunzip Alyrata_384_v1.fa.gz

  # Change the name to something more manageable
  mv Alyrata_384_v1.fa alyrata.fasta

.. code-block:: bash

  # Index with samtools
  samtools faidx alyrata.fasta

  # Index with BWA
  bwa index alyrata.fasta

  # Create a sequence dictionary with Picard.
  # We'll need this for running GATK. Be sure to substitute
  # the correct path to Picard before running
  java -jar path/to/picard.jar CreateSequenceDictionary \
    R=alyrata.fasta O=alyrata.dict

Getting and Processing FASTQ Files
----------------------------------



Read Mapping with BWA
---------------------



Preprocessing BAM Files with Sambamba
-------------------------------------



Variant Calling with GATK HaplotypeCaller
-----------------------------------------

.. code-block:: bash

  for sample in ArenosaSamples.txt
  do
    gatk HaplotypeCaller -R genome/alyrata.fasta \
      -I bam/${sample}.bam -O vcf/${sample}.g.vcf.gz \
      -ERC GVCF
  done

----

**References**

Monnahan, P. *et al*. 2019. Pervasive population genomic consequences of genome
duplication in *Arabidopsis arenosa*. *Nature Ecology & Evolution* 3:457--468.
