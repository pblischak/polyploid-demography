.. _BnapusSNP:

Variant Calling in *Brassica napus*
===================================



Obtaining the Reference Genome
------------------------------

.. code-block:: bash

  # Obtain the reference sequence from the BRAD database
  wget http://brassicadb.org/brad/datasets/pub/Genomes/Brassica_napus/Brassica_napus_v4.1.chromosomes.fa.gz

  # Unzip the file
  gunzip Brassica_napus_v4.1.chromosomes.fa.gz

  # Change the name to something more manageable
  mv Brassica_napus_v4.1.chromosomes.fa bnapus.fasta

.. code:: bash

  # Index the genome using samtools
  samtools faidx bnapus.fasta

.. code-block:: bash

  # Use samtools to make separate references for the A and C subgenomes.
  # This will also remove the sequences with unknown placement
  # (chromosomes 'chrAnn', 'chrCnn', and 'chrUnn')

  # Get the chromosomes from the A subgenome
  samtools faidx bnapus.fasta \
    $(for i in {1..9}; do printf "chrA0${i} chrA0${i}_random "; done) \
    chrA10 chrA10_random > A_bnapus.fasta

  # Now getting the chromosomes from the C subgenome
  samtools faidx bnapus.fasta \
    $(for i in {1..9}; do printf "chrC0${i} chrC0${i}_random "; done) > C_bnapus.fasta

.. code-block:: bash

  # Now we'll recombine the A and C references to create a joint reference with
  # the unknown sequences removed. Get rid of original reference first:
  rm bnapus.fasta*

  # Now we'll `cat` the A and C references together to create the new total reference
  cat A_bnapus.fasta C_bnapus.fasta > bnapus.fasta

.. code-block:: bash

  # Index everything with samtools
  samtools faidx A_bnapus.fasta
  samtools faidx C_bnapus.fasta
  samtools faidx bnapus.fasta

  # Index everything with BWA
  bwa faidx A_bnapus.fasta
  bwa faidx C_bnapus.fasta
  bwa faidx bnapus.fasta

  # Create a sequence dictionary for everything with Picard.
  # We'll need this for running GATK. Be sure to substitute
  # the correct path to Picard before running
  java -jar path/to/picard.jar CreateSequenceDictionary \
    R=A_bnapus.fasta O=A_bnapus.dict
  java -jar path/to/picard.jar CreateSequenceDictionary \
    R=C_bnapus.fasta O=C_bnapus.dict
  java -jar path/to/picard.jar CreateSequenceDictionary \
    R=bnapus.fasta O=bnapus.dict

Getting and Processing FASTQ Files
----------------------------------

Read Mapping with BWA
---------------------

Preprocessing BAM Files with GATK
---------------------------------

Variant Calling with GATK HaplotypeCaller
-----------------------------------------
