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

.. code-block:: bash

  mkdir tmp
  for sample in $(cat NapusSamples.txt)
  do
    parallel-fastq-dump --threads 4 --sra-id ${sample} --tmpdir ./tmp \
      --gzip --split-files --split-technical --read-filter pass \
      --dumpbase --clip

    # CLEAN UP: Run this to remove the cache built by ncbi
    # Make sure you uncomment the line first
    # rm -f $HOME/ncbi/public/sra/${sample}.sra.cache

    # Trim the sequences using sickle
    sickle pe -g -t sanger -f ${sample}_pass_1.fastq.gz -r ${sample}_pass_2.fastq.gz \
      -o trim_${sample}_1.fastq.gz -p trim_${sample}_2.fastq.gz \
      -s singles_${sample}.fastq.gz -q 20 -l 75

    # CLEAN UP: Uncomment the line below to remove the original, untrimmed
    # fastq files and the singletons. You can always run this later.
    # rm -f ${sample}_pass_1.fastq.gz ${sample}_pass_2.fastq.gz singles_${sample}.fastq.gz
  done

  rm -rf tmp

Read Mapping with BWA
---------------------

.. code-block:: bash

  for sample in $(cat NapusSamples.txt)
  do
    bwa mem -t 3 -R "@RG\tID:${sample}\tSM:${sample}" genome/bnapus.fasta \
      fastq/trim_${sample}_1.fastq.gz fastq/trim_${sample}_2.fastq.gz | \
      samtools view -b - > bam/${sample}.raw.bam
  done

Preprocessing BAM Files with Sambamba
-------------------------------------

.. code-block:: bash

  mkdir tmp
  for sample in $(cat NapusSamples.txt)
  do
    sambamba sort -t 3 --tmpdir=./tmp -p -o ${sample}.sorted.bam ${sample}.raw.bam
    sambamba markdup -t 3 --tmpdir=./tmp -p ${sample}.sorted.bam ${sample}.bam
  done

  rm -rf tmp

Variant Calling with GATK HaplotypeCaller
-----------------------------------------

.. code-block:: bash

  for sample in $(cat NapusSamples.txt)
  do
    gatk HaplotypeCaller -R genome/bnapus.fasta \
      -I bam/${sample}.bam -O vcf/${sample}.g.vcf.gz \
      -ERC GVCF
  done


.. code-block:: bash

  for c in $(cat C-chromosomes.txt)
  do
    gatk CombineGVCFs -R genome/bnapus.fasta \
      -L $c -V vcfs_bnapus.list \
      -O vcf/${c}.g.vcf.gz

    gatk GenotypeGVCFs -R genome/bnapus.fasta \
      -L $c -V vcf/${c}.g.vcf.gz \
      -O vcf/${c}.vcf.gz
  done

  for a in $(cat A-chromosomes.txt)
  do
    gatk CombineGVCFs -R genome/bnapus.fasta \
      -L $a -V vcfs_bnapus.list \
      -O vcf/${a}.g.vcf.gz

    gatk GenotypeGVCFs -R genome/bnapus.fasta \
      -L $a -V vcf/${a}.g.vcf.gz \
      -O vcf/${a}.vcf.gz
  done

----

**References**

Wu, D. *et al*. 2019. Whole-Genome Resequencing of a Worldwide Collection of Rapeseed
Accessions Reveals the Genetic Basis of Ecotype Divergence. *Molecular Plant*
12:30--43.
