.. _Simulations:

Simulation-Based Model Validation
=================================

We used `SLiM 3 <https://messerlab.org/slim/>`__ and `SFS_CODE <http://sfscode.sourceforge.net/SFS_CODE/index/index.html>`__
to perform forward simulations.

Autotetraploid Simulations
--------------------------

.. code-block:: bash

  perl optimizeLL.pl 1000000 500 -t 0.005 -r 0.005 -B 0.2 -v

.. code-block:: bash

  for r in {1..50}
  do
    sfs_code 1 100 -t 0.005 -r 0.005 -L 200 5000 \
      -P 4 -p 1 > autotetraploid_${r}.txt
  done

Allotetraploid Simulations
--------------------------

.. code-block:: c

  // set up a simple neutral simulation
  // theta (4 * Na * mu * L) = 10^4
  // L = 10^6 bp
  // Na = 500
  // mu = 5,000 / (4 * 500 * 10^6) = 2.5e-6
  initialize() {
    setwd(".");

    initializeMutationRate(2.5e-6);

    // m1 mutation type: neutral
    initializeMutationType("m1", 0.5, "f", 0.0);

    // g1 genomic element type: uses m1 for all mutations
    initializeGenomicElementType("g1", m1, 1.0);

    // uniform chromosome of length 100 kb with no recombination
    initializeGenomicElement(g1, 0, 999999);
    initializeRecombinationRate(2.5e-6);
  }

  1 {
    sim.addSubpop("p1", 500);
  }

  5000 {
    sim.addSubpopSplit("p2", 500, p1);
    m1.convertToSubstitution=F;
  }

  6000 late() {
    for(i in 1:100){
      pop1_g = sample(p1.genomes, 20, T);
      pop2_g = sample(p2.genomes, 20, T);
      pop1_m = sortBy(unique(pop1_g.mutations), "position");
      pop2_m = sortBy(unique(pop2_g.mutations), "position");
      m = setUnion(pop1_m,pop2_m);
      print(size(m));
      pop1_mutSum = rep(0,size(m));
      pop2_mutSum = rep(0,size(m));
      sfs = rep(0,41);
      for(genome in pop1_g){
        hasMuts = (match(m, genome.mutations) >= 0);
        pop1_mutSum = pop1_mutSum + asInteger(hasMuts);
      }
      for(genome in pop2_g){
        hasMuts = (match(m, genome.mutations) >= 0);
        pop2_mutSum = pop2_mutSum + asInteger(hasMuts);
      }
      mutSum = pop1_mutSum + pop2_mutSum;

      for(mut in mutSum){
        sfs[mut] = sfs[mut] + 1;
      }
      print(sfs);
      writeFile("allotetraploid_" + Rep + "_" + i + ".fs", "41 unfolded\n" + paste(sfs, " ") + "\n");
    }
  }

.. code-block:: bash

  for r in {1..50}; do slim -d Rep=$r allotetraploid.slim; done


Segmental Allotetraploid Simulations
------------------------------------

.. code-block:: c

  // set up a simple neutral simulation
  // theta (4 * Na * mu * L) = 10^4
  // L = 10^6 bp
  // Na = 500
  // mu = 5,000 / (4 * 500 * 10^6) = 2.5e-6
  initialize() {
    setwd(".");

    initializeMutationRate(2.5e-6);

    // m1 mutation type: neutral
    initializeMutationType("m1", 0.5, "f", 0.0);

    // g1 genomic element type: uses m1 for all mutations
    initializeGenomicElementType("g1", m1, 1.0);

    // uniform chromosome of length 100 kb with no recombination
    initializeGenomicElement(g1, 0, 999999);
    initializeRecombinationRate(2.5e-6);
  }

  //...

----

**References**

Haller, B. C. and P. W. Messer. 2019. SLiM 3: Forward genetic simulations beyond the Wright–Fisher model.
*Molecular Biology and Evolution* 36:632–-637.

Hernandez, R. D. 2008. A flexible forward simulator for populations subject to selection and demography.
*Bioinformatics* 24:2786--2787.
