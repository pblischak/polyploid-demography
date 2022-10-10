#!/usr/bin/env python3

import dadi


def main():
    """
    Example for reading in a VCF file to make a SFS.

    For this example, assume we have two populations, pop1 and pop2. Pop1 is
    diploid and pop2 is autotetraploid; both have 10 individuals sampled. We
    have their genotype calls in 'example.vcf.gz'. We also have a mapping file,
    'popinfo.txt', that maps individuals in the VCF to their respective
    populations. To do this we use a two column format for the population info
    mapping: The first column has the names of the individual samples and the
    second column has the name of their corresponding population (see example
    format below).


    Popinfo file
    ------------

    diploid1    pop1
    diploid2    pop1
    diploid3    pop1
    .
    .
    .
    autotet1    pop2
    autotet2    pop2
    autotet3    pop2
    .
    .
    .    


    For the genotypes themselves, dadi can handle parsing polyploid genotype
    calls (e.g., '0/0/0/1', '0/0/1/1') so there's nothing special we need to do
    there.
    """

    # Put the name of your VCF file here
    vcf_file = "example.vcf.gz"

    # Put the name of your population info file here
    popinfo_file = "popinfo.txt"

    # 1. Create the data dictionary
    data_dict = dadi.Misc.make_data_dict_vcf(vcf_file, popinfo_file)

    # 2. Use the data dictionary to create the Spectrum object
    # 
    # For the pop_ids, use the population names from the popinfo file in the same order
    # For the projections, use the number of individuals times the ploidy for each population (# of chromosomes)
    # If you can determine an outgroup allele, you can polarize the direction of mutations from ancestral
    #   to derived and use an unfolded SFS. In this case, set polarized=True. Otherwise, just set polarized
    #   to False use a folded spectrum, which doesn't assume any information about the direction of mutations.
    sfs = dadi.Spectrum.from_data_dict(
        data_dict=data_dict,
        pop_ids=["pop1", "pop2"],
        projections=(20, 40),
        polarized=False
    )

    # If you want to write the SFS to a file, use the `to_file` method
    sfs.to_file("example_sfs.fs")


if __name__ == "__main__":
    main()
