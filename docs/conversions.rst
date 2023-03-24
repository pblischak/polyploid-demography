.. _Conversions:

Converting Parameter Values
===========================

To best understand the relationship between the units used in forward simulators
like SLiM, where we simulate the actual number of individuals, generations, etc.,
the original description of dadi's genetic units in the manual is a good place
to start: `link <https://dadi.readthedocs.io/en/latest/user-guide/specifying-a-model/#units>`_
(see the subsection titled **Units**). In general, the conversion rules described
below apply.

Population sizes
----------------

Population sizes in dadi are estimated as a multiplier of the reference
population size and are typically assigned the parameter :math:`\nu_p`
for population :math:`p`. 

To start, we need to calculate the reference population size (typically the
size of the ancestral population) since dadi scales everything relative to
this quantity. After we infer the parameters of a demographic model with dadi,
we also can calculate :math:`\theta`, which is equal to :math:`4N_{ref}\mu L`
where :math:`\mu` is the mutation rate in *number of expected mutations per generation*
and :math:`L` is the total sequence length.

Times
-----

Times (denoted :math:`T`) in dadi are measured as the number of :math:`2N_{ref}`
generations. This is regardless of the ploidy level. So, if :math:`N_{ref}` is 
10,000 diploid individuals and the generation time is two years, 100,000 years 
would be 50,000 generations and :math:`T` would be 50,000 / (2 * 10,000) = 2.5.

Migration/Exchange rates
------------------------

Exchange rates are calculated analogously to migration rates, but with some additional
scaling to account for potential differences in the ploidy of each subgenome.
For migration rates in diploids, :math:`M_{i \leftarrow j}` represents the population-scaled
influx of alleles from population :math:`j` into population :math:`i` and is the
parameter that dadi reports when estimating migration rates in a demographic model.
To get the per-generation proportion of allelic influx, :math:`m_{i \leftarrow j}`,
we can divide :math:`M_{i \leftarrow j}` by the number of chromosomes in the
reference population, which would be :math:`2N_{ref}`.

For homoeologous exchanges, denoted :math:`E_{i \leftrightarrow j}`, the only
difference in converting from the parameter estimated by dadi to the non-scaled
version is to take the ploidy level of the subgenomes (or populations) into account.
To do this, we include and additional factor of :math:`\frac{2}{k_i}`, where :math:`k_i`
is the ploidy level of subgenome :math:`i`.
This means that homoeologous exchange rates (and migration rates, too) can still be
calculated by dividing :math:`E_{i \leftrightarrow j} = 2N_{ref}e_{i \leftrightarrow j}`
by :math:`2N_{ref}`. However, for populations that are not diploid, the effect of
the allelic influx must be scaled additionally by :math:`\frac{2}{k_i}`. For diploids,
this fraction is simply equal to 1, but for a tetraploid, it would be equal to 0.5,
meaning that the influx of alleles into the tetraploid subgenomes changes the allele
frequency half as much as it would in a diploid. To put it another way, the same
proportion of allelic exchange can be happening, so :math:`e_{i \leftrightarrow j}`
is still the same between two subgenomes, but if one of the subgenomes has a higher
ploidy level, then the influx is scaled by an additional factor of :math:`\frac{2}{k_i}`.
