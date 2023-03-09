.. _Conversions:

Converting Parameter Values
===========================

To best understand the relationship between the units used in forward simulators
like SLiM, where we simulate the actual number of individuals, generations, etc.,
the original description of dadi's genetic units in the manual is a good place
to start: `link <https://dadi.readthedocs.io/en/latest/user-guide/specifying-a-model/>`_
(see the subsection titled **Units**). In general, the following conversion rules
described below apply.

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
generations for diploids. So, if :math:`N_{ref}` is 10,000 diploid individuals,
and the generation time is two years, 100,000 years would be 50,000 generations,
and :math:`T` would be 50,000 / (2 * 10,000) = 2.5. For polyploids, the only thing
that changes is the multiplier of :math:`N_{ref}`. For example, if the reference
population is tetraploid, then time is measured as the number of :math:`4N_{ref}`
generations. In general, if the reference population is a K-ploid, then time in
dadi units is measured as the number of :math:`KN_{ref}` generations.

Migration/Exchange rates
------------------------


