.. _Dadi:

Analyses with dadi
==================

Below we provide code and parameter descriptions for the implementation of
each of the demographic models we tested using dadi. Each of the models assumes
the that the dadi Python module has been imported. All input parameters are
also assumed to be scaled in dadi units (typically a multiple of the reference
population size :math:`N_{ref}`).

**allotetraploid_iso model**

.. code-block:: python
   
   def allotetraploid_iso(params, ns, pts):
       """
       params = (nu,T)
       ns = (n1,n2)
       Split into two populations of specifed size.
       nu: Size of populations after split.
       T: Time in the past of split (in units of 2*Na generations)
       n1,n2: Sample sizes of resulting Spectrum
       pts: Number of grid points to use in integration.
       """
       nu,T = params
       new_ns = [int(ns[0]/2),int(ns[0]/2)]

       xx = dadi.Numerics.default_grid(pts)

       phi = dadi.PhiManip.phi_1D(xx)
       phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

       phi = dadi.Integration.two_pops(phi, xx, T, nu, nu)

       fs_2D = dadi.Spectrum.from_phi(phi, new_ns, (xx,xx), pop_ids=['sub1','sub2'])
       fs_1D = fs_2D.combine_pops([1,2])
       return fs_1D

**allotetraploid_bottleneck model**

.. code-block:: python
   
   def allotetraploid_bottleneck(params, ns, pts):
       """
       params = (nu0,nuBot,T1,T2)
       ns = (n1,n2)
       Split into two populations of specifed size.
       nu0: Size of populations after split.
       nuBot: Size of populations after bottleneck.
       T1: Time in the past between split and polyploid formation
           (in units of 2*Na generations)
       T2: Time since polyploid formation and bottleneck
       n1,n2: Sample sizes of resulting Spectrum
       pts: Number of grid points to use in integration.
       """
       nu0,nuBot,T1,T2 = params
       new_ns = [int(ns[0]/2),int(ns[0]/2)]

       xx = dadi.Numerics.default_grid(pts)

       phi = dadi.PhiManip.phi_1D(xx)
       phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

       phi = dadi.Integration.two_pops(phi, xx, T1, nu0, nu0)
       phi = dadi.Integration.two_pops(phi, xx, T2, nuBot, nuBot)

       fs_2D = dadi.Spectrum.from_phi(phi, new_ns, (xx,xx), pop_ids=['sub1','sub2'])
       fs_1D = fs_2D.combine_pops([1,2])
       return fs_1D

**segtetraploid_iso model**

.. code-block:: python
   
   def segtetraploid_iso(params, ns, pts):
       """
       params = (nu,T1,T2,dij)
       ns = (n1,n2)
       Split into two populations of specifed size.
       nu: Size of populations after split.
       T1: Time in the past between split and polyploid formation
           (in units of 2*Na generations)
       T2: Time since polyploid formation and bottleneck
       dij: homoeologous exchange rate.
       n1,n2: Sample sizes of resulting Spectrum.
       pts: Number of grid points to use in integration.
       """
       nu,T1,T2,dij = params
       new_ns = [int(ns[0]/2),int(ns[0]/2)]

       xx = dadi.Numerics.default_grid(pts)

       phi = dadi.PhiManip.phi_1D(xx)
       phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

       phi = dadi.Integration.two_pops(phi, xx, T1, nu, nu)
       phi = dadi.Integration.two_pops(phi, xx, T2, nu, nu, m12=dij, m21=dij)

       fs_2D = dadi.Spectrum.from_phi(phi, new_ns, (xx,xx), pop_ids=['sub1','sub2'])
       fs_1D = fs_2D.combine_pops([1,2])
       return fs_1D

.. code-block:: python

   def segtetraploid_bottleneck(params, ns, pts):
       """
       params = (nu,nuBot,T1,T2,dij)
       ns = (n1,n2)
       Split into two populations of specifed size.
       nu: Size of populations after split.
       nuBot: Size of populations after bottleneck
       T1: Time in the past between split and polyploid formation
           (in units of 2*Na generations)
       T2: Time since polyploid formation and bottleneck
       n1,n2: Sample sizes of resulting Spectrum
       pts: Number of grid points to use in integration.
       """
       nu,nuBot,T1,T2,dij = params
       new_ns = [int(ns[0]/2),int(ns[0]/2)]

       xx = dadi.Numerics.default_grid(pts)

       phi = dadi.PhiManip.phi_1D(xx)
       phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

       phi = dadi.Integration.two_pops(phi, xx, T1, nu, nu)
       phi = dadi.Integration.two_pops(phi, xx, T2, nuBot, nuBot, m12=dij, m21=dij)

       fs_2D = dadi.Spectrum.from_phi(phi, new_ns, (xx,xx), pop_ids=['sub1','sub2'])
       fs_1D = fs_2D.combine_pops([1,2])
       return fs_1D
