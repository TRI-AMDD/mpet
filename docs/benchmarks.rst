Benchmarks
==========
An important aspect of software development is to check the accuracy of results against accepted standards. To this end, mpet has been benchmarked against the following published results in literature, and the mpet configs for reproducing these results are included in the source code repository.


Fuller, Doyle, Newman (1994)
----------------------------
T. F. Fuller, M. Doyle, and J. Newman, `Simulation and Optimization of the Dual Lithium Ion Insertion Cell <https://iopscience.iop.org/article/10.1149/1.2054684>`_, J. Electrochem. Soc. 141, 1 (1994).

MPET config: `configs/params_system_Fuller94.cfg <https://github.com/TRI-AMDD/mpet/blob/master/configs/params_system_Fuller94.cfg>`_

.. image:: benchmarks/Fuller94-Fig2.svg
  :width: 325
  :alt: MPET benchmark, Fuller 1994, Figure 2.
.. image:: benchmarks/Fuller94-Fig3.svg
  :width: 325
  :alt: MPET benchmark, Fuller 1994, Figure 3.

Doyle (1996)
----------------------------
M. Doyle, J. Newman, A. S. Gozdz, C. N. Schmutz, and J.-M. Tarascon, `Comparison of Modeling Predictions with Experimental Data from Plastic Lithium Ion Cells <https://iopscience.iop.org/article/10.1149/1.1836921>`_, J. Electrochem. Soc. 143, 1890 (1996).

MPET config: `configs/params_system_Doyle96-cell1.cfg <https://github.com/TRI-AMDD/mpet/blob/master/configs/params_system_Doyle96-cell1.cfg>`_

.. image:: benchmarks/Doyle96-Fig7.svg
  :width: 325
  :alt: MPET benchmark, Doyle 1996, Figure 7.
.. image:: benchmarks/Doyle96-Fig11.svg
  :width: 325
  :alt: MPET benchmark, Doyle 1996, Figure 11.

LIONSIMBA (2016)
----------------------------
M. Torchio, L. Magni, R. B. Gopaluni, R. D. Braatz, and D. M. Raimondo, `LIONSIMBA: A Matlab Framework Based on a Finite Volume Model Suitable for Li-Ion Battery Design, Simulation, and Control <https://iopscience.iop.org/article/10.1149/2.0291607jes>`_, J. Electrochem. Soc. 163, A1192 (2016).

MPET config: `configs/params_system_LIONSIMBA.cfg <https://github.com/TRI-AMDD/mpet/blob/master/configs/params_system_LIONSIMBA.cfg>`_ 

.. image:: benchmarks/LIONSIMBA-voltage.svg
  :width: 325
  :alt: MPET benchmark, LIONSIMBA, Figure 5a.
.. image:: benchmarks/LIONSIMBA-electrolyte.svg
  :width: 325
  :alt: MPET benchmark, LIONSIMBA, Figure 5b.
.. image:: benchmarks/LIONSIMBA-electrolyte-potential.svg
  :width: 325
  :alt: MPET benchmark, LIONSIMBA, Figure 5c.
.. image:: benchmarks/LIONSIMBA-Cs.svg
  :width: 325
  :alt: MPET benchmark, LIONSIMBA, Figure 5d.