Using `WaveBlocks` for performing simulations
=============================================

In this chapter we show how to use the `WaveBlocks` framework for performing
simulations. The process is always the same and consists of a preprocessing
step, a main step and a post-processing step. The preprocessing step is where
we configure the simulations we want to perform. Then there is the main step
where the simulations are run. Finally, there follows a postprocessing step where
we evaluate the data and (optionally) create visualization. We will see that the
post processing step consists of many small and independent sub-steps reflecting
the various options of what to do with the data obtained.

.. include:: ch1_s1_run_sim.rst
.. include:: ch1_s2_compute_data.rst
.. include:: ch1_s3_visualization.rst
.. include:: ch1_s4_batch_loop.rst
.. include:: ch1_s5_comparing_sims.rst
