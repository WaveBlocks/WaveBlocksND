Ready made Potentials
---------------------

The following sections contain some potentials that are implemented in the potential
library. The plots show the eigenvalues or energy surfaces. Some potentials
have additional parameters, the default values for these are also Name.


Potential ``cos_osc``
^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = a \left(- \cos{\left (b x \right )} + 1\right)`

* Default values:

  * :math:`a = 0.07`
  * :math:`b = 1.0`

.. image:: fig/cos_osc.png

Potential ``cosh_osc``
^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = a \cosh{\left (b x \right )}`

* Default values:

  * :math:`a = 1`
  * :math:`b = 1`

.. image:: fig/cosh_osc.png

Potential ``delta_gap``
^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \left[\begin{smallmatrix}\frac{1}{2} \tanh{\left (x \right )} & \delta\\\delta & - \frac{1}{2} \tanh{\left (x \right )}\end{smallmatrix}\right]`


.. image:: fig/delta_gap.png

Potential ``delta_gap_diag``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \left[\begin{smallmatrix}\sqrt{\delta^{2} + \frac{1}{4} \tanh^{2}{\left (x \right )}} & 0\\0 & - \sqrt{\delta^{2} + \frac{1}{4} \tanh^{2}{\left (x \right )}}\end{smallmatrix}\right]`


.. image:: fig/delta_gap_diag.png

Potential ``double_well``
^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \sigma \left(x^{2} - 1\right)^{2}`

* Default values:

  * :math:`\sigma = 1.0`

.. image:: fig/double_well.png

Potential ``eckart``
^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \frac{\sigma}{\cosh^{2}{\left (\frac{x}{a} \right )}}`

* Default values:

  * :math:`a = 0.944858082316`
  * :math:`\sigma = 0.038088`

.. image:: fig/eckart.png

Potential ``five_quadratic``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \left[\begin{smallmatrix}\frac{\sigma x^{2}}{2} & 0 & 0 & 0 & 0\\0 & \frac{\sigma x^{2}}{2} & 0 & 0 & 0\\0 & 0 & \frac{\sigma x^{2}}{2} & 0 & 0\\0 & 0 & 0 & \frac{\sigma x^{2}}{2} & 0\\0 & 0 & 0 & 0 & \frac{\sigma x^{2}}{2}\end{smallmatrix}\right]`

* Default values:

  * :math:`\sigma = 0.05`

.. image:: fig/five_quadratic.png

Potential ``four_powers``
^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \left[\begin{smallmatrix}\frac{\sigma x^{2}}{2} & 0 & 0 & 0\\0 & \frac{\sigma x^{4}}{4} & 0 & 0\\0 & 0 & \frac{\sigma x^{6}}{6} & 0\\0 & 0 & 0 & \frac{\sigma x^{8}}{8}\end{smallmatrix}\right]`

* Default values:

  * :math:`\sigma = 0.05`

.. image:: fig/four_powers.png

Potential ``four_quadratic``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \left[\begin{smallmatrix}\frac{\sigma x^{2}}{2} & 0 & 0 & 0\\0 & \frac{\sigma x^{2}}{2} & 0 & 0\\0 & 0 & \frac{\sigma x^{2}}{2} & 0\\0 & 0 & 0 & \frac{\sigma x^{2}}{2}\end{smallmatrix}\right]`

* Default values:

  * :math:`\sigma = 0.05`

.. image:: fig/four_quadratic.png

Potential ``free_particle``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = c`

* Default values:

  * :math:`c = 0`

.. image:: fig/free_particle.png

Potential ``kratzer``
^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \frac{b \left(b - 1\right)}{2 x^{2}} + \frac{x^{2}}{2}`

* Default values:

  * :math:`b = 2.0`

.. image:: fig/kratzer.png

Potential ``morse``
^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = D \left(1 - e^{- a \left(x - x_{0}\right)}\right)^{2}`

* Default values:

  * :math:`a = 0.5`
  * :math:`x_{0} = 0.0`
  * :math:`D = 3.0`

.. image:: fig/morse.png

Potential ``morse_zero``
^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = D \left(- 2 e^{- a \left(x - x_{0}\right)} + e^{- 2 a \left(x - x_{0}\right)}\right)`

* Default values:

  * :math:`a = 0.5`
  * :math:`x_{0} = 0.0`
  * :math:`D = 3.0`

.. image:: fig/morse_zero.png

Potential ``morse_zero_2``
^^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = l^{2} \left(e^{- 2 x + 2 x_{0}} - 2 e^{- x + x_{0}}\right)`

* Default values:

  * :math:`x_{0} = 0.0`
  * :math:`l = 1.0`

.. image:: fig/morse_zero_2.png

Potential ``pert_quadratic``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \frac{\delta^{2} x^{2}}{2} + \frac{\sigma x^{2}}{2}`

* Default values:

  * :math:`\sigma = 0.05`
  * :math:`\delta = 0.2`

.. image:: fig/pert_quadratic.png

Potential ``quadratic``
^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \frac{\sigma x^{2}}{2}`

* Default values:

  * :math:`\sigma = 1/2`

.. image:: fig/quadratic.png

Potential ``quartic``
^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \frac{\sigma x^{4}}{4}`

* Default values:

  * :math:`\sigma = 0.05`

.. image:: fig/quartic.png

Potential ``three_levels``
^^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \left[\begin{smallmatrix}\tanh{\left (- \rho + x \right )} + \tanh{\left (\rho + x \right )} & \delta_{1} & \delta_{2}\\\delta_{1} & - \tanh{\left (\rho + x \right )} & 0\\\delta_{2} & 0 & - \tanh{\left (- \rho + x \right )} + 1\end{smallmatrix}\right]`

* Default values:

  * :math:`\rho = 3.0`

.. image:: fig/three_levels.png

Potential ``three_quadratic``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \left[\begin{smallmatrix}\frac{\sigma x^{2}}{2} & 0 & 0\\0 & \frac{\sigma x^{2}}{2} & 0\\0 & 0 & \frac{\sigma x^{2}}{2}\end{smallmatrix}\right]`

* Default values:

  * :math:`\sigma = 0.05`

.. image:: fig/three_quadratic.png

Potential ``two_crossings``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \left[\begin{smallmatrix}\frac{1}{2} \tanh{\left (- \rho + x \right )} \tanh{\left (\rho + x \right )} & \frac{\delta}{2}\\\frac{\delta}{2} & - \frac{1}{2} \tanh{\left (- \rho + x \right )} \tanh{\left (\rho + x \right )}\end{smallmatrix}\right]`

* Default values:

  * :math:`\rho = 3.0`

.. image:: fig/two_crossings.png

Potential ``two_quadratic``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \left[\begin{smallmatrix}\frac{\sigma x^{2}}{2} & 0\\0 & \frac{\sigma x^{2}}{2}\end{smallmatrix}\right]`

* Default values:

  * :math:`\sigma = 0.05`

.. image:: fig/two_quadratic.png

Potential ``two_quartic``
^^^^^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \left[\begin{smallmatrix}\frac{\sigma x^{4}}{4} & 0\\0 & \frac{\sigma x^{4}}{8}\end{smallmatrix}\right]`

* Default values:

  * :math:`\sigma = 1`

.. image:: fig/two_quartic.png

Potential ``v_shape``
^^^^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \frac{1}{2} \sqrt{4 \delta^{2} + \tanh^{2}{\left (x \right )}}`

* Default values:

  * :math:`\delta = 0.2`

.. image:: fig/v_shape.png

Potential ``wall``
^^^^^^^^^^^^^^^^^^

* Formula: :math:`V(x) = \operatorname{atan}{\left (\sigma x \right )} + \frac{\pi}{2}`

* Default values:

  * :math:`\sigma = 10.0`

.. image:: fig/wall.png
