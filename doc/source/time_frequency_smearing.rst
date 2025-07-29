Time-frequency smearing
=======================

Visibilities are an average over time and frequency. It is assumed that the phase is constant over the integration interval.
For sources further away from the phase centre this assumption no longer holds.
The varying phase over the integration interval causes decorrelation and hence a lower amplitude.
This effect is direction and baseline dependent.

Approximation
-------------

The computation of the phase change over frequency is relatively straightforward
because in the uv-plane the frequency dependence translates to a scaling,
so the track in the uv plane is a straight line.

The time dependence is more complicated because this is a rotation causing
an arc in the uv-plane

To get a closed form expression for time frequency smearing we make two 
simplifying assumptions:

* Over the integration interval the rotation can be approximated by a straight line,
* The smearing over time and frequency are independent from each other, i.e.
  frequency smearing is not time dependent and time smearing is not frequency dependent.

Derivation
----------

Now define the source direction as vector :math:`\mathbf{d} = [l, m, n]`, and the baseline as vector :math:`\mathbf{bl} = [u, v, w]`.
The change of the baseline over time and frequency is then given as 

.. math::
  \mathbf{delta\_bl\_time} (\lambda) & = \mathbf{bl} (m) \times \mathbf{ncp\_uvw} / speed\_of\_light * frequency * rotation\_angle \\
  \mathbf{delta\_bl\_freq} (\lambda) & = \mathbf{bl} (m) / speed\_of\_light * delta\_frequency

where :math:`\times` denotes the cross product, and

.. math::
  rotation\_angle = T\_int * angular\_speed

.. math::
  angular\_speed = \frac{2 \pi}{length\_of\_sidereal\_day}

and :math:`\mathbf{ncp\_uvw}` is the vector in the direction of the celestial pole, i.e. earth axis in the u,v,w, coordinate system.
The result of the cross product is the change of the baseline when it is rotated about the vector
in the direction of the ncp.

The smearing effect is then approximately given by

.. math::
  smearing\_factor = \operatorname{sinc}( \langle \mathbf{delta\_bl\_time},\mathbf{d} \rangle ) * \operatorname{sinc}( <\mathbf{delta\_bl\_freq}, \mathbf{d}>)

where :math:`\langle a,b \rangle` denotes the inner product.
