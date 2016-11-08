.. _theory:

Background theory
=================

Introduction
------------

.. _mag3d: http://mag3d.readthedocs.io

The MAGES suite of algorithms is based on the mag3d_ suit, developed at the UBC Geophysical Inversion Facility, is used to invert magnetic responses over a three dimensional distribution of anomalous susceptibility. This manual is designed so that geophysicists who are familiar with the magnetic experiment, but who are not necessarily versed in the details of inverse theory, can use the codes and invert their data. In the following, we describe the basics of the algorithm, but readers are referred to :cite:`LiOldenburg10` for an in-depth discussion of various aspects of the algorithm. Note that an understanding of these components is necessary for the user to have a global view of the algorithms and to use the program library to its fullest extent. 

A magnetic experiment involves measuring the anomalous magnetic field
produced by magnetically susceptible materials beneath the surface,
which have been magnetized by the earth's main magnetic field. The
material with susceptibility :math:`\kappa(x,y,z)` is magnetized when
the earth's main field impinges upon the subsurface formation with flux
intensity :math:`\mathbf{B}_o`. The magnetized material gives rise to a
magnetic field, :math:`\mathbf{B}_a`, which is superimposed on the
inducing field to produce a total, or resultant, field. By measuring the
resultant field and removing the inducing field from the measurements
through numerical processing, one obtains the distribution of the
anomalous field due to the susceptible material. In this program library
we make the assumption that no remanent magnetization is present and
restrict our attention to induced magnetization.

The data from a typical magnetic survey are a set of magnetic field
measurements acquired over a 2D grid above the surface or along a number
of boreholes within the volume of interest. These data are first
processed to yield an estimate of the anomalous field due to the
susceptible material in the area. The goal of the magnetic inversion is
to obtain, from the extracted anomaly data, quantitative information
about the distribution of the magnetic susceptibility in the ground.
Thus, it is assumed that the input data to the inversion program is the
extracted residual anomaly and the programs in the library are developed
accordingly.

Forward modelling
-----------------

General formulation
~~~~~~~~~~~~~~~~~~~

For a given inducing field :math:`\mathbf{B}_o`, the magnetization
:math:`\mathbf{J}` depends upon the susceptibility through a
differential equation. In most mineral exploration cases, the actual
susceptibility is very small. Thus, we use the first order
approximation, where the magnetization is proportional to the
susceptibility and is given by the product of susceptibility and the
inducing magnetic field :math:`\mathbf{H}_o`,

.. math::
    \mathbf{J}={\kappa}\mathbf{H}_o,
    :label: magnetization


where :math:`\mathbf{H}_o=\mathbf{B}_o / \mu_o` and :math:`\mu_o` is the free-space magnetic permeability. This ignores the self-demagnetization effect by which the secondary field reduces the total inducing field within the susceptible region and results in a weaker magnetization than that given by equation :eq:`magnetization`. The anomalous field produced by the distribution of magnetization :math:`\mathbf{J}` given by the following integral equation with a dyadic Green's function: 

.. math::
   \mathbf{B}_a(\mathbf{r})=\frac{\mu_0}{4\pi}\int\limits_V \nabla \nabla \frac{1}{\left |\mathbf{r}-\mathbf{r}_o\right |}\cdot\mathbf{J} dv,
   :label: greensf

where :math:`\mathbf{r}` is the position of the observation point and :math:`V` represents the volume of magnetization at its position, :math:`\mathbf{r}_o`. The above equation is valid for observation locations above the earth's surface. We assume the magnetic permeability is :math:`\mu_o` at borehole locations to allow equation :eq:`greensf` to be valid everywhere.  When the susceptibility is constant within a volume of source region, the above equation can be written in matrix form as:

.. math::
    \mathbf{B}_a=\mu_o\left( \begin{array}{ccc}
    T_{11} & T_{12} & T_{13} \\
    T_{21} & T_{22} & T_{23}\\
    T_{31} & T_{32} & T_{33}
    \end{array} \right)\kappa\mathbf{H}_o\equiv \mu_o\kappa\mathbf{T}\mathbf{H}_o.
    :label: bmatrix


The tensor :math:`\mathbf{T}_{ij}` is given by

.. math::
   \mathbf{T}_{ij}=\frac{1}{4\pi}\int\limits_V\frac{\partial}{\partial x_i}\frac{\partial}{\partial x_j}\frac{1}{\left |\mathbf{r}-\mathbf{r}_o\right |}dv, \mbox{  for }i=1,3 ; j=1,3,
   :label: tij


and where :math:`x_1`, :math:`x_2`, and :math:`x_3` represent :math:`x-, y-`, and :math:`z-`\ directions, respectively. The expressions of :math:`\mathbf{T}_{ij}` for a cuboidal source volume can be found in :cite:`Bhattacharyya64` and :cite:`Sharma66`. Since :math:`\mathbf{T}` is symmetric and its trace is equal to :math:`-1` when the observation is inside the cell and is :math:`0` when the observation is outside the cell, only five independent elements need to be calculated.

Once :math:`\mathbf{T}` is formed, the magnetic anomaly :math:`\mathbf{B}_a` and its projection onto any direction of measurement are easily obtained by the inner product with the directional vector. The projection of the field :math:`\mathbf{B}_a` onto different directions yields different anomalies commonly obtained in the magnetic survey. For instance, the vertical anomaly is simply :math:`B_{az}`, the vertical component of :math:`\mathbf{B}_a`, whereas the total field anomaly is, to first order, the projection of :math:`\mathbf{B}_a` onto the direction of the inducing field :math:`\mathbf{B}_o`.

Numerical implementation of forward modelling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We divide the region of interest into a set of 3D prismatic cells by using a 3D orthogonal mesh and assume a constant susceptibility within each cell. By equation eq:`magnetization`, we have a uniform magnetization within each cell and its field anomaly can be calculated using equation :eq:`bmatrix`. The actual anomaly that would be measured at an observation point is the sum of field produced by all cells having a non-zero susceptibility value. The calculation involves the evaluation of equation :eq:`bmatrix` in a 3D rectangular domain define by each cell. The program that performs this calculation is MAGSENES. As input parameters, the coordinates of the observation points and the inclination and declination of the anomaly direction must be specified for each datum. For generality, each component in a multi-component data set is specified as a separate datum with its own location and direction of projection.

Inversion methodology
---------------------

Let the set of extracted anomaly data be :math:`\mathbf{d} = (d_1,d_2,...,d_N)^T` and the magnetic susceptibility of cells in the model be :math:`\mathbf{\kappa} = (\kappa_1,\kappa_2,...,\kappa_M)^T`. The two are related by the sensitivity matrix

.. _sens_:
.. math::
     \mathbf{d}=\mathbf{G}\mathbf{\kappa}.
     :label: sens

The matrix has elements :math:`g_{ij}` which quantify the contribution to the :math:`i^{th}` datum due to a unit susceptibility in the :math:`j^{th}` cell. The program performs the calculation of the sensitivity matrix, which is to be used by the subsequent inversion. The sensitivity matrix provides the forward mapping from the model to the data during the entire inverse process. We will discuss its efficient representation via the wavelet transform in a separate section. 

For the inversion, the first question that arises concerns definition of the "model". We choose susceptibility in unitless SI units, :math:`\kappa`, as the model for since the anomalous field is directly proportional to the susceptibility. The inverse problem is formulated as an optimization problem where a global objective function, :math:`\phi`, is minimized subject to the constraints in equation :eq:`sens`. The global objective functions consists of two components: a model objective function, :math:`\phi_m`, and a data misfit function, :math:`\phi_d`, such that

.. _globphi_:
.. math::
    \begin{aligned}
    \min \phi = \phi_d+\beta\phi_m \\
    \mbox{s. t. } \kappa^l\leq \kappa \leq \kappa^u, \nonumber
    \end{aligned}
    :label: globphi

where :math:`\beta` is a trade off parameter that controls the relative importance of the model smoothness through the model objective function and data misfit function. When the standard deviations of data errors are known, the acceptable misfit is given by the expected value :math:`\phi_d` and we will search for the value of :math:`\beta` via an L-curve criterion :cite:`Hansen00` that produces the expected misfit. Otherwise, a user-defined value is used. Bound are imposed through the projected gradient method so that the recovered model lies between imposed lower (:math:`\kappa^l`) and upper (:math:`\kappa^u`) bounds. Note that negative bounds are acceptible due to the nature of the equivalent source processing and attempting to re-produce the data given a smooth model. 

We next discuss the construction of a model objective function which, when minimized, produces a model that is geophysically interpretable. However, here we stress that the functional is to create a smooth model since the layer of sources being solved for is completely fictional and a means to create a model that re-produces the data. In general, the objective function gives the flexibility to incorporate as little or as much information as possible. At the very minimum, this function drives the solution towards a reference model :math:`\kappa_o` and requires that the model be relatively smooth in the three spatial directions. Here we adopt a right handed Cartesian coordinate system with positive north and positive down. Let the model objective function be

.. _mof:
.. math::
    \phi_m(\kappa) &=& \alpha_s\int\limits_V \left\{\kappa(\mathbf{r}) \right\}^2dv + \alpha_x\int\limits_V \left\{\frac{\partial \kappa(\mathbf{r})}{\partial x}\right\}^2dv + \alpha_y\int\limits_V \left\{\frac{\partial \kappa(\mathbf{r})}{\partial y}\right\}^2dv
    :label: mof

where :math:`\alpha_s`, :math:`\alpha_x`, and :math:`\alpha_y` are coefficients, which affect the relative importance of different components in the objective function. Numerically, the model objective function in equation :eq:`mof` is discretized onto the mesh defining the susceptibility model using a finite difference approximation. This yields:

.. _modobjdiscr_:
.. math::
    \begin{aligned}
    \phi_m(\mathbf{\kappa}) = \mathbf{\kappa}^T(\alpha_s \mathbf{W}_s^T\mathbf{W}_s + \alpha_x \mathbf{W}_x^T\mathbf{W}_x+\alpha_y \mathbf{W}_y^T\mathbf{W}_y)\mathbf{\kappa}, \nonumber\\
    \equiv \mathbf{\kappa}^T\mathbf{W}_m^T\mathbf{W}_m\mathbf{\kappa}, \nonumber\\
    =\left \| \mathbf{W}_m \mathbf{\kappa} \right \|^2,\end{aligned}
    :label: modobjdiscr


where :math:`\mathbf{\kappa}` is an :math:`M`-length vector representing the recovered model. The next step in setting up the inversion is to define a misfit measure. Here we use the :math:`l_2`-norm measure

.. _phid_:
.. math::
    \phi_d = \left\| \mathbf{W}_d(\mathbf{G}\mathbf{\kappa}-\mathbf{d})\right\|^2.
    :label: phid

For the work here, we assume that the contaminating noise on the data is independent and Gaussian with zero mean. Specifying :math:`\mathbf{W}_d` to be a diagonal matrix whose :math:`i^{th}` element is :math:`1/\sigma_i`, where :math:`\sigma_i` is the standard deviation of the :math:`i^{th}` datum makes :math:`\phi_d` a chi-squared distribution with :math:`N` degrees of freedom. The optimal data misfit for data contaminated with independent, Gaussian noise has an expected value of :math:`E[\chi^2]=N`, providing a target misfit for the inversion. We now have the components to solve the inversion as defined in equation :eq:`globphi`.

To solve the optimization problem when constraints are imposed we use the projected gradients method :cite:`CalamaiMore87,Vogel02`. This technique forces the gradient in the Krylov sub-space minimization (in other words a step during the conjugate gradient process) to zero if the proposed step would make a model parameter exceed the bound constraints. The result is a model that reaches the bounds, but does not exceed them. This method is computationally faster than the log-barrier method because (1) model parameters on the bounds are neglected for the next iteration and (2) the log-barrier method requires the calculation of a barrier term. Previous versions of potential-field codes from UBC-GIF used the logarithmic barrier method :cite:`Wright97,NocedalWright99`.

.. _waveletSection:

Wavelet Compression of Sensitivity Matrix
-----------------------------------------

The major obstacle to the solution of a large-scale magnetic processing problem is the CPU time required for the application of the sensitivity matrix to model vectors. The MAGES program library overcomes these difficulties by forming a sparse representation of the sensitivity matrix using a wavelet transform based on compactly supported, orthonormal wavelets. For more details, the users are referred to :cite:`LiOldenburg03,LiOldenburg10`. In the following, we give a brief description of the method necessary for the use of the MAGES library. 

Each row of the sensitivity matrix in a 3D magnetic inversion can be treated as a 2D image and a 2D wavelet transform can be applied to it. By the properties of the wavelet transform, most transform coefficients are nearly or identically zero. When coefficients of small magnitudes are discarded (the process of thresholding), the remaining coefficients still contain much of the necessary information to reconstruct the sensitivity accurately. These retained coefficients form a sparse representation of the sensitivity in the wavelet domain. The need to store only these large coefficients means that the memory requirement is reduced. Further, the multiplication of the sensitivity with a vector can be carried out by a sparse multiplication in the wavelet domain. This greatly reduces the CPU time. Since the matrix-vector multiplication constitutes the core computation of the inversion, the CPU time for the inverse solution is reduced accordingly. The use of this approach increases the size of solvable problems by nearly two orders of magnitude. However, the model subtleties from the wavelet compression will have an impact on the forward modelled data (typically these are not worried about when 3D inversion is performed as the minute differences in the models won't affect intpertation). Therefore, the full sensitivity matrix is also written to file and can be used at the end of the inversion allowing the inversion code to solve the problem almost the entire way in the wavelet domain and then finish with the accuracy needed to create an accurate representation of the data.

Let :math:`\mathbf{G}` be the sensitivity matrix and :math:`\mathcal{W}` be the symbolic matrix-representation of the 3D wavelet transform. Then applying the transform to each row of :math:`\mathbf{G}` and forming a new matrix consisting of rows of transformed sensitivity is equivalent to the following operation:

.. math::
     \widetilde{\mathbf{G}}=\mathbf{G}\mathcal{W}^T,
     :label: senswvt

where :math:`\widetilde{\mathbf{G}}` is the transformed matrix. The thresholding is applied to individual rows of :math:`\mathbf{G}` by the following rule to form the sparse representation :math:`\widetilde{\mathbf{G}}^S`,

.. math::
     \widetilde{g}_{ij}^{s}=\begin{cases} \widetilde{g}_{ij} & \mbox{if } \left|\widetilde{g}_{ij}\right| \geq \delta _i \\
     0 & \mbox{if } \left|\widetilde{g}_{ij}\right| < \delta _i
     \end{cases}, ~~ i=1,\ldots,N,
     :label: elemg


where :math:`\delta _i` is the threshold level, and :math:`\widetilde{g}_{ij}` and :math:`\widetilde{g}_{ij}^{s}` are the elements of :math:`\widetilde{\mathbf{G}}` and :math:`\widetilde{\mathbf{G}}^S`, respectively. The threshold level :math:`\delta _i` are determined according to the allowable error of the reconstructed sensitivity, which is measured by the ratio of norm of the error in each row to the norm of that row, :math:`r_i(\delta_i)`. It can be evaluated directly in the wavelet domain by the following expression:

.. math:: 
    r_i(\delta_i)=\sqrt{\frac{\underset{\left | {\widetilde{g}_{ij}} \right |<\delta_i}\sum{\widetilde{g}_{ij}}^2}{\underset{j}\sum{\widetilde{g}_{ij}^2}}}, ~~i=1,\ldots,N,
    :label: kappai


Here the numerator is the norm of the discarded coefficients and the denominator is the norm of all coefficients. The threshold level :math:`\delta_{i_o}` is calculated on a representative row, :math:`i_o`. This threshold is then used to define a relative threshold :math:`\epsilon =\delta_{i_{o}}/ \underset{j}{\max}\left | {\widetilde{g}_{ij}} \right |`. The absolute threshold level for each row is obtained by

.. math::
     \delta_i = \epsilon \underset{j}{\max}\left | {\widetilde{g}_{ij}} \right|, ~~i=1,\ldots,N.
     :label: deltai

The program that implements this compression procedure is MAGSENES. The user is asked to specify the relative error :math:`r^*` and the program will determine the relative threshold level :math:`\delta_i`. Usually a value of a few percent is appropriate for :math:`r^*`. For experienced users and ones that are re-inverting the data, the program also allows the direct input of the relative threshold level.

