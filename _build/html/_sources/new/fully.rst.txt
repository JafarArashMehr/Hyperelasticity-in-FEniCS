2. Fully Incompressible Model
===========================================
 
2.1 Constitutive Equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the incompressible Neo-Hookean material, the strain energy function is defined as:

.. math:: 
  :name: eq.18 

  W=c_1(I_1-3)


For the uniaxial tension, the deformation gradient tensor is defined as:

.. math::
  :name: eq.19

   F=\begin{bmatrix} \lambda_1 & 0 & 0\\0 & \lambda_2 & 0\\0 & 0 & \lambda_3\end{bmatrix}


Considering the :ref:`equation.18 <eq.18>` and with respect to :ref:`equation.16 <eq.16>` and :ref:`equation.17 <eq.17>` then we can calculate the :math:`\beta_1`, :math:`\beta_2` and :math:`\beta_3`: 


.. math::
  :name: eq.20

   \beta_1 = 0

   \beta_2 = 2J^{-1}c_1

   \beta_3 = 0

In the above equation :math:`c_1= \frac {\mu}{2}` where and :math:`\mu` is the shear modulus. In addition for a fully compressible material, :math:`J= 1` so :math:`\beta_2=\mu` . According to the :ref:`equation.19 <eq.19>` and :ref:`equation.12 <eq.12>` we can calculate the :math:`\textbf{b}` as following: 

.. math::
  :name: eq.21

   \textbf{b}=\begin{bmatrix} {\lambda_1}^2 & 0 & 0\\0 & {\lambda_2}^2  & 0\\0 & 0 & {\lambda_3}^2 \end{bmatrix}



Then we can write the stress as:

.. math:: 
  :name: eq.22

   \begin{bmatrix} \sigma_{11} & 0 & 0\\0 & 0& 0\\0 & 0 & 0\end{bmatrix} = \beta_{2}\begin{bmatrix} \lambda_{1}^2 & 0 & 0\\0 & \lambda_{2}^2& 0\\0 & 0 & \lambda_{3}^2\end{bmatrix}


It results in 2 equations:


.. math:: 
  :name: eq.23

  \sigma_{11}=\beta_2 \lambda_{1}^2

  0=\beta_2 \lambda_{2}^2
	
                                   
By subtracting the above equations:



.. math:: 
  :name: eq.24 

   \sigma_{11}=(\lambda_{1}^2-\lambda_{2}^2)\beta_{2}


It should be noted that in the uniaxial stretch, the stretches in different directions are:


.. math:: 
  :name: eq.25 

   \lambda_1=\lambda\\
   \lambda_2=\lambda_3= \frac{1}{\sqrt{\lambda}}


Then the final form of the stress is presented in this form:

.. math:: 
  :name: eq.26 

   \sigma_{11}=\mu(\lambda^2-\frac{1}{\lambda})


For the Biaxial tension, the deformation gradients tensor is defined as:

.. math:: 
  :name: eq.27
 
   F=\begin{bmatrix} \lambda & 0 & 0\\0 & \lambda & 0\\0 & 0 & \frac {1}{\lambda^2}\end{bmatrix}



As :math:`\beta_1` and :math:`\beta_3` are equal to zero, after calculating the :math:`\textbf{b}` tensor we can write the :ref:`Equation.16 <eq.16>` in this form:



.. math:: 
  :name: eq.28

   \begin{bmatrix} \sigma_{11} & 0 & 0\\0 & \sigma_{22} & 0\\0 & 0 & 0\end{bmatrix}= \beta_1 \begin{bmatrix} \lambda^2 & 0 & 0\\0 & \lambda^2 & 0\\0 & 0 & \frac{1}{\lambda^4}\end{bmatrix} 


It results in 2 equations:

.. math:: 
  :name: eq.29 

   \sigma_{11}= \beta_2 \lambda^2\\
   0= \beta_2 \frac{1}{\lambda^4}


By subtracting the above equations:

.. math:: 
  :name: eq.30

   \sigma_{11}= \mu (\lambda^2-\frac{1}{\lambda^4})


.. note:: **Alternative Way**

   The stress is defined as follows:

   .. math:: 
     :name: eq.31 

      \sigma= \alpha_1 \lambda_i^{2} + \alpha_{-1} \lambda_i^{-2}-p

  
   Where:

   .. math:: 
     :name: eq.32

      \alpha_1= 2 \frac{\partial W}{\partial I_1}\\
      \alpha_{-1}= -2 \frac{\partial W}{\partial I_2}


   The parameter :math:`p` is defined as hydrostatic pressure:

   .. math:: 
     :name: eq.33

      p= \alpha_1 \frac{1}{\lambda_1^{2} \lambda_2^{2}}+\alpha_{-1}\lambda_1^{2} \lambda_2^{2}


   We can find :math:`\alpha_1` and :math:`\alpha_{-1}` : :math:`\alpha_1=2c_1=\mu` and :math:`\alpha_{-1}=0`

   By combining the :ref:`equation.31 <eq.31>` and :ref:`equation.33 <eq.33>`:

   .. math:: 
     :name: eq.34 

      \sigma_{ii}= \mu \lambda_{i}^2-\mu (\frac {1}{\lambda_i^{2}\lambda_i^{2}})


   Which is same as :ref:`equation.30 <eq.30>`.


2.2. Finite Element Implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the fully incompressible Neo-Hookean material the Jacobian of the deformation gradient
tensor is unity (e.g. :math:`J=1`). The strain energy function is defined as follows:

.. math:: 
  :name: eq.35
 
   W= c_1(I_1-3)+p(J-1)


.. note:: In the above equation, the parameter :math:`Lagrange Multiplier` enforcing the fully incompressibility condition

The stress is defined:

.. math:: 
  :name: eq.36 

   \sigma=\alpha_1 b + \alpha_{-1} b^{-1} - p I


Then the stress term is reduced to:

.. math:: 
  :name: eq.37
 
   \sigma=\mu b - p I


.. note:: When we solve for a fully incompressible material, we should define our problem on a mixed space including a scalar space (to solve for the :math:`p`) and a vector space (to solve for the displacement)


The analytical solution for the uniaxial stretch for different shear modulus could be obtained using the code: 

.. code-block:: python

	lamda = [0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.47,5]

	# Shear modulus = 0.5 MPa
	sigma_05 = []
	# Shear modulus = 1.5 MPa
	sigma_15 = []
	# Shear modulus = 3.5 MPa
	sigma_35 = []

	for i in range (len(lamda)):

		a = 0.5 * (pow(lamda[i], 2) - 1. /(lamda[i]))
		b = 1.5 * (pow(lamda[i], 2) - 1. /(lamda[i]))
		c = 3.5 * (pow(lamda[i], 2) - 1. /(lamda[i]))

		sigma_05.append(a)
		sigma_15.append(b)
		sigma_35.append(c)

	import matplotlib.pyplot as plt
	print (sigma_05)
	plt.xlabel(r'$\mathrm{Stretch}$', fontsize=20)
	plt.ylabel(r'$\mathrm{\sigma_{xx}(MPa)}$', fontsize=20)

	plt.plot(lamda,sigma_05,  linestyle='-', linewidth=4, color='maroon',label=r'$(\mu=0.5)$')
	plt.plot(lamda,sigma_15,  linestyle='-', linewidth=4, color='r',label=r'$(\mu=1.5)$')
	plt.plot(lamda,sigma_35,  linestyle='-', linewidth=4, color='teal',label=r'$(\mu=3.5)$')
	lg=plt.legend(ncol=1, loc=2, fontsize=15)
	axes = plt.gca()
	axes.set_xlim([0,5])
	axes.set_ylim([-15, 30])

	axes.set_yticks([-15,-10,-5,0,5,10,15,20,25,30,35])
	axes.set_xticks([0,1,2,3,4,5])

	plt.tick_params(axis='both', which='major', labelsize=15)
	plt.grid()
	plt.show()

The implmentation in FEniCS is presented in this code:  

.. code-block:: python

	from dolfin import *

	# Defining Stretches
	stretch = [0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.47,5]
	BC = []

	for x in range(len(stretch)):
		      N = stretch[x] - 1.0
		      BC.append(N)

	tol = 1E-14
	# Define boundaries
	def FRONT(x, on_boundary):
		      return on_boundary and abs(x[2] - 1.0) < tol


	# Defining the mesh which is a single hexahedron element
	mesh = UnitCubeMesh.create(1,1,1,CellType.Type.hexahedron)

	############################################
	#element for pressure field
	Element1 = FiniteElement("CG", mesh.ufl_cell(), 1)
	#element for displacement field
	Element2 = VectorElement("CG", mesh.ufl_cell(), 2)

	# Defining the mixed function space
	W_elem = MixedElement([Element1, Element2])
	W = FunctionSpace(mesh, W_elem)
	#############################################

	boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
	subdomains = MeshFunction('size_t', mesh, mesh.topology().dim())


	# Defining integration symbols
	dx = Measure('dx', domain=mesh, subdomain_data=subdomains, metadata={'quadrature_degree': 15})
	ds = Measure('ds', domain=mesh, subdomain_data=boundaries, metadata={'quadrature_degree': 15})


	## Define variational problem
	dw = TrialFunction(W)            # Incremental displacement
	v  = TestFunction(W)
	w = Function(W)
	p,u = split(w)

	#Definig some continuum mechanics relations
	d = u.geometric_dimension()
	I = Identity(d)             # Identity tensor
	F = I + grad(u)             # Deformation gradient
	C = F.T*F                   # Right cauchy stress tensor
	b = F*F.T                   # Left cauchy stress tensor
	Ib = tr(b)                  # First invariant of b tensor
	J  = det(F)                 # Jacobian of deformation gradient tensor

	# Shear modulus
	mu = 0.5e6

	# Strain energy function for fully incompressible model
	psi = mu/2.*(Ib - 3.)*dx - p*(J - 1)*dx

	# Gateaux derivative in the direction of the test function
	F1 = derivative(psi, w, TestFunction(W))

	# Compute Jacobian of F
	Jac = derivative(F1, w, TrialFunction(W))

	sigma_11 = []

	def border(x, on_boundary):
		      return on_boundary

	bound_x =  Expression(("t*x[0]"), degree=1, t=0)

	for i in range(len(BC)):

		      bound_x.t = BC[i]

		      bc_x = DirichletBC(W.sub(1).sub(0), bound_x, border)
		      bc_front = DirichletBC(W.sub(1).sub(2), Constant((0)), FRONT)

		      bc_all = [bc_x,bc_front]

		      problem = NonlinearVariationalProblem(F1, w, bc_all, Jac)

		      solver = NonlinearVariationalSolver(problem)

		      solver.solve()

		      (p, u) = w.split(True)

	#Stress calculation
		      sig = mu * b - p * I
	# Defining a tenso function space
		      V = TensorFunctionSpace(mesh, 'Lagrange', 1)
	# Projection of the stress on the tensor function space
		      sig1 = project(sig, V)

		      sigma_11.append((sig1.vector().get_local()[0])*0.000001)

	print (sigma_11)

	#Obtained reaults for mu = 0.5 MPa#
	sigma_05 = [-3.3220833333333144, -2.4800000000022377, -1.968750000000018, -1.6216666666666604, -1.3673214285714232, -1.1699999999999993, -1.0098611111111104, -0.8749999999999992, -0.7578409090909061, -0.6533333334279298, -0.5579807692745585, -0.4692857143072567, -0.3854166666778288, -0.3050000000060496, -0.2269852941210564, -0.15055555555754435, -0.07506578947488027, -7.390848180406038e-13, 0.38124999999999554, 0.7916666666666644, 1.2455357142857126, 1.7499999999999942, 2.309027777777776, 2.924999999888631, 3.5994318181378744, 4.333333333314409, 5.12740384614507, 5.982142857138549, 6.897916666664428, 7.874999999998786, 8.913602941175803, 10.013888888888484, 9.878593176733789, 12.399999999906841]

	#Obtained reaults for mu = 1.5 MPa#
	sigma_15 = [-9.966249999999953, -7.4400000000052655, -5.906250000000058, -4.864999999999997, -4.101964285714497, -3.50999999999999, -3.029583333333331, -2.624999999999993, -2.2735227272727245, -1.9599999999999933, -1.6739423078102886, -1.4078571429217734, -1.1562500000334857, -0.9150000000181479, -0.6809558823631715, -0.45166666667263133, -0.2251973684246401, -2.217263493368278e-12, 1.1437499999999825, 2.3749999999999947, 3.7366071428571317, 5.250000000116358, 6.927083333333332, 8.774999999844844, 10.798295454413642, 12.999999999943254, 15.382211538435248, 17.94642857141561, 20.69374999999328, 23.624999999996415, 26.740808823527377, 30.04166666666549, 29.63577953020131, 37.19999999983032]

	#Obtained reaults for mu = 3.5 MPa#
	sigma_35 = [-23.254583333333265, -17.360000000015262, -13.781250000000021, -11.351666666666649, -9.571249999999983, -8.189999999999996, -7.069027777155056, -6.124999999999978, -5.304886363636268, -4.573333333994748, -3.9058653849219174, -3.285000000150805, -2.697916666744805, -2.1350000000423432, -1.5888970588473983, -1.0538888889028104, -0.5254605263241604, -5.173676343700641e-12, 2.6687499999999695, 5.5416666666666545, 8.718749999880098, 12.249999999999975, 16.163194443414437, 20.47499999922037, 25.19602272696517, 30.333333333200915, 35.89182692301562, 41.874999999969944, 48.28541666665088, 55.1249999999915, 62.39522058823055, 70.09722222221926, 69.15015223713617, 86.79999999993669]

	#Plotting the Results#
	import matplotlib.pyplot as plt
	plt.xlabel(r'$\mathrm{Stretch}$', fontsize=20)
	plt.ylabel(r'$\mathrm{\sigma_{xx}(MPa)}$', fontsize=20)

	plt.plot(stretch,sigma_05,  linestyle='-', linewidth=4, color='y',label=r'$(\mu=0.5\/\/MPa)$')
	plt.plot(stretch,sigma_15,  linestyle='-', linewidth=4, color='c',label=r'$(\mu=1.5\/\/MPa)$')
	plt.plot(stretch,sigma_35,  linestyle='-', linewidth=4, color='k',label=r'$(\mu=3.5\/\/MPa)$')
	lg=plt.legend(ncol=1, loc=2, fontsize=15)
	axes = plt.gca()
	axes.set_xlim([0,5])
	axes.set_ylim([-15, 30])

	axes.set_yticks([-15,-10,-5,0,5,10,15,20,25,30,35])
	axes.set_xticks([0,1,2,3,4,5])

	plt.tick_params(axis='both', which='major', labelsize=15)
	plt.grid()
	plt.show()

The results of analytical solution and finite element implementation are compared in
figure.1:


.. figure:: PNG/2.png
   :align: center
	
   The stress results for the fully incompressible material model in uniaxial loading - Finite Element (Right) vs Analytical (Left)

The analytical solution for the Biaxial stretch for different shear modulus could be obtained using the code: 

.. code-block:: python


	lamda = [0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.47,5]

	# Shear modulus = 0.5 MPa
	sigma_05 = []
	# Shear modulus = 1.5 MPa
	sigma_15 = []
	# Shear modulus = 3.5 MPa
	sigma_35 = []

	for i in range (len(lamda)):

		a = 0.5 * (pow(lamda[i],2) - 1. /(pow(lamda[i],4)))
		b = 1.5 * (pow(lamda[i], 2) - 1. / (pow(lamda[i], 4)))
		c = 3.5 * (pow(lamda[i], 2) - 1. / (pow(lamda[i], 4)))

		sigma_05.append(a)
		sigma_15.append(b)
		sigma_35.append(c)


	import matplotlib.pyplot as plt
	print (sigma_05)
	plt.xlabel(r'$\mathrm{Stretch}$', fontsize=20)
	plt.ylabel(r'$\mathrm{\sigma_{xx}=\sigma_{yy}\/\/(MPa)}$', fontsize=20)

	plt.plot(lamda,sigma_05,  linestyle='-', linewidth=4, color='lime',label=r'$(\mu=0.5)$')
	plt.plot(lamda,sigma_15,  linestyle='-', linewidth=4, color='m',label=r'$(\mu=1.5)$')
	plt.plot(lamda,sigma_35,  linestyle='-', linewidth=4, color='orange',label=r'$(\mu=3.5)$')
	lg=plt.legend(ncol=1, loc=2, fontsize=15)
	axes = plt.gca()
	axes.set_xlim([0,5])
	axes.set_ylim([-15, 35])

	axes.set_yticks([-15,-10,-5,0,5,10,15,20,25,30,35])
	axes.set_xticks([0,1,2,3,4,5])

	plt.tick_params(axis='both', which='major', labelsize=15)
	plt.grid()
	plt.show()

And similarly the implmentation in FEniCS is presented here:

.. code-block:: python

	from dolfin import *

	# Defining Stretches
	stretch = [0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.47,5]
	BC = []

	for x in range(len(stretch)):
		 N = stretch[x] - 1.0
		 BC.append(N)

	tol = 1E-14
	# Define boundaries
	def LEFT(x, on_boundary):
		 return on_boundary and x[0] < tol

	def RIGHT(x, on_boundary):
		 return on_boundary and abs(x[0] - 1.0) < tol

	def BOTTOM(x, on_boundary):
		 return on_boundary and x[1] < tol

	def TOP(x, on_boundary):
		 return on_boundary and abs(x[1] - 1.0) < tol

	def FRONT(x, on_boundary):
		 return on_boundary and abs(x[2] - 1.0) < tol

	def BACK(x, on_boundary):
		 return on_boundary and (x[2]) < tol

	# Defining the mesh which is a single hexahedron element
	mesh = UnitCubeMesh.create(1,1,1,CellType.Type.hexahedron)


	#element for pressure field
	Element1 = FiniteElement("CG", mesh.ufl_cell(), 1)
	#element for displacement field
	Element2 = VectorElement("CG", mesh.ufl_cell(), 2)

	# Defining the mixed function space
	W_elem = MixedElement([Element1, Element2])
	W = FunctionSpace(mesh, W_elem)
	#############################################

	boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
	subdomains = MeshFunction('size_t', mesh, mesh.topology().dim())


	# Defining integration symbols
	dx = Measure('dx', domain=mesh, subdomain_data=subdomains, metadata={'quadrature_degree': 10})
	ds = Measure('ds', domain=mesh, subdomain_data=boundaries, metadata={'quadrature_degree': 10})


	## Define variational problem
	dw = TrialFunction(W)            # Incremental displacement
	v  = TestFunction(W)
	w = Function(W)
	p,u = split(w)

	#Definig some continuum mechanics relations
	d = u.geometric_dimension()
	I = Identity(d)             # Identity tensor
	F = I + grad(u)             # Deformation gradient
	C = F.T*F                   # Right cauchy stress tensor
	b = F*F.T                   # Left cauchy stress tensor
	Ib = tr(b)                  # First invariant of b tensor
	J  = det(F)                 # Jacobian of deformation gradient tensor

	# Shear modulus
	mu = 3.5e6

	# Strain energy function for fully incompressible model
	psi = mu/2.*(Ib - 3.)*dx - p*(J - 1)*dx

	# Gateaux derivative in the direction of the test function
	F1 = derivative(psi, w, TestFunction(W))

	# Compute Jacobian of F
	Jac = derivative(F1, w, TrialFunction(W))

	sigma_11 = []
	sigma_22 = []


	def border(x, on_boundary):
		 return on_boundary

	bound_x =  Expression(("t*x[0]"), degree=1, t=0)
	bound_y =  Expression(("t*x[1]"), degree=1, t=0)

	for i in range(len(BC)):

		 bound_x.t = BC[i]
		 bound_y.t = BC[i]

		 bc_x = DirichletBC(W.sub(1).sub(0), bound_x, border)
		 bc_y = DirichletBC(W.sub(1).sub(1), bound_y, border)
		 bc_front = DirichletBC(W.sub(1).sub(2), Constant((0)), FRONT)

		 bc_all = [bc_x,bc_y,bc_front]

		 problem = NonlinearVariationalProblem(F1, w, bc_all, Jac)

		 solver = NonlinearVariationalSolver(problem)

		 solver.solve()

		 (p, u) = w.split(True)

	#Stress calculation
		 sig = mu * b - p * I
	# Defining a tensor function space
		 Vv = TensorFunctionSpace(mesh, 'Lagrange', 1)
	# Projection of the stress on the tensor function space
		 sig1 = project(sig, Vv)

		 sigma_11.append((sig1.vector().get_local()[0])*0.000001)

	print (sigma_11)

	#Obtained reaults for mu = 0.5 MPa#
	sigma_05 = [-987.6430709876522, -312.4800000000002, -127.96874999999937, -61.68339506172817, -33.258200229071136, -19.451249999999895, -12.092013222069799, -7.8750000000000036, -5.31285764292055, -3.678024691358014, -2.589772373166203, -1.8374656393169522, -1.2989969135802488, -0.9007031249999968, -0.5965929377042906, -0.35707895137936096, -0.16261883157741178, 8.145274061386054e-17, 0.5764499999999999, 1.0262345679012317, 1.477938879633484, 1.9687500000000013, 2.5117407788446906, 3.1121999999999943, 3.772507427771319, 4.4938271604938125, 5.276768364202937, 6.121668054977084, 7.028721604938255, 7.99804687499999, 9.029717451299652, 10.12378067367779, 9.989197609075557, 12.499199999999991]
	
	#Obtained reaults for mu = 0.5 MPa#
	sigma_15 = [-2962.92921296296, -937.4399999999958, -383.90625000000057, -185.0501851851853, -99.77460068721344, -58.3537499999997, -36.27603966620935, -23.624999999999954, -15.938572928761642, -11.034074074074056, -7.769317119498598, -5.512396917950832, -3.8969907407407285, -2.702109374999993, -1.7897788131128693, -1.0712368541380837, -0.4878564947322365, 3.805590898063224e-16, 1.7293499999999953, 3.0787037037036926, 4.433816638900447, 5.906249999999991, 7.53522233653405, 9.336599999999983, 11.317522283314, 13.48148148148149, 15.830305092608764, 18.365004164931293, 21.08616481481487, 23.994140624999996, 27.089152353899063, 30.371342021033435, 29.967592827226596, 37.49759999999998]
	
	#Obtained reaults for mu = 0.5 MPa#
	sigma_35 = [-6913.501496913527, -2187.3599999999847, -895.7812499999939, -431.78376543209447, -232.8074016034988, -136.15874999999932, -84.6440925544885, -55.12499999999978, -37.19000350044388, -25.746172839506176, -18.128406612163456, -12.862259475218625, -9.092978395061715, -6.304921874999984, -4.176150563930013, -2.499552659655526, -1.1383318210418847, 6.882290512141593e-16, 4.035149999999988, 7.183641975308622, 10.345572157434399, 13.781249999999941, 17.582185451912803, 21.78539999999999, 26.407551994399167, 31.456790123456745, 36.93737854942055, 42.85167638483961, 49.20105123456792, 55.986328124999986, 63.2080221590976, 70.8664647157445, 69.92438326352887, 87.49440000000017]
	
	#Plotting#
	import matplotlib.pyplot as plt
	#print (sigma_05)
	plt.xlabel(r'$\mathrm{Stretch}$', fontsize=20)
	plt.ylabel(r'$\mathrm{\sigma_{xx}=\sigma_{yy}\/\/(MPa)}$', fontsize=20)

	plt.plot(stretch,sigma_05,  linestyle='-', linewidth=4, color='b',label=r'$(\mu=0.5)$')
	plt.plot(stretch,sigma_15,  linestyle='-', linewidth=4, color='g',label=r'$(\mu=1.5)$')
	plt.plot(stretch,sigma_35,  linestyle='-', linewidth=4, color='crimson',label=r'$(\mu=3.5)$')
	lg=plt.legend(ncol=1, loc=2, fontsize=15)
	axes = plt.gca()
	axes.set_xlim([0,5])
	axes.set_ylim([-15, 35])

	axes.set_yticks([-15,-10,-5,0,5,10,15,20,25,30,35])
	axes.set_xticks([0,1,2,3,4,5])

	plt.tick_params(axis='both', which='major', labelsize=15)
	plt.grid()
	plt.show()

.. figure:: PNG/3.png
   :align: center
	
   The stress results for the fully incompressible material model in biaxial loading - Finite Element (Right) vs Analytical (Left)


.. note:: 

   Considering the :math:`\textbf{Mooney-Rivlin}` hyperelastic model, the strain energy function is defined as:

   .. math:: 
     :name: eq.38 

      W= c_1(I_1-3)+c_2(I_2-3)+p(J-1)

   The stress term is defined as:  

   .. math:: 
     :name: eq.39 

      \sigma= 2c_1 \textbf{b}-2c_2 \textbf{b}^{-1} - p\textbf{I}

   The relation between the shera modulus and :math:`c_1` and :math:`c_2` constants is expressed as:  :math:`\mu=\frac{c_1+c_2}{2}`

