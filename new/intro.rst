1. Introduction
===========================================


The strain energy function :math:`W` is defined as:

.. math:: 
  :name: eq.1

   W=W(I_1 (\textbf {C}),I_2 (\textbf {C}),I_3 (\textbf {C}))

In the above equations, :math:`\textbf {C}` is the right Cauchy deformation tensor:

.. math:: 
  :name: eq.2

  \textbf {C}=F^TF


Where :math:`F` is the deformation gradient tensor: :math:`F=I+\nabla u`. The Jacobian of :math:`F` is defined as its determinant:  :math:`J= \lvert F \rvert`.

The :math:`I_1 (\textbf {C})`, :math:`I_2 (\textbf {C})` and :math:`I_3 (\textbf {C})` are the invariants of the :math:`\textbf {C}` tensor. 

.. math:: 
  :name: eq.3

  I_1=tr(\textbf{C}) 

  I_2=0.5[(tr \textbf{C})^2-tr(\textbf{C}^2)] 

  I_3=det (\textbf{C})


The second Piola-Kirchhoff stress is defined as: 

.. math:: 
  :name: eq.4

  S=2 \frac{\partial W}{\partial \textbf {C}} 

After applying the chain rule, we can write the term :math:`\frac{\partial W}{\partial \textbf {C}}` as: 


.. math:: 
  :name: eq.5

  \frac{\partial W}{\partial \textbf {C}}=\frac{\partial W}{I_1} \frac{\partial I_1}{\textbf {C}} + \frac{\partial W}{I_2} \frac{\partial I_2}{\textbf {C}}+\frac{\partial W}{I_3} \frac{\partial I_3}{\textbf {C}}

The first term on the right hand side, we could apply the chain rul to calculate the :math:`\frac{\partial I_1}{\partial C}`:

.. math:: 
  :name: eq.6

  \frac{\partial I_1}{\partial C}=\frac{\partial (I_1)_{11}}{C_{ij}} e_i \otimes e_j + \frac{\partial (I_1)_{22}}{C_{ij}} e_i \otimes e_j + \frac{\partial (I_1)_{33}}{C_{ij}} e_i \otimes e_j 

  \frac{\partial I_1}{\partial C}=e_1 \otimes e_1 +e_2 \otimes e_2 +e_3 \otimes e_3 = \textbf{I}

Similary we can can derive the :math:`\frac{\partial I_2}{\partial C}` and :math:`\frac{\partial I_3}{\partial C}`:


.. math:: 
  :name: eq.7

  \frac{\partial I_2}{\partial \textbf {C}}= I_1 \textbf {I} - \textbf {C}

  \frac{\partial I_3}{\partial \textbf {C}}= I_2 \textbf {I} - I_1 \textbf {C} + \textbf {C}^2

By plugging back the  :ref:`equation.6 <eq.6>` and :ref:`equation.7 <eq.7>` into the :ref:`equation.5 <eq.5>`, we can find the second piola kirchhoff stress as: 

.. math:: 
  :name: eq.8

  \textbf {S}= 2 \alpha_1 \textbf {I} + 2 \alpha_2 \textbf {C} + 2 \alpha_3 \textbf {C}^2

In :ref:`Equation.8 <eq.8>`: 


.. math:: 
  :name: eq.9

  \alpha_1 =  \frac{\partial W}{\partial I_1} +  \frac{\partial W}{\partial I_2} I_1 +  \frac{\partial W}{\partial I_3} I_2

  \alpha_2 =  -\frac{\partial W}{\partial I_2} -  \frac{\partial W}{\partial I_3} I_1 

  \alpha_3 =  \frac{\partial W}{\partial I_3} 

The Cauchy stress is defined as: 

.. math:: 
  :name: eq.10

  \sigma = J^{-1} F F^T \textbf {S}

Similar to the :ref:`Equation.1 <eq.1>` we can define: 


.. math:: 
  :name: eq.11


   W=W(I_1 (\textbf {b}),I_2 (\textbf {b}),I_3 (\textbf {b}))



Where :math:`\textbf {b}` is the left Cauchy deformation tensor:


.. math:: 
  :name: eq.12

  \textbf {b}=FF^T

It should be noted that the invariants of the \textbf {b} are obtained similar to the :ref:`equation.3 <eq.3>`. By taking the same steps shown in the :ref:`equation.6 <eq.6>` and :ref:`equation.7 <eq.7>`

In addition, we take advantage of Cayley-Hamilton theorem:

.. math:: 
  :name: eq.13

  \textbf{b}^3-I_1 \textbf{b}^2+I_2 \textbf{b} - I_3 \textbf{I} =0


So the term :math:`\frac{\partial I_3}{\partial \textbf {b}}` could be rewritten again: 

.. math:: 
  :name: eq.14

  \frac{\partial I_3}{\partial \textbf {b}}= I_3 \textbf {b}^{-1}

We can write the :math:`\frac{\partial W}{\partial \textbf {b}}`: 

.. math:: 
  :name: eq.15

  \frac{\partial W}{\partial \textbf {b}}= (\frac{\partial W}{\partial I_1} + \frac{\partial W}{\partial I_2} I_1) \textbf {I} - \frac{\partial W}{\partial I_2}\textbf {b} +\frac{\partial W}{\partial I_3} I_3 \textbf {b}^{-1}

After multiplying the :math:`\textbf {b}` tensor into the :ref:`equation.15 <eq.15>` and then by substitution in :ref:`equation.10 <eq.10>`, the stress could be defined in the form form:

    
.. math:: 
  :name: eq.16

   \sigma=\beta_1 \textbf {I}+\beta_2 \textbf {b}+\beta_3 \textbf {b}^2

In the above equation: 


.. math:: 
  :name: eq.17

  \beta_1 = 2J^{-1} (\frac{\partial W}{\partial I_3}I_3) 

  \beta_2 = 2J^{-1} (\frac{\partial W}{\partial I_1} +  \frac{\partial W}{\partial I_2} I_1) 

  \beta_3 = 2J^{-1} (-\frac{\partial W}{\partial I_2}) 


