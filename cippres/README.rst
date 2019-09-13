=======
cippres
=======

Authors : Nico/Manu
Date : 12/09/2019

CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.

ORMAS-CI calculations are performed to compute bound and approximate continuum states. Hamiltonian and Dipole coupling matrix elements can then be computed between these states. Finally, Stieltjes imaging technique is applied to recover the correct continuum coupling matrix elements.

REQUIRMENTS:
* QP (obviously)
* python
* python libs : xml.dom, itertools, collections

USERS GUIDE:

1) set the xml input file in EZFIO
 qp set cippres finput_cippres test.xml

1') Initialize the variables n_det_max_csf, n_ciruns_cippres, n_csf_max 
Note:  Their values are read from parser.txt and set correctly in ezfio afterwards, but the code crashes if it is not initialized (ask Manu why?)
 qp set cippres n_csf_max 0
 qp set cippres n_ciruns_cippres 0
 qp set cippres n_det_max_csf 0
 qp set cippres ifcsf 0
 
 

2) generate CSFs in header.txt and list.txt  (read in EZFIO:finput_cippres)
 qp run cippres
