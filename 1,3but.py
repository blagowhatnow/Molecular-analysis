import psi4
from psi4.driver.procrouting.response.scf_response import tdscf_excitations
from psi4.driver.p4util import spectrum

import psi4
from psi4.driver.procrouting.response.scf_response import tdscf_excitations
from psi4.driver.p4util import spectrum

molecule = psi4.geometry(""" 0 1
C          1.03200       -0.09230        0.07872
C          2.36955       -0.11330        0.08509
C          3.13157       -1.28640        0.43295
C          4.46929       -1.30768        0.44085
H          0.49647        0.81324       -0.18966
H          0.43892       -0.96365        0.33723
H          2.91446        0.78953       -0.18242
H          2.58577       -2.18923        0.69928
H          5.06210       -0.43574        0.18380
H          5.00436       -2.21308        0.70931
""", name= "1,3 butadiene")

psi4.set_options({
 "save_jk" : True,
})

e, wfn = psi4.energy("B3LYP/cc-pVDZ", return_wfn=True, molecule=molecule)

import numpy as np

HOMO = wfn.epsilon_a_subset("AO", "ALL").np[wfn.nalpha()]
LUMO = wfn.epsilon_a_subset("AO", "ALL").np[wfn.nalpha() + 1]

print("\nThe HOMO - LUMO gap is %16.8f hartree" % (LUMO - HOMO))
