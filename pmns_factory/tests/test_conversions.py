from sage.all import *
from pathlib import Path
import sys
from time import time

ROOT_DIR = Path(__file__).resolve().parents[2]
if str(ROOT_DIR) not in sys.path:
    sys.path.append(str(ROOT_DIR))

from pmns_factory.core.operations.convertions_gestion import montgomery_exact_conversion, montgomery_pseudo_fast_conversion, montgomery_fast_conversion
from pmns_generator.writers.format.PMNS_interface import PMNSContainer
import pmns_factory.pmns_E_type0_generic as type0

PR = PolynomialRing(ZZ,"X")
X = PR.gen()

m = 128
p = random_prime(2**m, lbound=2**(m-1))
k = 2
pmns = type0.gen_parameters(p, k)


container = PMNSContainer(None, pmns)
container.add_conversions_parameters()

gamma = container.get('gamma')
K = gamma.parent()

element = K.random_element()

exact_time_start = time()
exact_poly = montgomery_exact_conversion(element, container)
exact_time_end = time()

pf_time_start = time()
pseudo_fast_poly = montgomery_pseudo_fast_conversion(element, container)
pf_time_end = time()

fast_time_start = time()
fast_poly = montgomery_fast_conversion(element, container)
fast_time_end = time()

head = f"<{'='*24} CONVERSION TIMES {'='*24}>"
print(head)
print(f"Element: {element}")
print()
print(f"Exact conversion polynomial: {exact_poly}")
print(f"Montgomery exact conversion time: {exact_time_end - exact_time_start:.6f} seconds")
print()
print(f"Pseudo-fast conversion polynomial: {pseudo_fast_poly}")
print(f"Montgomery pseudo-fast conversion time: {pf_time_end - pf_time_start:.6f} seconds")
print()
print(f"Fast conversion polynomial: {fast_poly}")
print(f"Montgomery fast conversion time: {fast_time_end - fast_time_start:.6f} seconds")
print(head)