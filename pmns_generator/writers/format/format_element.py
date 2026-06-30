from sage.all import Integer, Rational

def format_to_int64(x:int) -> str:
    suffix = 'LL' if x < 0 else 'ULL'
    return str(x) + suffix


def format_list_to_int64(l:list) -> str:
    return '{' + ', '.join(map(format_to_int64, l)) + '}'

def format_matrix_to_int64(m, space:str =" "*4) -> str:
    lines = []
    for row in m:
        row_str = f"{space}{{{', '.join(str(x) + ('LL' if x < 0 else 'ULL') for x in row)}}}"
        lines.append(row_str)
    return '{\n' + ',\n'.join(lines) + '\n}'


def format_element_to_mpn(x, nb_limbs, phi_pow=64) -> str:
    decompose_x = []
    mask = (1 << phi_pow) - 1
    for _ in range(nb_limbs):
        part = x & mask
        decompose_x.append(f"{part}ULL")
        x >>= phi_pow
    body = ", ".join(decompose_x)
    return "{" + body + "}"

def format_extension_field_to_mpn(x, nb_limbs, space = " "*4, phi_pow=64):
    decompose_x = [space + format_element_to_mpn(int(c), nb_limbs, phi_pow) for c in x._vector_()]
    return "{\n" + ",\n".join(decompose_x) + "\n}"
    
def format_matrix_to_mpn(m, nb_limbs, space:str =" "*4, phi_pow=64) -> str:
    lines = []
    for row in m:
        formatted_elements = [format_element_to_mpn(int(x), nb_limbs, phi_pow) for x in row]
        row_str = space + "{\n" + space*2 + (",\n" + space*2).join(formatted_elements) + "\n" + space + "}"
        lines.append(row_str)

    return "{\n" + ",\n".join(lines) + "\n}"


def format_tensor_to_mpn(t, nb_limbs, space=" " * 4, phi_pow=64):
    lines = []
    for mat in t:
        formatted_mat = format_matrix_to_mpn(mat, nb_limbs, space * 2, phi_pow)
        indented_mat = "\n".join((space * 2 + line) if line else line for line in formatted_mat.splitlines())
        row_str = space + "{\n" + indented_mat + "\n" + space + "}"
        lines.append(row_str)

    return "{\n" + ",\n".join(lines) + "\n}"


def format_tensor_to_int64(t, space=" " * 4):
    lines = []
    for mat in t:
        formatted_mat = format_matrix_to_int64(mat, space)
        lines.append(formatted_mat)

    return "{\n" + ",\n".join(lines) + "\n}"

def _ensure_list(obj):
    """Convert element into a list if it has a 'rows' attribute (e.g., SageMath matrices)."""
    if hasattr(obj, 'rows'):
        return [list(row) for row in obj.rows()]
    return obj


def _pad_coefficients(poly, n: int):
    """Fill the coefficients of an object to length n with zeros."""
    if hasattr(poly, 'list'):
        poly = list(poly)
    if isinstance(poly, list) and len(poly) > 0 and isinstance(poly[0], list):
        return [_pad_coefficients(item, n) for item in poly]
    if isinstance(poly, list):
        return poly + [0] * (n - len(poly))
    return poly


def _get_dimension(val) -> int:
    """
    Recursively determine the dimension of a nested structure (list, matrix, etc.).
    Args:
        val (Any): The value whose dimension is to be determined.
    Returns:
        int: The dimension of the structure.
    """
    atomic_types = (int, float, Integer, Rational, str)
    
    # If the value is an atomic type or does not have a __getitem__ attribute, return 0 (base case).
    if isinstance(val, atomic_types) or not hasattr(val, '__getitem__'):
        return 0
        
    return 1 + _get_dimension(val[0])
