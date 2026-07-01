from sage.all import Integer, Rational

def format_to_int64(x:int) -> str:
    """
    Format an integer to a string representation suitable for C/C++ code, 
    appending 'LL' for negative values and 'ULL' for non-negative values.
    Args:
        x (int): The integer to format.
    Returns:
        str: The formatted string representation of the integer.
    """
    suffix = 'LL' if x < 0 else 'ULL'
    return str(x) + suffix


def format_list_to_int64(l:list) -> str:
    """
    Format a list of intergers to a string representation suitable for C/C++ code,
    using the format_to_int64 function for each element.
    Args:
        l (list): The list of integers to format.
    Returns:
        str: The formatted string representation of the list.
    """
    return '{' + ', '.join(map(format_to_int64, l)) + '}'


def format_matrix_to_int64(m, space:str =" "*4) -> str:
    """
    Format a matrix of intergers to a string representation suitable for C/C++ code,
    appending 'LL' for negative values and 'ULL' for non-negative values..
    Args:
        m (list): The matrix of integers to format.
        space (str): The indentation space for formatting.
    Returns:
        str: The formatted string representation of the matrix.
    """
    lines = []
    for row in m:
        row_str = f"{space}{{{', '.join(str(x) + ('LL' if x < 0 else 'ULL') for x in row)}}}"
        lines.append(row_str)
    return '{\n' + ',\n'.join(lines) + '\n}'


def format_element_to_mpn(x, nb_limbs:int, phi_pow=64) -> str:
    """
    Format an integer to a string representation suitable for gmp library, more specifically for mpn functions, 
    by decomposing it into a list of machine words (of size phi_pow) using nb_limbs registers.
    Args:
        x (int): The integer to format.
        nb_limbs (int): The number of machine words to use for the representation.
        phi_pow (int): The size of each machine word in bits (default is 64).
    Returns:
        str: The formatted string representation of the integer.
    """
    decompose_x = []
    mask = (1 << phi_pow) - 1
    for _ in range(nb_limbs):
        part = x & mask
        decompose_x.append(f"{part}ULL")
        x >>= phi_pow
    body = ", ".join(decompose_x)
    return "{" + body + "}"


def format_extension_field_to_mpn(x, nb_limbs:int, space = " "*4, phi_pow=64)-> str:
    """
    Format an extension field element to a string representation suitable for mpn representation, 
    by applying the format_element_to_mpn function to each coefficient of the extension field element.
    Args:
        x (int): The integer to format.
        nb_limbs (int): The number of machine words to use for the representation.
        phi_pow (int): The size of each machine word in bits (default is 64).
        space (str): The indentation space for formatting.
    Returns:
        str: The formatted string representation of the integer.
    """
    decompose_x = [space + format_element_to_mpn(int(c), nb_limbs, phi_pow) for c in x._vector_()]
    return "{\n" + ",\n".join(decompose_x) + "\n}"
    

def format_matrix_to_mpn(m, nb_limbs:int, space:str =" "*4, phi_pow=64) -> str:
    """
    Format a matrix of intergers to a string representation suitable for mpn representation,
    by applying the format_element_to_mpn function to each element of the matrix.
    Args:
        m (list): The matrix of integers to format.
        nb_limbs (int): The number of machine words to use for the representation.
        phi_pow (int): The size of each machine word in bits (default is 64).
        space (str): The indentation space for formatting.
    Returns:
        str: The formatted string representation of the matrix.
    """
    lines = []
    for row in m:
        formatted_elements = [format_element_to_mpn(int(x), nb_limbs, phi_pow) for x in row]
        row_str = space + "{\n" + space*2 + (",\n" + space*2).join(formatted_elements) + "\n" + space + "}"
        lines.append(row_str)

    return "{\n" + ",\n".join(lines) + "\n}"


def format_tensor_to_mpn(t, nb_limbs:int, space=" " * 4, phi_pow=64) -> str:
    """
    Format a tensor (3D array) of integers to a string representation suitable for mpn representation,
    by applying the format_matrix_to_mpn function to each matrix in the tensor.
    Args:
        t (list): The tensor of integers to format.
        nb_limbs (int): The number of machine words to use for the representation.
        phi_pow (int): The size of each machine word in bits (default is 64).
        space (str): The indentation space for formatting.
    Returns:
        str: The formatted string representation of the tensor.
    """
    lines = []
    for mat in t:
        formatted_mat = format_matrix_to_mpn(mat, nb_limbs, space * 2, phi_pow)
        indented_mat = "\n".join((space * 2 + line) if line else line for line in formatted_mat.splitlines())
        row_str = space + "{\n" + indented_mat + "\n" + space + "}"
        lines.append(row_str)

    return "{\n" + ",\n".join(lines) + "\n}"


def format_tensor_to_int64(t, space=" " * 4):
    """
    Format a tensor (3D array) of integers of less than 64 bits to a string representation suitable for C/C++ code, 
    by applying the format_matrix_to_int64 function to each matrix in the tensor.
    Args:
        t (list): The tensor of integers to format.
        space (str): The indentation space for formatting.
    Returns:
        str: The formatted string representation of the tensor.
    """
    lines = []
    for mat in t:
        formatted_mat = format_matrix_to_int64(mat, space)
        lines.append(formatted_mat)

    return "{\n" + ",\n".join(lines) + "\n}"

def _ensure_list(obj):
    """
    Convert element into a list if it has a 'rows' attribute.
    Args:
        obj (Any): The object to check and convert.
    Returns:
        list: The object as a list if it has 'rows', otherwise the original object.
    """
    if hasattr(obj, 'rows'):
        return [list(row) for row in obj.rows()]
    return obj


def _pad_coefficients(obj, n: int) -> list:
    """
    Fill the coefficients of an object to length n with zeros.
    Args:
        obj (Any): The object whose coefficients need to be padded. 
            The object must have have a method to have his length determined.
        n (int): The desired length of the coefficient list.
    Returns:
        list: The padded coefficient list.
    """
    if hasattr(obj, 'list'):
        obj = list(obj)

    if isinstance(obj, list) and isinstance(obj[0], list):
        return [_pad_coefficients(item, n) for item in obj]
    
    if isinstance(obj, list):
        return obj + [0] * (n - len(obj))
    return obj


def _get_dimension(obj) -> int:
    """
    Recursively determine the dimension of a nested structure (list, matrix, etc.).
    Args:
        obj (Any): The object whose dimension is to be determined.
    Returns:
        int: The dimension of the structure.
    """
    atomic_types = (int, float, Integer, Rational)
    
    # If the value is an atomic type or does not have a __getitem__ attribute, return 0 (base case).
    if isinstance(obj, atomic_types) or not hasattr(obj, '__getitem__'):
        return 0
        
    return 1 + _get_dimension(obj[0])
