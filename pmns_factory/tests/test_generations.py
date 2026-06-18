from sage.all import random_prime
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
import signal
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[2]
if str(ROOT_DIR) not in sys.path:
    sys.path.append(str(ROOT_DIR))

import pmns_factory.pmns_E_type0_generic as type0
import pmns_factory.pmns_E_type1_generic as type1
import pmns_factory.pmns_E_type0_specific as specific_type0
import pmns_factory.pmns_E_type0_sparse as sparse_type0


GOOD = 'GOOD'
BAD = 'BAD'
ERROR_FROM_TIME = "ERROR FROM TIME"
ERROR_FROM_CODE = "ERROR FROM CODE"
ERROR_FROM_UNKNOW = "ERROR FROM UNKNOW"

WRITE_SPACE = 30
ERROR_SPACE = 35

STATUS = 0
NORM = 1
ERROR_TIME = 2
ERROR_CODE = 3
ROUND = 4
ERROR_UNKNOW = 5
CATEGORIES = [STATUS, NORM, ERROR_TIME, ERROR_CODE, ROUND, ERROR_UNKNOW]

TYPE0 = type0.__name__.split('.')[-1]
TYPE1 = type1.__name__.split('.')[-1]
SPECIFIC_TYPE0 = specific_type0.__name__.split('.')[-1]
SPARSE_TYPE0 = sparse_type0.__name__.split('.')[-1]
TYPES = [type0, type1, specific_type0, sparse_type0]


def write_summarize_data(writer, k:int, results:list, timeout:int, ntest:int, range_test:list):
    """
    compute average data of pmns generation for each prime size and write them in a text file
    
    Args:
        writer: object used to write in a file
        k(int): extension degree
        results(list): list of pmns parameters
        timeout(int): seconds for maximum limit computation of a pmns
        ntest(int): number of test for each prime size
        range_test(list): list of bit size used for testing
    Return:
        void: write summarized data in result_generation text file
    """
    
    n = len(CATEGORIES)
    datas_register = {
        TYPE0: [[0]*n for _ in range(len(range_test))],
        TYPE1: [[0]*n for _ in range(len(range_test))],
        SPECIFIC_TYPE0: [[0]*n for _ in range(len(range_test))],
        SPARSE_TYPE0: [[0]*n for _ in range(len(range_test))]
    }

    # write specific data in registers
    for result in results:
        key = result['type'].split('.')[-1]
        data = datas_register[key]
        bit_size = result['p'].nbits()
        idx = range_test.index(bit_size)
        
        data[idx][STATUS] += 1 if result['status'] == GOOD else 0
        data[idx][NORM] += result.get('norm' ,0)
        data[idx][ERROR_TIME] += 1 if result.get('error', -1) == ERROR_FROM_TIME else 0
        data[idx][ERROR_CODE] += 1 if result.get('error', -1) == ERROR_FROM_CODE else 0
        data[idx][ERROR_UNKNOW] += 1 if result.get('error', -1) == ERROR_FROM_UNKNOW else 0
        data[idx][ROUND] += result.get('it', 0)
    
    #write file
    writer.write(f"# WITH TIMEOUT = {timeout} AND K = {k}\n")

    header = (
        f"|{'BIT_SIZE':^{WRITE_SPACE}}|"
        f"{'TYPE':^{WRITE_SPACE}}|"
        f"{'ROUND':^{WRITE_SPACE}}|"
        f"{'NORM':^{WRITE_SPACE}}|"
        f"{GOOD:^{WRITE_SPACE}}|"
        f"{'ERRORS (CODE|TIME|UNKNOW)':^{ERROR_SPACE}}|\n"
    )
    
    break_line = "-"*len(header) + "\n"
    edge_line = '='*len(header) + "\n"

    writer.write(edge_line)
    writer.write(header)
    writer.write(edge_line)
    
    # fromat to txt parameters
    avg = lambda data, idx: float(data[idx]) / data[STATUS] if data[STATUS] else .0
    gen_text = lambda name, size, data: (
        f"|{size:^{WRITE_SPACE}}|"
        f"{name:^{WRITE_SPACE}}|"
        f"{avg(data, ROUND):^{WRITE_SPACE}.2f}|"
        f"{avg(data, NORM):^{WRITE_SPACE}.2f}|"
        f"{f'{data[STATUS]}/{ntest}':^{WRITE_SPACE}}|"
        f"{f'{data[ERROR_CODE]}/{ntest} | {data[ERROR_TIME]}/{ntest} | {data[ERROR_UNKNOW]}/{ntest}':^{ERROR_SPACE}}|\n"
    )
    
    for idx in range(len(range_test)):
        bit_size = range_test[idx]
        
        for Etype in datas_register.keys():
            result = datas_register[Etype][idx]
            writer.write(gen_text(Etype, bit_size, result))
    
        writer.write(break_line)
    

def alarm_handler(null1, null2):
    """
    Alarm function to raise time error during test generation
    """
    raise TimeoutError("generation timeout")


def run_single_generation(args:list)-> dict:
    """
    Run a test generation
    
    Args:
        p(int): prime used
        k(int): extension field used
        pmns_module(class): class type used to generate pmns parameters
        timeout(int): maximum time execution for generation
        
    Return:
        dict: register pmns data generation
    """
    p, k, pmns_module, timeout = args
    response = {"p": p, "type": pmns_module.__name__}
    
    signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(timeout)
    # try to generate parameters
    # generation is canceled if execution time is exceeded
    try:
        pmns = pmns_module.gen_parameters(p, k)
        response.update({        
                "status": GOOD,  
                "norm" : sum(abs(c) for c in pmns['E'].list()), 
                "it": pmns['it'], 
                })
    except TimeoutError:
        response.update({"status": BAD, "error": ERROR_FROM_TIME})

    # catch assertion error coming from code
    except AssertionError: 
        response.update({"status": BAD, "error": ERROR_FROM_CODE})

    # catch general exception
    except Exception:
        response.update({"status": BAD, "error": ERROR_FROM_UNKNOW})
        
    finally:
        signal.alarm(0)
    
    return response


def gen_prime(size:int):
    """
    gen random prime of a chosen bit size
    """
    return random_prime(2**size, lbound=2**(size-1))


def run_test(k:int, ntest:int, timeout:int, range_test:list):
    """
    execute pmns generation for given sizes and extension degree
    
    Args:
        k(int): extension degree used for pmns generation
        ntest(int): number of test run for every size of prime given
        timeout(int): maximum seconds used for pmns parameters generation
    Return:
        void: write data about pmns generation
    """
    # primes generation 
    primes_slots = []
    for size in range_test:
        primes_slots.extend([size] * ntest)

    primes = []
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        primes = list(executor.map(gen_prime, primes_slots))    
    
    # generate parameters for pmns geenration
    tasks = []
    for p in primes:
        for Etype in TYPES:
            tasks.append((p, k, Etype, timeout))
    
    # generate result
    results = []
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        results = list(executor.map(run_single_generation, tasks))
    
    # write data in text file
    with open(f"result_generation.txt", "a") as f:
        write_summarize_data(f, k, results, timeout, ntest, range_test)
    
    
if __name__ == "__main__":
    timeout = 60
    ntest = 1000
    range_test = [64, 128, 256, 512]
    k = 2
    
    run_test(k, ntest, timeout, range_test)