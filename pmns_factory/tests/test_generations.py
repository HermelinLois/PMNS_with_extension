from sage.all import random_prime
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, TimeoutError
import sys
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[2]
if str(ROOT_DIR) not in sys.path:
    sys.path.append(str(ROOT_DIR))

import pmns_factory.pmns_E_type0_generic as type0
import pmns_factory.pmns_E_type1_generic as type1
import pmns_factory.pmns_E_type0_specific as specific_type0
import pmns_factory.pmns_E_type0_double_sparse as sparse_type0

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
    n = len(CATEGORIES)
    datas_register = {
        TYPE0: [[0]*n for _ in range(len(range_test))],
        TYPE1: [[0]*n for _ in range(len(range_test))],
        SPECIFIC_TYPE0: [[0]*n for _ in range(len(range_test))],
        SPARSE_TYPE0: [[0]*n for _ in range(len(range_test))]
    }

    for result in results:
        key = result['type'].split('.')[-1]
        data = datas_register[key]
        bit_size = result['bits']
        idx = range_test.index(bit_size)
        
        data[idx][STATUS] += 1 if result['status'] == GOOD else 0
        data[idx][NORM] += result.get('norm' ,0)
        data[idx][ERROR_TIME] += 1 if result.get('error', -1) == ERROR_FROM_TIME else 0
        data[idx][ERROR_CODE] += 1 if result.get('error', -1) == ERROR_FROM_CODE else 0
        data[idx][ERROR_UNKNOW] += 1 if result.get('error', -1) == ERROR_FROM_UNKNOW else 0
        data[idx][ROUND] += result.get('it', 0)
    
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


def run_single_generation(args: tuple) -> dict:
    """
    Run a single PMNS generation task using the prime bit size.
    """
    bit_size, k, pmns_module = args
    response = {"bits": bit_size, "type": pmns_module.__name__}
    
    try:        
        pmns = pmns_module.gen_parameters(bit_size, k)
        
        response.update({        
            "status": GOOD,  
            "norm" : sum(abs(c) for c in pmns['E'].list()), 
            "it": pmns['it'], 
        })
    except AssertionError: 
        response.update({"status": BAD, "error": ERROR_FROM_CODE})
    except Exception:
        response.update({"status": BAD, "error": ERROR_FROM_UNKNOW})
        
    return response


def run_test(k: int, ntest: int, timeout: int, range_test: list):
    # task generation
    tasks = []
    for size in range_test:
        for _ in range(ntest):
            for Etype in TYPES:
                tasks.append((size, k, Etype))
    
    results = []
    workers = multiprocessing.cpu_count()
    
    # task execution with timeout handling
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(run_single_generation, task): task for task in tasks}
        
        for future in futures:
            task = futures[future]
            bit_size, _, pmns_module = task
            try:
                res = future.result(timeout=timeout)
                results.append(res)

            except TimeoutError:
                results.append({"bits": bit_size, "type": pmns_module.__name__, "status": BAD, "error": ERROR_FROM_TIME})

            except Exception:
                results.append({"bits": bit_size, "type": pmns_module.__name__, "status": BAD, "error": ERROR_FROM_UNKNOW})
    
    with open("result_generation.txt", "a") as f:
        write_summarize_data(f, k, results, timeout, ntest, range_test)


if __name__ == "__main__":
    timeout = 60
    ntest = 10  
    range_test = [64, 128, 256, 512]
    k = 2
    
    run_test(k, ntest, timeout, range_test)