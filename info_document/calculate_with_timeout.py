


import concurrent.futures
import sympy as sp
import multiprocessing
import time


def worker(func, args):
    return func(*args)

def calculate_with_timeout(func, args, timeout_in_s, msg:str = "Calculation timed out."):
    with multiprocessing.Pool(processes=1) as pool:
        result = pool.apply_async(worker, (func, args))
        try:
            return result.get(timeout=timeout_in_s)
        except:# multiprocessing.TimeoutError:
            print(msg)
            pool.terminate()
            pool.join()  # Ensure the pool is properly cleaned up
            return None
        # except AttributeError:# because result is None and thereby has not attribute "get"
        #     print(msg)
        #     pool.terminate()
        #     pool.join()  # Ensure the pool is properly cleaned up
        #     return None



if __name__ == "__main__":

    def test_function():
        while True:
            pass

    calculate_with_timeout(test_function, (), 5)