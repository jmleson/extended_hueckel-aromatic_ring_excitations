

import multiprocessing


def worker(func, args):
    return func(*args)

def calculate_with_timeout(func, args, timeout_in_s, msg:str = "Calculation timed out."):
    with multiprocessing.Pool(processes=1) as pool:
        result = pool.apply_async(worker, (func, args))
        try:
            return result.get(timeout=timeout_in_s)
        except:
            print(msg, flush=True)
            pool.terminate()
            pool.join()  # Ensure the pool is properly cleaned up
            return None



# if __name__ == "__main__":
