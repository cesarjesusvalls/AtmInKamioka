import fire
import time
from app import *

def main():

    tags = [None, 'A', 'B']
    opt_param_list = ['t13', 't23', 'mAtm', 'dcp', 'kappa_IC', 'kappa_OC', 'kappa_MT', 'ICR', 'OCR']
    
    scaling_flag = [False, True]

    for p in opt_param_list:
        for tag in tags:
            for flag in scaling_flag:
                opt = optManager(scale_to_HK=flag, tag=tag)
                stime = time.time()
                outname = opt.calc_1D_profile_DeltaChi2(p)
                print(outname)
                print(f"Execution time: {time.time()-stime:.0f} seconds.")

if __name__ == '__main__':
    fire.Fire(main)