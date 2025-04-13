# import fire
# from app import *

# def main():
#     opt = optManager()
#     opt.calc_1D_profile_DeltaChi2('kappa')
#     #opt.calc_1D_profile_DeltaChi2('deltaR')
#     #opt.calc_2D_profile_DeltaChi2()
#     print('Done')

# if __name__ == '__main__':
#     fire.Fire(main)



import fire
import time
from app import *

def main():

    # tags = [None, 'A', 'B']
    # opt_param_list = ['kappa_IC', 'kappa_OC', 'kappa_MT', 'ICR', 'OCR', 't13', 't23', 'mAtm', 'dcp']
    # # opt_param_list = ['ICR', 'OCR', 't13', 't23', 'mAtm', 'dcp']
    # scaling_flag = [True, False]

    tags = ['E']
    opt_param_list = ['t13', 't23', 'mAtm', 'dcp', 'kappa_IC', 'kappa_OC', 'kappa_MT', 'ICR', 'OCR']
    #opt_param_list = ['kappa_IC', 'kappa_OC', 'kappa_MT', 'ICR', 'OCR']
    
    #scaling_flag = [False, True]
    #opt_param_list = ['kappa_IC']
    scaling_flag = [False]

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


# import fire
# import time
# from app import *



# def main():

#     tags = [None, 'A', 'B']
#     # opt_param_list = ['t13', 't23', 'mAtm', 'dcp', 'kappa', 'deltaR']
#     # there seems to be some issue with deltaR...
#     #opt_param_list = ['t13', 't23', 'mAtm', 'dcp', 'kappa']
#     opt_param_list = ['mAtm']
#     scaling_flag = [False, True]

#     for tag in tags:
#         for flag in scaling_flag:
#             opt = optManager(scale_to_HK=flag, tag=tag)
#             for p in opt_param_list:
#                 stime = time.time()
#                 outname = opt.calc_1D_profile_DeltaChi2(p)
#                 print(outname)
#                 print(f"Execution time: {time.time()-stime:.0f} seconds.")

#             # df = pd.read_csv(outname)
#             # plt.plot(df.x, df.chi2, color='black', zorder=-1)
#             # plt.scatter(df.x, df.chi2, marker='o', color='red')
#             # plt.

#             # ../profile_LL_files/'

#             break

# if __name__ == '__main__':
#     fire.Fire(main)