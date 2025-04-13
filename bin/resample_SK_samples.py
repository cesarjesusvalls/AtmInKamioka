import fire
from app import *

def main():
    anaMaster = AnaMaster(dont_replace_names=True, mask_bins=False)
    anaMaster.resample_from_quantiles()
    anaMaster.fill_histograms()
    
    print('---->  anaMaster is Ready.\n')

    print('Creating Fake Data "A".\n')
    #### FDS A ####
    process_FDS (anaMaster, "A", create_bins_using_narrow_binning(anaMaster))
    print('Done!\n')

    print('Creating Fake Data "B".\n')
    #### FDS B ####
    process_FDS (anaMaster, "B", create_bins_using_narrow_binning(anaMaster))
    print('Done!\n')

    print('Creating Fake Data "C".\n')
    #### FDS C ####
    process_FDS (anaMaster, "C", create_bins_using_narrow_binning(anaMaster))
    print('Done!\n')

    print('Creating Fake Data "D".\n')
    #### FDS D ####
    process_FDS (anaMaster, "D", create_bins_using_narrow_binning(anaMaster))
    print('Done!\n')

    print('Creating Fake Data "E".\n')
    #### FDS E ####
    process_FDS (anaMaster, "E", create_bins_using_narrow_binning(anaMaster))
    print('Done!\n')

    print('Creating Fake Data "F".\n')
    #### FDS F ####
    process_FDS (anaMaster, "F", create_bins_using_narrow_binning(anaMaster))
    print('Done!\n')

    print('Creating Fake Data "G".\n')
    #### FDS G ####
    process_FDS (anaMaster, "G", create_bins_using_narrow_binning(anaMaster))
    print('Done!\n')

    print('Creating Fake Data "H".\n')
    #### FDS H ####
    process_FDS (anaMaster, "H", create_bins_using_narrow_binning(anaMaster))
    print('Done!\n')

    print('Creating Fake Data "I".\n')
    #### FDS I ####
    process_FDS (anaMaster, "I", create_bins_using_narrow_binning(anaMaster))
    print('Done!\n')

    print('Creating Fake Data "J".\n')
    #### FDS J ####
    process_FDS (anaMaster, "J", create_bins_using_narrow_binning(anaMaster))
    print('Done!\n')

    print('---->  All Fake Data is Ready.\n')

    # Let's store plots for all the new data:
    tags = ["J"]
    for tag in tags:
        ana_new = AnaMaster("../fake_data/FDS_"+tag+"/unoscillated/",  binning_file="../fake_data/FDS_"+tag+"/binning.txt", mask_bins=False)
        ana_new.fill_histograms()
        for i in range(len(ana_new.samples)):
            s = ana_new.samples[i]
            c = ana_new.samples[i].plot()
            c.Draw()
            c.SaveAs("../fake_data/FDS_"+tag+"/figures/"+s.name+'.png')

if __name__ == '__main__':
    fire.Fire(main)