import os
import sys
from argparse import ArgumentParser

import ROOT
from ROOT import TFile, TDirectory

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--package-path')
    parser.add_argument('--source-dir')
    args = parser.parse_args()

    structure = {
        'v9': {
            2019: [
                '700',
                '712.5',
                '725',
                '737.5',
                '750',
                '762',
                '775',
                '787.5',
                '800',
                '812.5',
                '825',
                '837.5',
                '850',
                '862.5',
                '875',
                '887.5',
                '900',
                '912.5',
                '925',
                '936',
                '950',
                '962.5',
                '978',
                '975',
                '987.5',
                '955',
                '951.1',
                '945',
            ],
            2020: [
                '935',
                '945',
                '950',
                '960',
                '970',
            ],
            2021: [
                '970',
                '980',
                '990',
                '1003.5',
            ],
	    2022: [
                '937.5',
                '938.3',
                '938.9',
                '939.6',
                '940.2',
                '941',
                '942',
                '943.5',
                '945',
                '947.5',
                '952',
                '953',
                '950',
                '954',
                '951',
                '948.75',
                '935',
                '930',
                '920',
                '910',
                '900',
                '890',
                '880',
                '870',
                '860',
                '850',
                '840',
                '830',
                '820',
                '810',
                '800',
                '790',
                '790_0',
            ],
            2023: [
                '780',
                '797',
                '787.5',
                '770',
                '760',
                '750',
                '740',
                '730',
                '720',
                '710',
                '700',
	    ],
        }
    }

    ## Finding a sample hists file
    hists_files_list = list(filter(
        lambda path: path.find('hists') != -1,
        os.listdir(args.source_dir)
    ))
    if hists_files_list != []:
        sample_file = TFile.Open(os.path.normpath(args.source_dir + '/' + hists_files_list[0]), 'read')
    else:
        print("No hists files present")
        sys.exit(0)
    ## Determining a list of used cuts
    cuts_dirs = map(
        lambda obj: obj.GetName(),
        filter(
            lambda obj: isinstance(obj, TDirectory),
            map(
                lambda key: sample_file.Get( key.GetName() ),
                sample_file.GetListOfKeys()
            )
        )
    )
    cuts_dirs = dict.fromkeys(cuts_dirs)

    ## Creating a package
    package_file = TFile.Open(args.package_path, 'recreate')
    ## Determining a list of hists built for each cut
    for cut_name in cuts_dirs.keys():
        package_file.mkdir(cut_name)
        ## names of hists built for 'cut_name' cut
        cuts_dirs[cut_name] = list(map(
                lambda key: key.GetName(),
                sample_file.GetDirectory(cut_name).GetListOfKeys()
        ))

        for hist_name in cuts_dirs[cut_name]:
            package_file.GetDirectory(cut_name).mkdir(hist_name)
    ## Building a cuts results graph (TO DO)
            
    sample_file.Close()
            
    for version, years in structure.items():
        for year, energies in years.items():
            for energy in energies:
                ## Putting histograms into directories of packages
                hists_file_basename = f'hists{year}_tr_ph_fc_e{energy}_{version}.root'
                if hists_file_basename not in hists_files_list: continue
                hists_file = TFile.Open(os.path.normpath(args.source_dir + '/' + hists_file_basename), 'read')
                for cut_name, hists_names in cuts_dirs.items():
                    for hist_name in hists_names:
                        curr_dir = package_file.GetDirectory(cut_name).GetDirectory(hist_name)
                        curr_dir.cd()
                        hist = hists_file.GetDirectory(cut_name).Get(hist_name).Clone(f"{hist_name}_{year}_e{energy}_{version}")
                        hist.Write()

                ## Joining cut statistics
                

                ## Closing hists file and saving result
                hists_file.Close()
                package_file.Save()
    package_file.Close()
