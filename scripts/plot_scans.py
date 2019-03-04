import chi2_scan

basedir = '/Users/James/Projects/H0LyA/'
scan_locations =   {'cf': basedir+'/data/BOSS_scans/BOSSDR11LyaF_k.scan',
                     'xcf': basedir+'/data/BOSS_scans/BOSSDR11QSOLyaF.scan'
                     }

interpolators = chi2_scan.chi2_interpolators(scan_locations,11.59,8.708)
interpolators.plot_interpolated_scan(corr_type='cf',save=True)
interpolators.plot_interpolated_scan(corr_type='xcf',save=True)
