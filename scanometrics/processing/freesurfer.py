"""
Wrapper for freesurfer scripts to be run on <subjects_dir>/<subjid>
"""


from shutil import rmtree, which
import subprocess
import os
from scanometrics.utils import logging

def run_recon_all(subjects_dir, subj_id, T1_file, compute_lgi=False):
    """Wrapper for freesurfer recon-all. Takes bids T1_file as input and overwrites previous recon-all output.
    Intended to use in bids freesurfer derivatives folder (bids/derivatives/freesurfer). Computation of LGI requires
    matlab to be installed.

    :param subjects_dir: path to freesurfer subject directory (eg <bids_database>/derivatives/freesurfer).
    :type subjects_dir: string
    :param subj_id: participant code of subject to run recon-all on.
    :type subj_id: string
    :param T1_file: path to T1 nifti file to be used as input for recon-all.
    :type T1_file: string
    :param compute_lgi: whether to compute the Local Gyrification Index from Schaer et al. (requires matlab, defaults to False)
    :type compute_lgi: bool
    """
    # Recover environmental variables, check that command is found
    os_env = os.environ.copy()
    os_env['SUBJECTS_DIR'] = subjects_dir
    if which('mri_convert') is None:
        logging.ERROR("""'mri_convert command not found. Please set $FREESURFER_HOME to your Freesurfer installation
        directory and run 'source $FREESURFER_HOME/SetUpFreeSurfer.sh'.""")
    # Check matlab availability
    if compute_lgi and which('matlab') is None:
        logging.ERROR("'matlab' command not found, but required to compute LGI. Please add matlab location to $PATH, or set 'compute_lgi' to False in scanometrics.processing.freesurfer.recon_all()")
    if os.path.exists(os.path.join(subjects_dir, subj_id)):
        rmtree(os.path.join(subjects_dir, subj_id))
    os.makedirs(os.path.join(subjects_dir, subj_id, 'mri', 'orig'))
    cmd = ["mri_convert", T1_file, os.path.join(subjects_dir, subj_id, 'mri', 'orig', '001.mgz')]
    p = subprocess.run(cmd, env=os_env)
    if p.returncode:
        logging.ERROR('Command failed : %s' % (' '.join(cmd)))
    os.makedirs(os.path.join(subjects_dir, subj_id, 'scripts'))
    with open(os.path.join(subjects_dir, subj_id, 'scripts', 'recon-all.log'), 'w') as fs_log_file:
        p = subprocess.run(['recon-all',
                            '-subjid', subj_id,
                            '-parcstats2',
                            '-all'], env=os_env, stdout=fs_log_file, stderr=fs_log_file)
        if p.returncode:
            logging.ERROR('recon-all failed with code %d (see %s for details).' % (p.returncode, os.path.join(subjects_dir, subj_id, 'scripts', 'recon-all.log')))
        if compute_lgi:
            # Run compute_lgi script
            p = subprocess.run(['recon-all',
                                '-subjid', subj_id,
                                '-localGI'], env=os_env, stdout=fs_log_file, stderr=fs_log_file)
            if p.returncode:
                logging.ERROR('recon-all localGI computation failed with code %d (see %s for details).' % (
                p.returncode, os.path.join(subjects_dir, subj_id, 'scripts', 'recon-all.log')))
        for hemi in ['lh','rh']:
            # WM-GM PCT
            cmd = ['mri_segstats', '--in', os.path.join(subjects_dir, subj_id, 'surf', '%s.w-g.pct.mgh' % hemi), '--annot', subj_id, hemi, 'aparc.a2009s',
                   '--sum', os.path.join(subjects_dir,subj_id,'stats','%s.w-g.pct.a2009s.stats' % hemi), '--snr']
            logging.PRINT('Running command: %s' % (' '.join(cmd)))
            p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
            if p.returncode:
                logging.ERROR('mri_segstats failed with code %d (see %s for details).' % (p.returncode, os.path.join(subjects_dir,subj_id,'scripts','recon-all.log')))
            if compute_lgi:
                # Desikan-Killiany LGI
                cmd = ['mri_segstats', '--in', os.path.join(subjects_dir, subj_id,'surf','%s.pial_lgi' % hemi), '--annot', subj_id, hemi, 'aparc', '--sum', os.path.join(subjects_dir,subj_id,'stats', '%s.pial_lgi.stats' % hemi)]
                logging.PRINT('Running command: %s' % (' '.join(cmd)))
                p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                if p.returncode:
                    logging.ERROR('mri_segstats failed with code %d (see %s for details).' % (p.returncode, os.path.join(subjects_dir,subj_id,'scripts','recon-all.log')))
                # Destrieux LGI
                cmd = ['mri_segstats', '--in', os.path.join(subjects_dir,subj_id,'surf','%s.pial_lgi' % hemi), '--annot', subj_id, hemi, 'aparc.a2009s', '--sum', os.path.join(subjects_dir,subj_id,'stats','%s.pial_lgi.a2009s.stats' % hemi)]
                logging.PRINT('Running command: %s' % (' '.join(cmd)))
                p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                if p.returncode:
                    logging.ERROR('mri_segstats failed with code %d (see %s for details).' % (p.returncode, os.path.join(subjects_dir,subj_id,'scripts','recon-all.log')))
    # Run qcache
    with open(os.path.join(subjects_dir, subj_id, 'scripts', 'recon-all_qcache.log'), 'w') as fs_log_file:
        cmd = ["recon-all", "-subject", subj_id, '-qcache']
        logging.PRINT('Running command: %s' % (' '.join(cmd)))
        p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
    # Process brain-stem structures
    with open(os.path.join(subjects_dir, subj_id, 'scripts', 'recon-all_brainstem.log'), 'w') as fs_log_file:
        cmd = ['recon-all', '-subject', subj_id, '-brainstem-structures']
        logging.PRINT('Running command: %s' % (' '.join(cmd)))
        p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
    # Process hippocampus
    with open(os.path.join(subjects_dir, subj_id, 'scripts', 'recon-all_hippocampus.log'), 'w') as fs_log_file:
        cmd = ["recon-all", '-subject', subj_id, '-hippocampal-subfields-T1']
        logging.PRINT('Running command: %s' % (' '.join(cmd)))
        p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)

def run_stats2tables(subjects_dir, subj_id,
                     seg_metrics=['volume'],
                     parc35_metrics=['volume', 'area', 'thickness', 'thicknessstd', 'meancurv', 'gauscurv', 'foldind', 'curvind', 'pctmean'],
                     parc75_metrics=['volume', 'area', 'thickness', 'thicknessstd', 'meancurv', 'gauscurv', 'foldind', 'curvind', 'pctmean'],
                     parcSTD_metrics=['mean', 'std', 'snr'],
                     computed_lgi=False):
    """Wrapper for asegstats2table and aparcstats2table. Assumes subject was processed with recon-all -parcstats2.
     Originally requires python2 but freesurfer/bin, aparcstats2table and asegstats2table were ran through 2to3."""
    if which('asegstats2table') is None:
        logging.ERROR("""asegstats2table command not found. Please set $FREESURFER_HOME to your Freesurfer installation
        directory and run 'source $FREESURFER_HOME/SetUpFreeSurfer.sh'.""")
    os_env = os.environ.copy()
    os_env['SUBJECTS_DIR'] = subjects_dir
    for metric in seg_metrics:
        cmd = ['asegstats2table',
               '--subjects', subj_id,
               '--meas', metric,
               '--all-segs',
               '--tablefile', os.path.join(subjects_dir, subj_id, 'stats2table', 'aseg_stats_%s.txt' % metric)
               ]
        logging.PRINT('Running command: %s' % (' '.join(cmd)))
        p = subprocess.run(cmd, env=os_env)
        if p.returncode:
            logging.ERROR('asegstats2table failed with code %d.' % p.returncode)
    for hemi in ['rh','lh']:
        for metric in parc35_metrics:
            if metric not in ['pctmean', 'lgimean']:
                cmd = ["aparcstats2table",
                                "--subjects", subj_id,
                                "--parc", "aparc",
                                "--hemi", hemi,
                                "--meas", metric,
                                "--parcs-from-file", os.path.join(os.path.dirname(__file__), 'resources', 'DesikanKilliany_ROIs_select.txt'),
                                "--tablefile", os.path.join(subjects_dir, subj_id, 'stats2table', '%s.aparc_stats_%s.txt' % (hemi, metric))]
                logging.PRINT('Running command: %s' % (' '.join(cmd)))
                p = subprocess.run(cmd, env=os_env)
                if p.returncode:
                    logging.ERROR("aparcstats2table failed with code %d." % p.returncode)
        for metric in parc75_metrics:
            if metric not in ['pctmean', 'lgimean']:
                cmd = ['aparcstats2table',
                    '--subjects', subj_id,
                    '--parc', 'aparc.a2009s',
                    '--hemi', hemi,
                    '--meas', metric,
                    '--parcs-from-file', os.path.join(os.path.dirname(__file__), 'resources', 'Destrieux_ROIs_select.txt'),
                    '--tablefile', os.path.join(subjects_dir, subj_id, 'stats2table', '%s.aparca2009s_stats_%s.txt' % (hemi, metric))]
                logging.PRINT('Running command: %s' % (' '.join(cmd)))
                p = subprocess.run(cmd, env=os_env)
                if p.returncode:
                    logging.ERROR("aparcstats2table failed with code %d." % p.returncode)
        # GM-WM contrast needs asegstats (table format a little different)
        for metric in parcSTD_metrics:
            cmd = ["asegstats2table",
                "--inputs", os.path.join(subjects_dir, subj_id, 'stats', '%s.w-g.pct.stats' % hemi),
                "--meas", metric,
                "--tablefile", os.path.join(subjects_dir, subj_id, 'stats2table', '%s.aparc_stats_pct%s.txt' % (hemi, metric))]
            logging.PRINT('Running command: %s' % (' '.join(cmd)))
            p = subprocess.run(cmd, env=os_env)
            cmd = ["asegstats2table",
                "--inputs", os.path.join(subjects_dir, subj_id, 'stats', '%s.w-g.pct.a2009s.stats' % hemi),
                "--meas", metric,
                "--tablefile", os.path.join(subjects_dir, subj_id, 'stats2table', '%s.aparca2009s_stats_pct%s.txt' % (hemi, metric))]
            logging.PRINT('Running command: %s' % (' '.join(cmd)))
            p = subprocess.run(cmd, env=os_env)
        # Check if lgi values where computed
        if computed_lgi:
            if not os.path.exists(os.path.join(subjects_dir,subj_id,'stats','%s.pial_lgi.stats' % hemi)):
                logging.ERROR('%s.pial_lgi.stats file not found for subject %s. Run recon-all with compute_lgi=True or disable lgi in scanometrics.processing.freesurfer.run_stats2table')
            cmd = ["asegstats2table",
                "--inputs", os.path.join(subjects_dir, subj_id, 'stats', '%s.pial_lgi.stats' % hemi),
                "--meas", "mean",
                "--tablefile", os.path.join(subjects_dir, subj_id, 'stats2table', '%s.aparc_stats_lgimean.txt' % hemi)]
            logging.PRINT('Running command: %s' % (' '.join(cmd)))
            p = subprocess.run(cmd, env=os_env)
            cmd = ["asegstats2table",
                "--inputs", os.path.join(subjects_dir, subj_id, 'stats', '%s.pial_lgi.a2009s.stats' % hemi),
                "--meas", "mean",
                "--tablefile", os.path.join(subjects_dir, subj_id, 'stats2table', '%s.aparca2009s_stats_lgimean.txt' % hemi)]
            logging.PRINT('Running command: %s' % (' '.join(cmd)))
            p = subprocess.run(cmd, env=os_env)
