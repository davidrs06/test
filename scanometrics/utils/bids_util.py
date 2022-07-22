"""
Several methods to automate edits to bids database.
"""


import os
from csv import DictReader
from glob import glob

def load_participants(bids_database):
    participants = {}
    with open(os.path.join(bids_database, 'participants.tsv'), 'r') as f:
        reader = DictReader(f, delimiter='\t')
        for row in reader:
            ID = row.pop('participant_id')
            participants[ID] = row
    return participants

def check_tsv_folder_congruency(bids_database):
    tsv_participants = load_participants(bids_database)
    for ID in tsv_participants.keys():
        if not os.path.exists(os.path.join(bids_database, ID)):
            print('Subject %s listed in participants.tsv but no folder found in %s' % (ID, bids_database))
    for ID in glob(os.path.join(bids_database, 'sub-*')):
        if ID not in tsv_participants.keys():
            print('Folder %s found in bids database but not in participants.tsv file.' % ID)

def add_session_files(bids_database, repeated_variables=['age']):
    for ID in load_participants().keys():
        subj_dir = os.path.join(bids_database, ID)
        sessions = glob(os.path.join(subj_dir, 'ses-*'))
        if len(sessions) > 1:
            with open(os.path.join(subj_dir, '%s_sessions.tsv' % ID), 'w') as sess_file:
                sess_file.write('session_id\t%s')
        # To be implemented

def CS_to_session(CS_database, sess_database):
    if not os.path.exists(sess_database):
        os.makedirs(sess_database)
    CS_participants = load_participants(CS_database)
    with open(os.path.join(sess_database, 'participants.tsv'), 'w') as participants_file:
        covariates = list(CS_participants[list(CS_participants.keys())[0]].keys())
        participants_file.write('participant_id\tsession_id\t%s\n' % ("\t".join(covariates)))
        ses_file = open('/tmp/test.txt', 'w')
        ses_file.close()
        for ID in CS_participants.keys():
            split_ID = ID.split('-')  # len is 1 if non-repeat or 1st reapeat, 2 for repeats
            subj_dir = os.path.join(sess_database, split_ID[0])
            if not os.path.exists(subj_dir):
                os.makedirs(subj_dir)
            if len(split_ID) == 1:  # non-repeat or 1st repeat: open session file
                ses_dir = os.path.join(subj_dir, 'ses-1')
                if not os.path.exists(ses_dir):
                    os.makedirs(ses_dir)
                if not ses_file.closed:
                    ses_file.close()
                ses_file = open(os.path.join(subj_dir, 'sessions.tsv'), 'w')
                ses_file.write('session_id\tage\n')
                ses_file.write('ses-1\t%1.2f\n' % (float(CS_participants[ID]['age'])))
                participants_file.write('%s\tses-1\t%s\n' % (split_ID[0], "\t".join(list(CS_participants[ID].values()))))
            elif len(split_ID) == 2:
                ses_file.write('ses-%d\t%1.2f\n' % (int(split_ID[1])+1,float(CS_participants[ID]['age'])))
                ses_dir = os.path.join(subj_dir, 'ses-%d' % (int(split_ID[1])+1))
                if not os.path.exists(ses_dir):
                    os.makedirs(ses_dir)
                participants_file.write('%s\tses-%d\t%s\n' % (split_ID[0],int(split_ID[1])+1,"\t".join(list(CS_participants[ID].values()))))



