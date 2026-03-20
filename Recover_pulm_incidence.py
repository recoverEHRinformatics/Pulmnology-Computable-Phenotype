# pulmonary incidence
import pandas as pd
import numpy as np
import pickle as pkl
import time
import re

BASELINE_LEFT = -1095 # day
BASELINE_RIGHT = -7 # day

FOLLOWUP_LEFT = 31 # day
FOLLOWUP_RIGHT = 180 # day

OP_apart_threshold = 30 # day
flag_exclusion_duration = '-7_-3_years' # '-7-3years','0-3years'


def _is_in_baseline(event_time, index_time):
    return BASELINE_LEFT <= (event_time - index_time).days <= BASELINE_RIGHT

def _is_in_followup(event_time, index_time):
    return FOLLOWUP_LEFT <= (event_time - index_time).days <= FOLLOWUP_RIGHT

def ndc_normalization(x):
    # https://www.nlm.nih.gov/research/umls/rxnorm/docs/techdoc.html#s6_0
    # "You can download the algorithm used for normalizing NDC codes (RTF document)."
    # https://www.nlm.nih.gov/research/umls/rxnorm/docs/techdoc.html#sat
    if len(x) - len(x.replace('-', '')) == 2:  # If NDC string contains 2 dashes
        a, b, c = x.split('-')
        na, nb, nc = [len(t) for t in x.split('-')]  # a-b-c format
        if (na == 6) and (nb == 4) and (nc == 2):
            ndc = a[1:] + b + c
        elif (na == 6) and (nb == 4) and (nc == 1):
            ndc = a[1:] + b + '0' + c
        elif (na == 6) and (nb == 3) and (nc == 2):
            ndc = a[1:] + '0' + b + c
        elif (na == 6) and (nb == 3) and (nc == 1):
            ndc = a[1:] + '0' + b + '0' + c
        elif (na == 5) and (nb == 4) and (nc == 2):  # MMX, GS, CVX
            ndc = a + b + c
        elif (na == 5) and (nb == 4) and (nc == 1):
            ndc = a + b + '0' + c
        elif (na == 5) and (nb == 3) and (nc == 2): # MTHSPL
            ndc = a + '0' + b + c
        elif (na == 4) and (nb == 4) and (nc == 2):
            ndc = '0' + a + b + c
        else:
            ndc = ''
    elif (len(x) == 11) and ('-' not in x):  # rxnorm, NDDF, MMSL
        ndc = x
    elif (len(x) == 12) and ('-' not in x):  # and (x[0] == '0'): #VANDF, some initial with 9
        ndc = x[1:]
    else:
        ndc = ''

    ndc = ndc.replace('*', '0')

    if re.sub(r"[0123456789]", '', ndc):
        ndc = ''

    return ndc

def read_code_list(pulmonary_CL_file):

    pulm_CL = pulmonary_CL_file + 'RECOVERY Pulmonary CP Code Lists_8.22.22-NDC-crosswalk.xlsx' #

    # dx code
    pulm_dx_1a = pd.read_excel(pulm_CL, sheet_name='1a. Specific DX - DX1', dtype={'DX Code': str})
    dx_1a = list(pulm_dx_1a['DX Code'])
    dx_1a = [i.strip(' ') for i in dx_1a]

    # med
    pulm_med = pd.read_excel(pulm_CL, sheet_name='Medications', dtype= str)

    med_rx1_cui = list(set(pulm_med.loc[pulm_med['Code List'] == 'RX1', 'rxcui'])) # using rxcui ? RXNORM_CODE ?  len(med_rx1_cui) 885   len(RXNORM_CODE) 352, '828284','828248'
    med_rx2_cui = list(set(pulm_med.loc[pulm_med['Code List'] == 'RX2', 'rxcui']))
    med_rx1_norm = list(set(pulm_med.loc[pulm_med['Code List'] == 'RX1', 'RXNORM_CODE']))
    med_rx2_norm = list(set(pulm_med.loc[pulm_med['Code List'] == 'RX2', 'RXNORM_CODE']))

    med_rx1_ndc = list(set(pulm_med.loc[pulm_med['Code List'] == 'RX1', 'ndc_code']))
    med_rx2_ndc = list(set(pulm_med.loc[pulm_med['Code List'] == 'RX2', 'ndc_code']))

    # remove nan
    med_rx1_cui = [item for item in med_rx1_cui if str(item) != 'nan']
    med_rx2_cui = [item for item in med_rx2_cui if str(item) != 'nan']
    med_rx1_norm = [item for item in med_rx1_norm if str(item) != 'nan']
    med_rx2_norm = [item for item in med_rx2_norm if str(item) != 'nan']

    med_rx1_ndc = [item for item in med_rx1_ndc if str(item) != 'nan']
    med_rx2_ndc = [item for item in med_rx2_ndc if str(item) != 'nan']

    # clean
    med_rx1_cui = [i.strip(' ') for i in med_rx1_cui]
    med_rx2_cui = [i.strip(' ') for i in med_rx2_cui]
    med_rx1_norm = [i.strip(' ') for i in med_rx1_norm]
    med_rx2_norm = [i.strip(' ') for i in med_rx2_norm]

    med_rx1_ndc = [i.strip(' ') for i in med_rx1_ndc]
    med_rx2_ndc = [i.strip(' ') for i in med_rx2_ndc]

    # normalize ndc code 73521-0030-01 -> 73521003001
    med_rx1_ndc_normalized = [ndc_normalization(i) for i in med_rx1_ndc]
    med_rx2_ndc_normalized = [ndc_normalization(i) for i in med_rx2_ndc]

    med_rx1 = med_rx1_cui + med_rx1_norm + med_rx1_ndc_normalized
    med_rx2 = med_rx2_cui + med_rx2_norm + med_rx2_ndc_normalized

    return dx_1a, med_rx1, med_rx2

def _dx_clean_and_translate_any_ICD9_to_ICD10(dx_list, icd9_icd10, icd_ccsr):
    dx_list_new = []
    n_icd9 = 0
    for records in dx_list:
        dx_t, icd, dx_type, enc_type = records
        icd = icd.replace('.', '').upper().strip()
        if icd == '':
            # e.g. in mshs, many dx colum is '', not icd10 encoded
            continue
        if ('9' in dx_type) or ((icd.isnumeric() or (icd[0] in ('E', 'V'))) and (icd not in icd_ccsr)):
            icd10_translation = icd9_icd10.get(icd, [])
            if len(icd10_translation) > 0:
                n_icd9 += 1
                for x in icd10_translation:
                    new_records = (dx_t, x, 'from' + icd, enc_type)
                    dx_list_new.append(new_records)
                # print(icd, icd10_translation)
            else:
                dx_list_new.append((dx_t, icd, dx_type, enc_type))
        else:
            dx_list_new.append((dx_t, icd, dx_type, enc_type))

    return dx_list_new

def read_dx_mapping_function(code_path):
    icd9_icd10_f = open(code_path + "icd9_icd10.pkl", "rb")
    icd9_icd10 = pkl.load(icd9_icd10_f)
    icd_ccsr_f = open(code_path + "icd_ccsr_mapping.pkl", "rb")
    icd_ccsr = pkl.load(icd_ccsr_f)
    return icd9_icd10, icd_ccsr

def func_read_ori_data_server(candidate_folder_path, candidate_site, icd9_icd10, icd_ccsr):
    # read procedure, diagnosis files;
    covid_data = pd.read_csv(candidate_folder_path +candidate_site + '\\'+'matrix_cohorts_covid_4manuNegNoCovidV2age18_boolbase-nout-withAllDays-withPreg_' + candidate_site + '.csv',  dtype={'patid': str, 'zip': str}, parse_dates=['index date','dob'], index_col=False)
    covid_data['index date'] = pd.to_datetime(covid_data['index date']).dt.date

    id_diagnosis_f = open(candidate_folder_path +candidate_site + '\\' + "diagnosis_" + candidate_site + ".pkl", "rb")
    id_diagnosis = pkl.load(id_diagnosis_f) # id_diagnosis is dictionary, and patid is key; value is a list with each item: such assample: (Timestamp('2020-08-14 00:00:00'), 'R07.2', '10', 'AV'),

    id_medication_f = open(candidate_folder_path +candidate_site + '\\' + "medication_" + candidate_site + ".pkl", "rb")
    id_medication = pkl.load(id_medication_f)  # id_medication is dictionary, and patid is key; value is a list with each item: such as (Timestamp('2006-12-22 00:00:00'), '757969', 364)

    print('converting ICD9 to ICD10...')
    # for id_diagnosis, we convert ICD9 to ICD10 code
    for key, value in id_diagnosis.items():
        patid = key
        dx_list_old = value
        # print('dx_list_old:',dx_list_old)
        dx_list = _dx_clean_and_translate_any_ICD9_to_ICD10(dx_list_old, icd9_icd10, icd_ccsr)
        # print('dx_list_new:', dx_list)
        id_diagnosis[patid] = dx_list

    return covid_data, id_diagnosis, id_medication

def func_read_ori_data(candidate_folder_path, candidate_site, icd9_icd10, icd_ccsr):
    # read procedure, diagnosis files;
    covid_data = pd.read_csv(candidate_folder_path +candidate_site + '/'+'matrix_cohorts_covid_simplified_' + candidate_site + '.csv', dtype={'patid': str, 'zip': str}, parse_dates=['index date'], index_col=False)
    covid_data['index date'] = pd.to_datetime(covid_data['index date']).dt.date

    id_diagnosis_f = open(candidate_folder_path +candidate_site + '/' + "diagnosis_simplified_" + candidate_site + ".pkl", "rb")
    id_diagnosis = pkl.load(id_diagnosis_f) # id_diagnosis is dictionary, and patid is key; value is a list with each item: such assample: (Timestamp('2020-08-14 00:00:00'), 'R07.2', '10', 'AV'),

    id_medication_f = open(candidate_folder_path +candidate_site + '/' + "medication_simplified_" + candidate_site + ".pkl", "rb")
    id_medication = pkl.load(id_medication_f)  # id_medication is dictionary, and patid is key; value is a list with each item: such as (Timestamp('2006-12-22 00:00:00'), '757969', 364)

    print('converting ICD9 to ICD10...')
    # for id_diagnosis, we convert ICD9 to ICD10 code
    for key, value in id_diagnosis.items():
        patid = key
        dx_list_old = value
        # print('dx_list_old:',dx_list_old)
        dx_list = _dx_clean_and_translate_any_ICD9_to_ICD10(dx_list_old, icd9_icd10, icd_ccsr)
        # print('dx_list_new:', dx_list)
        id_diagnosis[patid] = dx_list

    return covid_data, id_diagnosis, id_medication

def func_set_pulm_flag(candidate_df):
    add_flag_baseline = ['flag_baseline_pulmonary',
                         'flag_baseline_dx', 'flag_baseline_med']  # flag is used to exclude pulmonary at baseline

    add_flag_followup = ['flag_followup_pulmonary',
                         'flag_pulmonary_dx', 'flag_pulmonary_dx_t2e',
                         'flag_pulmonary_med', 'flag_pulmonary_med_t2e',
                         'flag_pulmonary_dx_med', 'flag_pulmonary_dx_med_t2e'
                         ]  # flag is used to include cardiology at followup based on different criteria and their time to event

    candidate_df[add_flag_baseline] = candidate_df.apply(lambda x: (0, 0, 0), axis=1, result_type='expand')
    candidate_df[add_flag_followup] = candidate_df.apply(lambda x: (0, 0, 9999, 0, 9999, 0, 9999), axis=1,result_type='expand')  # set 9999 as default time to event
    candidate_df['flag_pulmonary'] = 0  # flag is used for final pulmonary based on baseline pulmonary flag and followup pulmonary flag
    candidate_df['flag_pulmonary_t2e'] = 9999 # time to event based on min (flag_pulmonary_dx_t2e, flag_pulmonary_med_t2e, flag_pulmonary_dx_med_t2e) and maxfollowup, FOLLOWUP_RIGHT

    return candidate_df

def func_find_specific_diagnosis_code_date(candidate_dx_code_list, candidate_dx_date_list, candidate_type_list, specific_dx_list):
    Inpatient_IP_type = ['EI','IP','OS'] # inpatient
    Outpatient_OP_type = ['AV','OA','TH'] + ['ED'] # outpatient and emergency
    IP_code, OP_code = [],[]
    IP_date, OP_date = [],[]

    for i in range(len(candidate_dx_code_list)):
        i_code = candidate_dx_code_list[i]
        i_code = i_code.replace(".", '') # remove "."
        i_type = candidate_type_list[i]
        i_date = candidate_dx_date_list[i]

        if i_code in specific_dx_list:
            if i_type in Inpatient_IP_type:  # for inpatients
                IP_code.append(i_code)
                IP_date.append(i_date)

            if i_type in Outpatient_OP_type: # for outpatients
                OP_code.append(i_code)
                OP_date.append(i_date)

    return IP_code, IP_date, OP_code, OP_date

def func_find_specific_medication_code_date(candidate_med_code_list, candidate_med_date_list, specific_medication_rxcui_list):
    specific_rxcui = specific_medication_rxcui_list
    object_med_code = []
    object_med_date = []
    for i in range(len(candidate_med_code_list)):
        if candidate_med_code_list[i] in specific_rxcui:
            object_med_code.append(candidate_med_code_list[i])
            object_med_date.append(candidate_med_date_list[i])

    return object_med_code, object_med_date

def _eligibility_dx_only(id_indexrecord, id_dx, func_is_in_followup, dx_CL2_list):
    print("Step: applying _eligibility_followup_dx", 'input cohorts size:', len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # patid, site, covid, index date, hospitalized, ventilation, criticalcare, maxfollowup, death
        pid, covid_flag, index_date = row['patid'],  row['covid'], row['index date']

        v_dx = id_dx.get(pid, []) # # obtain all diagnosis code of a patient,
        if v_dx:
            dx_code_list = []
            dx_date_list = []
            dx_type_list = []
            for r in v_dx: # r:  time, dx code, ICD10, encounter type, such as (Timestamp('2020-03-02 00:00:00'), 'I35.0', '10', 'AV')
                dx_date, dx_code, dx_type = r[0], r[1], r[3]
                if func_is_in_followup(dx_date, index_date):
                    dx_code_list.append(dx_code)
                    dx_date_list.append(dx_date)
                    dx_type_list.append(dx_type)

            if dx_code_list:
                IP_code_list, IP_date_list, OP_code_list, OP_date_list = func_find_specific_diagnosis_code_date(dx_code_list, dx_date_list, dx_type_list, dx_CL2_list)

                OP_date_list = list(set(OP_date_list)) # remove, same time, for checking OP at different days (30 day apart)
                OP_date_list.sort() # rank

                IP_flag = False
                OP_flag = False

                if IP_code_list:  # at least 1 IP
                    IP_flag = True
                    id_indexrecord.loc[i, 'flag_pulmonary_dx'] = 1
                    id_indexrecord.loc[i,'flag_pulmonary_dx_t2e'] = (min(IP_date_list) - index_date).days # add time to event

                if len(OP_date_list)>=2:# OP_date_list has ranked, at least 2 OP at different days
                    OP_apart = (OP_date_list[-1] - OP_date_list[0]).days
                    if OP_apart>OP_apart_threshold: # (30 day apart)
                        OP_flag = True
                        id_indexrecord.loc[i, 'flag_pulmonary_dx'] = 1
                        OP_time2event = OP_date_list[-1]
                        # search time2event
                        j_start = OP_date_list[0]
                        for j in OP_date_list:
                            if (j - j_start).days>OP_apart_threshold:
                                OP_time2event = j
                                id_indexrecord.loc[i, 'flag_pulmonary_dx_t2e'] = (OP_time2event - index_date).days # add time to event, Note that, OP_date_list is ranked, and use the second date as time to event; Note that, OP_date_list[-1] should be adjusted.
                                break

                if IP_flag and OP_flag:
                    id_indexrecord.loc[i, 'flag_pulmonary_dx'] = 1
                    id_indexrecord.loc[i, 'flag_pulmonary_dx_t2e'] = min([(min(IP_date_list) - index_date).days, (OP_time2event - index_date).days])  # add min time of IP and OP as time to event
    # Summary
    total_pulmonary = id_indexrecord.loc[id_indexrecord['flag_pulmonary_dx']==1]
    covid_pos_pulmonary = id_indexrecord.loc[(id_indexrecord['flag_pulmonary_dx']==1)&(id_indexrecord['covid']==1)]
    covid_neg_pulmonary = id_indexrecord.loc[(id_indexrecord['flag_pulmonary_dx'] == 1) & (id_indexrecord['covid'] == 0)]
    print('...After  EC, total-pulmonary: {}\tcovid pos pulmonary: {}\tcovid neg pulmonary: {}'.format(len(total_pulmonary), len(covid_pos_pulmonary), len(covid_neg_pulmonary)))

    return id_indexrecord

def _eligibility_med_only(id_indexrecord, id_med, func_is_in_followup, med_rx1_list):
    print("Step: applying _eligibility_followup_med_only",'input cohorts size:', len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # patid, site, covid, index date, hospitalized, ventilation, criticalcare, maxfollowup, death
        pid, covid_flag, index_date = row['patid'],  row['covid'], row['index date']

        # considering RX1 prescription
        v_med = id_med.get(pid, [])
        if v_med:
            med_code_list = []
            med_date_list = []
            for r in v_med:
                med_date = r[0]
                med_code = r[1]
                if func_is_in_followup(med_date, index_date):
                    med_code_list.append(med_code)
                    med_date_list.append(med_date)
            if med_code_list:
                p_med_rx1_list, p_med_rx1_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_rx1_list)

                if len(p_med_rx1_list) >= 1: # med criteria
                    id_indexrecord.loc[i, 'flag_pulmonary_med'] = 1
                    id_indexrecord.loc[i, 'flag_pulmonary_med_t2e'] = (min(p_med_rx1_date) - index_date).days  # add time to event
    # Summary
    total_pulmonary = id_indexrecord.loc[id_indexrecord['flag_pulmonary_med'] == 1]
    covid_pos_pulmonary = id_indexrecord.loc[(id_indexrecord['flag_pulmonary_med'] == 1)&(id_indexrecord['covid'] == 1)]
    covid_neg_pulmonary = id_indexrecord.loc[(id_indexrecord['flag_pulmonary_med'] == 1) & (id_indexrecord['covid'] == 0)]
    print('...After  EC, total-pulmonary: {}\tcovid pos pulmonary: {}\tcovid neg pulmonary: {}'.format(len(total_pulmonary), len(covid_pos_pulmonary), len(covid_neg_pulmonary)))

    return id_indexrecord

def _eligibility_dx_med(id_indexrecord, id_dx, id_med, func_is_in_followup, dx_CL2_list, med_rx1_list, med_rx2_list):
    print("Step: applying _eligibility_followup_dx", 'input cohorts size:', len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # patid, site, covid, index date, hospitalized, ventilation, criticalcare, maxfollowup, death
        pid, covid_flag, index_date = row['patid'],  row['covid'], row['index date']

        flag_dx = False
        flag_med = False

        v_dx = id_dx.get(pid, []) # # obtain all diagnosis code of a patient,
        if v_dx:
            dx_code_list = []
            dx_date_list = []
            dx_type_list = []
            for r in v_dx: # r:  time, dx code, ICD10, encounter type, such as (Timestamp('2020-03-02 00:00:00'), 'I35.0', '10', 'AV')
                dx_date, dx_code, dx_type = r[0], r[1], r[3]
                if func_is_in_followup(dx_date, index_date):
                    dx_code_list.append(dx_code)
                    dx_date_list.append(dx_date)
                    dx_type_list.append(dx_type)

            if dx_code_list:
                IP_code_list, IP_date_list, OP_code_list, OP_date_list = func_find_specific_diagnosis_code_date(dx_code_list, dx_date_list, dx_type_list, dx_CL2_list)
                if len(IP_code_list+OP_code_list)>0:  # at least 1 IP
                    flag_dx = True
                    min_dx_date = min(IP_date_list + OP_date_list)

        # considering RX1 or RX2 prescription
        v_med = id_med.get(pid, [])
        if v_med:
            med_code_list = []
            med_date_list = []
            for r in v_med:
                med_date = r[0]
                med_code = r[1]
                if func_is_in_followup(med_date, index_date):
                    med_code_list.append(med_code)
                    med_date_list.append(med_date)
            if med_code_list:
                p_med_rx1_rx2_list, p_med_rx1_rx2_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_rx1_list+med_rx2_list)

                if len(p_med_rx1_rx2_list) >0: # med criteria
                    flag_med = True
                    min_med_date = min(p_med_rx1_rx2_date)

        if flag_dx and flag_med:
            id_indexrecord.loc[i, 'flag_pulmonary_dx_med'] = 1
            id_indexrecord.loc[i, 'flag_pulmonary_dx_med_t2e'] = (max([min_dx_date, min_med_date]) - index_date).days  # add time to event

    # Summary
    total_pulmonary = id_indexrecord.loc[id_indexrecord['flag_pulmonary_dx_med']==1]
    covid_pos_pulmonary = id_indexrecord.loc[(id_indexrecord['flag_pulmonary_dx_med']==1)&(id_indexrecord['covid']==1)]
    covid_neg_pulmonary = id_indexrecord.loc[(id_indexrecord['flag_pulmonary_dx_med'] == 1) & (id_indexrecord['covid'] == 0)]
    print('...After  EC, total-pulmonary: {}\tcovid pos pulmonary: {}\tcovid neg pulmonary: {}'.format(len(total_pulmonary), len(covid_pos_pulmonary), len(covid_neg_pulmonary)))

    return id_indexrecord

def label_dx_baseline(id_indexrecord, id_dx, func_is_in_baseline, dx_CL2_list):
    print("Step: label baseline_dx", 'input cohorts size:', len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # patid, site, covid, index date, hospitalized, ventilation, criticalcare, maxfollowup, death
        pid, covid_flag, index_date = row['patid'],  row['covid'], row['index date']

        v_dx = id_dx.get(pid, []) # obtain all diagnosis code of a patient,
        if v_dx:
            dx_code_list = []
            dx_date_list = []
            dx_type_list = []
            for r in v_dx: # r:  time, dx code, ICD10, encounter type, such as (Timestamp('2020-03-02 00:00:00'), 'I35.0', '10', 'AV')
                dx_date, dx_code, dx_type = r[0], r[1], r[3]
                if func_is_in_baseline(dx_date, index_date):
                    dx_code_list.append(dx_code)
                    dx_date_list.append(dx_date)
                    dx_type_list.append(dx_type)

            if dx_code_list:
                IP_code_list, IP_date_list, OP_code_list, OP_date_list = func_find_specific_diagnosis_code_date(dx_code_list, dx_date_list, dx_type_list, dx_CL2_list)

                if len(IP_code_list + OP_code_list)>0:  # at least 1 IP or OP
                    id_indexrecord.loc[i, 'flag_baseline_dx'] = 1

    # Summary
    total_cardiology = id_indexrecord.loc[id_indexrecord['flag_baseline_dx']==1]
    covid_pos_cardiology = id_indexrecord.loc[(id_indexrecord['flag_baseline_dx']==1)&(id_indexrecord['covid']==1)]
    covid_neg_cardiology = id_indexrecord.loc[(id_indexrecord['flag_baseline_dx'] == 1) & (id_indexrecord['covid'] == 0)]
    print('...baseline after  EC, total-cardiology: {}\tcovid pos cardiology: {}\tcovid neg cardiology: {}'.format(len(total_cardiology), len(covid_pos_cardiology), len(covid_neg_cardiology)))

    return id_indexrecord

def label_med_baseline(id_indexrecord, id_med, func_is_in_baseline, med_rx1_list):
    print("Step: applying label med baseline",'input cohorts size:', len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # patid, site, covid, index date, hospitalized, ventilation, criticalcare, maxfollowup, death
        pid, covid_flag, index_date = row['patid'],  row['covid'], row['index date']

        # considering RX1 prescription
        v_med = id_med.get(pid, [])
        if v_med:
            med_code_list = []
            med_date_list = []
            for r in v_med:
                med_date = r[0]
                med_code = r[1]
                if func_is_in_baseline(med_date, index_date):
                    med_code_list.append(med_code)
                    med_date_list.append(med_date)
            if med_code_list:
                p_med_rx1_list, p_med_rx1_date = func_find_specific_medication_code_date(med_code_list, med_date_list, med_rx1_list)

                if len(p_med_rx1_list) >= 1: # med criteria
                    id_indexrecord.loc[i, 'flag_baseline_med'] = 1

    # Summary
    total_pulmonary = id_indexrecord.loc[id_indexrecord['flag_baseline_med'] == 1]
    covid_pos_pulmonary = id_indexrecord.loc[(id_indexrecord['flag_baseline_med'] == 1)&(id_indexrecord['covid'] == 1)]
    covid_neg_pulmonary = id_indexrecord.loc[(id_indexrecord['flag_baseline_med'] == 1) & (id_indexrecord['covid'] == 0)]
    print('...  baseline after EC, total-pulmonary: {}\tcovid pos pulmonary: {}\tcovid neg pulmonary: {}'.format(len(total_pulmonary), len(covid_pos_pulmonary), len(covid_neg_pulmonary)))

    return id_indexrecord

def determine_time2event(id_indexrecord):
    print("Step: adding time to event:", len(id_indexrecord))
    for i, row in id_indexrecord.iterrows():
        # patid, site, covid, index date, hospitalized, ventilation, criticalcare, maxfollowup, death
        pid = row['patid']

        maxfollowup = id_indexrecord.loc[i,'maxfollowup']
        flag_pulmonary = id_indexrecord.loc[i,'flag_pulmonary']
        flag_pulmonary_dx_t2e = id_indexrecord.loc[i,'flag_pulmonary_dx_t2e']
        flag_pulmonary_med_t2e = id_indexrecord.loc[i,'flag_pulmonary_med_t2e']
        flag_pulmonary_dx_med_t2e = id_indexrecord.loc[i,'flag_pulmonary_dx_med_t2e']

        if flag_pulmonary ==1:
            id_indexrecord.loc[i, 'flag_pulmonary_t2e'] = min([flag_pulmonary_dx_t2e, flag_pulmonary_med_t2e, flag_pulmonary_dx_med_t2e])
        else:
            id_indexrecord.loc[i, 'flag_pulmonary_t2e'] = min([maxfollowup, FOLLOWUP_RIGHT])

    return id_indexrecord

def print_incidence_table(site, subcohort, incident_cohort, result_path, flag_disease):
    # print incidence tables for each site
    for i in [0, 1, 2]:
        print('**************************************************************')
        if i == 0:
            print('the result of non-hospitalized:')
            cohort = subcohort.loc[subcohort['hospitalized'] == 0]
            cohort_name = 'non-hospitalized'
        if i == 1:
            print('the result of hospitalized:')
            cohort = subcohort.loc[subcohort['hospitalized'] == 1]
            cohort_name = 'hospitalized'
        if i == 2:
            print('the result of all patients:')
            cohort = subcohort
            cohort_name = 'all_patients'

        total = len(cohort)
        inc_neg = len(cohort.loc[(cohort[flag_disease] == 1) & (cohort['covid'] == 0)])
        inc_pos = len(cohort.loc[(cohort[flag_disease] == 1) & (cohort['covid'] == 1)])
        inc = len(cohort.loc[cohort[flag_disease] == 1])

        non_inc_neg = len(cohort.loc[(cohort[flag_disease] == 0) & (cohort['covid'] == 0)])
        non_inc_pos = len(cohort.loc[(cohort[flag_disease] == 0) & (cohort['covid'] == 1)])
        non_inc = len(cohort.loc[cohort[flag_disease] == 0])

        neg = len(cohort.loc[cohort['covid'] == 0])
        pos = len(cohort.loc[cohort['covid'] == 1])
        inc_neg_per_1 = float('%.1f' % ((inc_neg / neg) * 100))
        inc_pos_per_1 = float('%.1f' % ((inc_pos / pos) * 100))
        inc_per_1 = float('%.1f' % ((inc / total) * 100))
        burden = inc_pos_per_1 - inc_neg_per_1

        print('total', total)
        print('inc_neg_per_100', inc_neg_per_1)
        print('inc_pos_per_100', inc_pos_per_1)
        print('inc_per_100', inc_per_1)

        print('inc_neg', inc_neg)
        print('inc_pos', inc_pos)
        print('inc', inc)

        print('non_inc_neg', non_inc_neg)
        print('non_inc_pos', non_inc_pos)
        print('non_inc', non_inc)

        print('pos', pos)
        print('neg', neg)

        incident_result = pd.DataFrame(index=['Negative', 'Positive', 'Total','Burden'],
                                       columns=['No', 'Yes', 'Total', 'Per 100'])
        incident_result.loc['Negative'] = [non_inc_neg, inc_neg, neg, inc_neg_per_1]
        incident_result.loc['Positive'] = [non_inc_pos, inc_pos, pos, inc_pos_per_1]
        incident_result.loc['Total'] = [non_inc, inc, total, inc_per_1]
        incident_result.loc['Burden'] = [np.nan, np.nan, np.nan, burden]

        print('incident_result in site:', site)
        print(incident_result)

        incident_result.to_csv(result_path + site + '_' + cohort_name + '_incidence.csv', index=True)

    # print incidence tables for all sites
    if site == 'vumc':
        site_name = 'all_sites'

        for i in [0, 1, 2]:
            print('**************************************************************')
            if i == 0:
                print('the result of non-hospitalized:')
                cohort = incident_cohort.loc[incident_cohort['hospitalized'] == 0]
                cohort_name = 'non-hospitalized'
            if i == 1:
                print('the result of hospitalized:')
                cohort = incident_cohort.loc[incident_cohort['hospitalized'] == 1]
                cohort_name = 'hospitalized'
            if i == 2:
                print('the result of all patients:')
                cohort = incident_cohort
                cohort_name = 'all_patients'

            total = len(cohort)
            inc_neg = len(cohort.loc[(cohort[flag_disease] == 1) & (cohort['covid'] == 0)])
            inc_pos = len(cohort.loc[(cohort[flag_disease] == 1) & (cohort['covid'] == 1)])
            inc = len(cohort.loc[cohort[flag_disease] == 1])

            non_inc_neg = len(cohort.loc[(cohort[flag_disease] == 0) & (cohort['covid'] == 0)])
            non_inc_pos = len(cohort.loc[(cohort[flag_disease] == 0) & (cohort['covid'] == 1)])
            non_inc = len(cohort.loc[cohort[flag_disease] == 0])

            neg = len(cohort.loc[cohort['covid'] == 0])
            pos = len(cohort.loc[cohort['covid'] == 1])
            inc_neg_per_1 = float('%.1f' % ((inc_neg / neg) * 100))
            inc_pos_per_1 = float('%.1f' % ((inc_pos / pos) * 100))
            inc_per_1 = float('%.1f' % ((inc / total) * 100))
            burden = inc_pos_per_1 - inc_neg_per_1

            print('total', total)
            print('inc_neg_per_100', inc_neg_per_1)
            print('inc_pos_per_100', inc_pos_per_1)
            print('inc_per_100', inc_per_1)

            print('inc_neg', inc_neg)
            print('inc_pos', inc_pos)
            print('inc', inc)

            print('non_inc_neg', non_inc_neg)
            print('non_inc_pos', non_inc_pos)
            print('non_inc', non_inc)

            print('pos', pos)
            print('neg', neg)

            incident_result = pd.DataFrame(index=['Negative', 'Positive', 'Total','Burden'],
                                           columns=['No', 'Yes', 'Total', 'Per 100'])
            incident_result.loc['Negative'] = [non_inc_neg, inc_neg, neg, inc_neg_per_1]
            incident_result.loc['Positive'] = [non_inc_pos, inc_pos, pos, inc_pos_per_1]
            incident_result.loc['Total'] = [non_inc, inc, total, inc_per_1]
            incident_result.loc['Burden'] = [np.nan, np.nan, np.nan, burden]

            print('incident_result in site:', site_name)  # using site name, not site
            print(incident_result)

            incident_result.to_csv(result_path + site_name + '_' + cohort_name + '_incidence.csv', index=True)

if __name__ == '__main__':
    start_time = time.time()

    flag_server = 0 # 1 run on server, 0 run on local

    if flag_server == 1:
        COVID_data_path = 'D:\\chz4001\\PycharmProjects\\pasc_phenotype\\data\\recover\\output\\'
        result_path = 'D:\\ZhenxingXu\\results\\pulmonary_incidence\\'
        pulmonary_CL_file = 'D:\\ZhenxingXu\\data\\PASC\\'
        ICD9_10_path = pulmonary_CL_file

        site_s = ['columbia', 'lsu', 'mcw', 'miami', 'michigan'] + \
                 ['montefiore', 'mshs', 'nch', 'nebraska', 'nyu', 'ochsner', 'pitt', 'psu'] + \
                 ['temple', 'ucsf', 'ufh', 'usf', 'utah', 'utsw', 'vumc', 'wcm']

    else:
        COVID_data_path = '/Users/xuzhenxing/Documents/PASC/data_old/RECOVER/'  # data for pulmonary
        result_path = '/Users/xuzhenxing/Documents/PASC/Cardiology_CP_and_Pulmonary_CP/pulmonary/pulmonary_incidence_results/'
        pulmonary_CL_file = '/Users/xuzhenxing/Documents/PASC/Cardiology_CP_and_Pulmonary_CP/'
        ICD9_10_path = '/Users/xuzhenxing/Documents/PASC/ICD9_10/'

        site_s = ['columbia']  # , 'mshs',  'nyu', 'wcm'

    dx_1a, med_rx1, med_rx2 = read_code_list(pulmonary_CL_file)
    icd9_icd10, icd_ccsr = read_dx_mapping_function(ICD9_10_path)
    # build a dataframe for saving results
    incident_cohort = pd.DataFrame()

    for site in site_s:
        print('site:', site)
        # read covid data
        if flag_server==1:
            id_indexrecord, id_diagnosis, id_medication = func_read_ori_data_server(COVID_data_path, site, icd9_icd10, icd_ccsr)
        else:
            id_indexrecord, id_diagnosis, id_medication = func_read_ori_data(COVID_data_path, site, icd9_icd10, icd_ccsr)

        # set flag
        id_indexrecord = func_set_pulm_flag(id_indexrecord)

        # pulmonary--3 criteria
        # ## EC--dx_only
        id_indexrecord = _eligibility_dx_only(id_indexrecord, id_diagnosis, _is_in_followup, dx_1a)
        # ## EC--med_only
        id_indexrecord = _eligibility_med_only(id_indexrecord, id_medication, _is_in_followup, med_rx1)
        # ## EC--dx_med
        id_indexrecord = _eligibility_dx_med(id_indexrecord, id_diagnosis, id_medication, _is_in_followup, dx_1a, med_rx1, med_rx2)

        # # ## baseline exclusion
        # # label dx
        id_indexrecord = label_dx_baseline(id_indexrecord, id_diagnosis, _is_in_baseline, dx_1a)
        # label med
        id_indexrecord = label_med_baseline(id_indexrecord, id_medication, _is_in_baseline, med_rx1)

        # baseline pulmonary
        id_indexrecord['flag_baseline_pulmonary'] = id_indexrecord['flag_baseline_dx'] + id_indexrecord['flag_baseline_med']
        id_indexrecord.loc[id_indexrecord['flag_baseline_pulmonary']>=1,'flag_baseline_pulmonary'] = 1

        # followup pulmonary
        id_indexrecord['flag_followup_pulmonary'] = id_indexrecord['flag_pulmonary_dx'] + id_indexrecord['flag_pulmonary_med'] + id_indexrecord['flag_pulmonary_dx_med']
        id_indexrecord.loc[id_indexrecord['flag_followup_pulmonary'] >= 1, 'flag_followup_pulmonary'] = 1
        #
        # final pulmonary based on baseline and followup pulmonary
        id_indexrecord.loc[(id_indexrecord['flag_followup_pulmonary']==1)&(id_indexrecord['flag_baseline_pulmonary']==0),'flag_pulmonary'] = 1
        #
        # determine time to event
        id_indexrecord = determine_time2event(id_indexrecord)

        # subcohort = id_indexrecord
        # incident_cohort = incident_cohort.append(id_indexrecord, ignore_index=True)
        # # print incidence tables for each site
        # flag_disease = 'flag_pulmonary'
        # print_incidence_table(site, subcohort, incident_cohort, result_path, flag_disease)

        # save final results to csv
        id_indexrecord.to_csv(result_path + 'pulmonary_incidence' + '_' + site + '.csv', index=False) # flag_disease = 'flag_pulmonary'
        print('len(pulmonary_incidence)', len(id_indexrecord))

    print('Done! Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

