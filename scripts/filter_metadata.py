# -*- coding: utf-8 -*-

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2020-03-24
# Last update: 2021-08-10


import pycountry_convert as pyCountry
import pycountry
from Bio import SeqIO
import pandas as pd
from epiweeks import Week
import time
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file genomes to be used")
    parser.add_argument("--metadata1", required=True, help="Metadata file from NextStrain")
    parser.add_argument("--metadata2", required=False, help="Custom lab metadata file")
    parser.add_argument("--filter", required=False, nargs='+', type=str,  help="List of filters for tagged rows in lab metadata")
    parser.add_argument("--output1", required=True, help="Filtered metadata file")
    parser.add_argument("--output2", required=True, help="Reformatted, final FASTA file")
    args = parser.parse_args()

    genomes = args.genomes
    metadata1 = args.metadata1
    metadata2 = args.metadata2
    filterby = args.filter
    output1 = args.output1
    output2 = args.output2

    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_variants/nextstrain/runX_20210617_filter/'
    # genomes = path + 'pre-analyses/temp_sequences.fasta'
    # metadata1 = path + 'pre-analyses/metadata_nextstrain.tsv'
    # metadata2 = path + 'pre-analyses/GLab_SC2_sequencing_data.xlsx'
    # filterby = ['caribe', 'test']
    # output1 = path + 'pre-analyses/metadata_filtered.tsv'
    # output2 = path + 'pre-analyses/sequences.fasta'

    # temporal boundaries
    today = time.strftime('%Y-%m-%d', time.gmtime())
    min_date = '2019-12-15'

    variants = {'B.1.1.7': 'Alpha (B.1.1.7)',
                'B.1.351': 'Beta (B.1.351)',
                'B.1.351.2': 'Beta (B.1.351.2)',
                'B.1.351.3': 'Beta (B.1.351.3)',
                'P.1': 'Gamma (P.1)',
                'P.1.1': 'Gamma (P.1.1)',
                'P.1.2': 'Gamma (P.1.2)',
                'B.1.617.2': 'Delta (B.1.617.2)',
                'AY.1': 'Delta (AY.1)',
                'AY.2': 'Delta (AY.2)',
                'AY.3': 'Delta (AY.3)',
                'AY.3.1': 'Delta (AY.3.1)',
                'AY.4': 'Delta (AY.4)',
                'AY.4.1': 'Delta (AY.4.1)',
                'AY.4.2': 'Delta (AY.4.2)',              
                'AY.4.3': 'Delta (AY.4.3)',
                'AY.4.4': 'Delta (AY.4.4)',
                'AY.4.5': 'Delta (AY.4.5)',
                'AY.5': 'Delta (AY.5)',
                'AY.5.1': 'Delta (AY.5.1)',
                'AY.5.2': 'Delta (AY.5.2)',
                'AY.5.3': 'Delta (AY.5.3)',
                'AY.5.4': 'Delta (AY.5.4)',
                'AY.6': 'Delta (AY.6)',
                'AY.7': 'Delta (AY.7)',
                'AY.7.1': 'Delta (AY.7.1)',
                'AY.7.2': 'Delta (AY.7.2)',
                'AY.8': 'Delta (AY.8)',
                'AY.9': 'Delta (AY.9)',
                'AY.9.1': 'Delta (AY.9.1)',
                'AY.9.2': 'Delta (AY.9.2)',
                'AY.9.2.1': 'Delta (AY.9.2.1)',               
                'AY.10': 'Delta (AY.10)',
                'AY.11': 'Delta (AY.11)',
                'AY.12': 'Delta (AY.12)',
                'AY.13': 'Delta (AY.13)',
                'AY.14': 'Delta (AY.14)',
                'AY.15': 'Delta (AY.15)',
                'AY.16': 'Delta (AY.16)',
                'AY.16.1': 'Delta (AY.16.1)',
                'AY.17': 'Delta (AY.17)',
                'AY.18': 'Delta (AY.18)',
                'AY.19': 'Delta (AY.19)',
                'AY.20': 'Delta (AY.20)',
                'AY.21': 'Delta (AY.21)',
                'AY.22': 'Delta (AY.22)',
                'AY.23': 'Delta (AY.23)',
                'AY.23.1': 'Delta (AY.23.1)',
                'AY.24': 'Delta (AY.24)',
                'AY.25': 'Delta (AY.25)',
                'AY.26': 'Delta (AY.26)',
                'AY.27': 'Delta (AY.27)',
                'AY.28': 'Delta (AY.28)',
                'AY.29': 'Delta (AY.29)',
                'AY.29.1': 'Delta (AY.29.1)',
                'AY.30': 'Delta (AY.30)',
                'AY.31': 'Delta (AY.31)',
                'AY.32': 'Delta (AY.32)',
                'AY.33': 'Delta (AY.33)',
                'AY.34': 'Delta (AY.34)',
                'AY.35': 'Delta (AY.35)',
                'AY.36': 'Delta (AY.36)',
                'AY.37': 'Delta (AY.37)',
                'AY.38': 'Delta (AY.38)',
                'AY.39': 'Delta (AY.39)',
                'AY.39.1': 'Delta (AY.39.1)',
                'AY.39.1.1': 'Delta (AY.39.1.1)',
                'AY.39.2': 'Delta (AY.39.2)',
                'AY.40': 'Delta (AY.40)',
                'AY.41': 'Delta (AY.41)',
                'AY.42': 'Delta (AY.42)',
                'AY.43': 'Delta (AY.43)',
                'AY.44': 'Delta (AY.44)',
                'AY.45': 'Delta (AY.45)',
                'AY.46': 'Delta (AY.46)',
                'AY.46.1': 'Delta (AY.46.1)',
                'AY.46.2': 'Delta (AY.46.2)',
                'AY.46.3': 'Delta (AY.46.3)',
                'AY.46.4': 'Delta (AY.46.4)',
                'AY.46.5': 'Delta (AY.46.5)',
                'AY.46.6': 'Delta (AY.46.6)',
                'AY.47': 'Delta (AY.47)',
                'AY.48': 'Delta (AY.48)',
                'AY.49': 'Delta (AY.49)',
                'AY.50': 'Delta (AY.50)',
                'AY.51': 'Delta (AY.51)',
                'AY.52': 'Delta (AY.52)',
                'AY.53': 'Delta (AY.53)',
                'AY.54': 'Delta (AY.54)',
                'AY.55': 'Delta (AY.55)',
                'AY.56': 'Delta (AY.56)',
                'AY.57': 'Delta (AY.57)',
                'AY.58': 'Delta (AY.58)',
                'AY.59': 'Delta (AY.59)',
                'AY.60': 'Delta (AY.60)',
                'AY.61': 'Delta (AY.61)',
                'AY.62': 'Delta (AY.62)',
                'AY.63': 'Delta (AY.63)',
                'AY.64': 'Delta (AY.64)',
                'AY.65': 'Delta (AY.65)',
                'AY.66': 'Delta (AY.66)',
                'AY.67': 'Delta (AY.67)',
                'AY.68': 'Delta (AY.68)',
                'AY.69': 'Delta (AY.69)',
                'AY.70': 'Delta (AY.70)',
                'AY.71': 'Delta (AY.71)',
                'AY.72': 'Delta (AY.72)',
                'AY.73': 'Delta (AY.73)',
                'AY.74': 'Delta (AY.74)',
                'AY.75': 'Delta (AY.75)',
                'AY.75.1': 'Delta (AY.75.1)',
                'AY.76': 'Delta (AY.76)',
                'AY.77': 'Delta (AY.77)',
                'AY.78': 'Delta (AY.78)',
                'AY.79': 'Delta (AY.79)',
                'AY.80': 'Delta (AY.80)',
                'AY.81': 'Delta (AY.81)',
                'AY.82': 'Delta (AY.82)',
                'AY.83': 'Delta (AY.83)',
                'AY.84': 'Delta (AY.84)',
                'AY.85': 'Delta (AY.85)',
                'AY.86': 'Delta (AY.86)',
                'AY.87': 'Delta (AY.87)',
                'AY.88': 'Delta (AY.88)',
                'AY.89': 'Delta (AY.89)',
                'AY.90': 'Delta (AY.90)',
                'AY.91': 'Delta (AY.91)',
                'AY.91.1': 'Delta (AY.91.1)',
                'AY.92': 'Delta (AY.92)',
                'AY.93': 'Delta (AY.93)',
                'AY.94': 'Delta (AY.94)',
                'AY.95': 'Delta (AY.95)',
                'AY.96': 'Delta (AY.96)',
                'AY.97': 'Delta (AY.97)',
                'AY.98': 'Delta (AY.98)',
                'AY.98.1': 'Delta (AY.98.1)',
                'AY.99': 'Delta (AY.99)',
                'AY.99.1': 'Delta (AY.99.1)',
                'AY.99.2': 'Delta (AY.99.2)',
                'AY.100': 'Delta (AY.100)',
                'AY.101': 'Delta (AY.101)',
                'AY.102': 'Delta (AY.102)',
                'AY.103': 'Delta (AY.103)',
                'AY.104': 'Delta (AY.104)',
                'AY.105': 'Delta (AY.105)',
                'AY.106': 'Delta (AY.106)',
                'AY.107': 'Delta (AY.107)',
                'AY.108': 'Delta (AY.108)',
                'AY.109': 'Delta (AY.109)',
                'AY.110': 'Delta (AY.110)',
                'AY.111': 'Delta (AY.111)',
                'AY.112': 'Delta (AY.112)',
                'AY.113': 'Delta (AY.113)',
                'AY.114': 'Delta (AY.114)',
                'AY.115': 'Delta (AY.115)',
                'AY.116': 'Delta (AY.116)',
                'AY.116.1': 'Delta (AY.116.1)',
                'AY.117': 'Delta (AY.117)', 
                'AY.118': 'Delta (AY.118)',    
                'AY.119': 'Delta (AY.119)',   
                'AY.120': 'Delta (AY.120)',   
                'AY.120.1': 'Delta (AY.120.1)', 
                'AY.120.2': 'Delta (AY.120.2)',
                'AY.120.2.1': 'Delta (AY.120.2.1)',
                'AY.121': 'Delta (AY.121)',
                'AY.121.1': 'Delta (AY.121.1)',
                'AY.122': 'Delta (AY.122)', 
                'AY.122.1': 'Delta (AY.122.1)',
                'AY.123': 'Delta (AY.123)',
                'AY.124': 'Delta (AY.124)', 
                'AY.125': 'Delta (AY.125)',
                'AY.126': 'Delta (AY.126)',
                'AY.127': 'Delta (AY.127)',
                'AY.128': 'Delta (AY.128)',
                'AY.129': 'Delta (AY.129',
                'AY.130': 'Delta (AY.130)',                                                                                                                                          
                'B.1.525': 'Eta (B.1.525)',
                'B.1.526': 'Iota (B.1.526)',
                'B.1.617.1': 'Kappa (B.1.617.1)',
                'C.37': 'Lambda (C.37)',
                'B.1.427': 'Epsilon (B.1.427/B.1.429)',
                'B.1.429': 'Epsilon (B.1.427/B.1.429)',
                'P.2': 'Zeta (P.2)',
                'B.1.621': 'Mu (B.1.621)',
                'B.1.621.1': 'Mu (B.1.621.1)',
                'BA.1': 'Omicron (BA.1)'
                }


    # get ISO alpha3 country codes

    isos =  {'US Virgin Islands':'VIR','Virgin Islands':'VIR','British Virgin Islands':'VGB','Curacao':'CUW','Northern Mariana Islands':'MNP',
            'Sint Maarten':'MAF','St Eustatius':'BES'}
    def get_iso(country):
        global isos
        if country not in isos.keys():
            try:
                isoCode = pyCountry.country_name_to_country_alpha3(country, cn_name_format="default")
                isos[country] = isoCode
            except:
                try:
                    isoCode = pycountry.countries.search_fuzzy(country)[0].alpha_3
                    isos[country] = isoCode
                except:
                    isos[country] = ''
        return isos[country]


    # create epiweek column
    def get_epiweeks(date):
        date = pd.to_datetime(date)
        epiweek = str(Week.fromdate(date, system="cdc"))  # get epiweeks
        epiweek = epiweek[:4] + '_' + 'EW' + epiweek[-2:]
        return epiweek


    # add state code
    us_state_abbrev = {
        'Alabama': 'AL',
        'Alaska': 'AK',
        'American Samoa': 'AS',
        'Arizona': 'AZ',
        'Arkansas': 'AR',
        'California': 'CA',
        'Colorado': 'CO',
        'Connecticut': 'CT',
        'Delaware': 'DE',
        'District of Columbia': 'DC',
        'Washington DC': 'DC',
        'Florida': 'FL',
        'Georgia': 'GA',
        'Guam': 'GU',
        'Hawaii': 'HI',
        'Idaho': 'ID',
        'Illinois': 'IL',
        'Indiana': 'IN',
        'Iowa': 'IA',
        'Kansas': 'KS',
        'Kentucky': 'KY',
        'Louisiana': 'LA',
        'Maine': 'ME',
        'Maryland': 'MD',
        'Massachusetts': 'MA',
        'Michigan': 'MI',
        'Minnesota': 'MN',
        'Mississippi': 'MS',
        'Missouri': 'MO',
        'Montana': 'MT',
        'Nebraska': 'NE',
        'Nevada': 'NV',
        'New Hampshire': 'NH',
        'New Jersey': 'NJ',
        'New Mexico': 'NM',
        'New York': 'NY',
        'North Carolina': 'NC',
        'North Dakota': 'ND',
        'Northern Mariana Islands': 'MP',
        'Ohio': 'OH',
        'Oklahoma': 'OK',
        'Oregon': 'OR',
        'Pennsylvania': 'PA',
        'Puerto Rico': 'PR',
        'Rhode Island': 'RI',
        'South Carolina': 'SC',
        'South Dakota': 'SD',
        'Tennessee': 'TN',
        'Texas': 'TX',
        'Utah': 'UT',
        'Vermont': 'VT',
        'Virgin Islands': 'VI',
        'Virginia': 'VA',
        'Washington': 'WA',
        'West Virginia': 'WV',
        'Wisconsin': 'WI',
        'Wyoming': 'WY'
    }

    # nextstrain metadata
    dfN = pd.read_csv(metadata1, encoding='utf-8', sep='\t', dtype='str')
    dfN['strain'] = dfN['strain'].replace('hCoV-19/', '')
    dfN.insert(4, 'iso', '')
    dfN.insert(1, 'category', '')
    dfN.fillna('', inplace=True)

    # add tag of variant category
    def variant_category(lineage):
        var_category = 'Other variants'
        for name in variants.keys():
            if lineage == name:
                var_category = variants[lineage]
        return var_category

    dfN['category'] = dfN['pango_lineage'].apply(lambda x: variant_category(x))


    list_columns = dfN.columns.values  # list of column in the original metadata file
    # print(dfN)

    # Lab genomes metadata
    dfE = pd.read_excel(metadata2, index_col=None, header=0, sheet_name=0,
                        converters={'Sample-ID': str, 'Collection-date': str})
    dfE.fillna('', inplace=True)

    dfE = dfE.rename(
        columns={'Sample-ID': 'id', 'Collection-date': 'date', 'Country': 'country', 'Division (state)': 'division',
                 'Location (county)': 'location', 'Country of exposure': 'country_exposure',
                 'State of exposure': 'division_exposure', 'Lineage': 'pango_lineage', 'Source': 'originating_lab',
                 'Filter': 'filter'})
    dfE['epiweek'] = ''

    # exclude rows with no ID
    if 'id' in dfE.columns.to_list():
        dfE = dfE[~dfE['id'].isin([''])]

    lab_sequences = dfE['id'].tolist()
    # exclude unwanted lab metadata row
    if len(filterby) > 0:
        print('\nFiltering metadata by category: ' + ', '.join(filterby) + '\n')
    dfL = pd.DataFrame(columns=dfE.columns.to_list())
    for value in filterby:
        dfF = dfE[dfE['filter'].isin([value])]  # batch inclusion of specific rows
        dfL = pd.concat([dfL, dfF]) # add filtered rows to dataframe with lab metadata

    # list of relevant genomes sequenced
    keep_only = dfL['id'].tolist()
    excluded = [id for id in lab_sequences if id not in keep_only]

    # create a dict of existing sequences
    sequences = {}
    for fasta in SeqIO.parse(open(genomes), 'fasta'):  # as fasta:
        id, seq = fasta.description, fasta.seq
        if id not in sequences.keys() and id not in excluded:
            sequences[id] = str(seq)

    # add inexistent columns
    for col in list_columns:
        if col not in dfL.columns:
            dfL[col] = ''

    # output dataframe
    outputDF = pd.DataFrame(columns=list_columns)
    found = []
    lab_label = {}
    metadata_issues = {}
    # process metadata from excel sheet
    for idx, row in dfL.iterrows():
        id = dfL.loc[idx, 'id'].replace('hCoV-19/', '')
        if id in sequences:
            dict_row = {}
            for col in list_columns:
                dict_row[col] = ''
                if col in row:
                    dict_row[col] = dfL.loc[idx, col].strip()  # add values to dictionary

            # check for missing geodata
            geodata = ['country'] # column
            for level in geodata:
                if len(dict_row[level]) < 1:
                    if id not in metadata_issues:
                        metadata_issues[id] = [level]
                    else:
                        metadata_issues[id].append(level)

            if dict_row['location'] in ['', None]:
                dict_row['location'] = dfL.loc[idx, 'location']

            collection_date = ''
            if len(str(dict_row['date'])) > 1 and 'X' not in dict_row['date']:
                collection_date = dict_row['date'].split(' ')[0].replace('.', '-').replace('/', '-')
                dict_row['date'] = collection_date
                # check is date is appropriate: not from the 'future', not older than 'min_date'
                if pd.to_datetime(today) < pd.to_datetime(collection_date) or pd.to_datetime(min_date) > pd.to_datetime(collection_date):
                    if id not in metadata_issues:
                        metadata_issues[id] = ['date']
                    else:
                        metadata_issues[id].append('date')
            else: # missing date
                if id not in metadata_issues:
                    metadata_issues[id] = ['date']
                else:
                    metadata_issues[id].append('date')

            # fix exposure
            columns_exposure = ['country_exposure', 'division_exposure']
            for level_exposure in columns_exposure:
                level = level_exposure.split('_')[0]
                dict_row[level_exposure] = dfL.loc[idx, level_exposure]
                if dict_row[level_exposure] in ['', None]:
                    if level_exposure == 'country_exposure':
                        dict_row[level_exposure] = dict_row[level]
                    else:
                        if dict_row['country_exposure'] != dfL.loc[idx, 'country']:
                            dict_row[level_exposure] = dict_row['country_exposure']
                        else:
                            dict_row[level_exposure] = dict_row[level]

            code = ''
            if dict_row['division'] in us_state_abbrev:
                code = us_state_abbrev[dict_row['division']] + '-'

            strain = dfL.loc[idx, 'country'].replace(' ', '') + '/' + code + dfL.loc[idx, 'id'] + '/' + collection_date.split('-')[0]

            # set the strain name
            dict_row['strain'] = strain
            dict_row['iso'] = get_iso(dict_row['country'])
            dict_row['originating_lab'] = dfL.loc[idx, 'originating_lab']
            dict_row['submitting_lab'] = 'Grubaugh Lab - Yale School of Public Health'
            dict_row['authors'] = 'GLab team'

            # add lineage
            lineage = ''
            if dfL.loc[idx, 'pango_lineage'] != '':
                lineage = dfL.loc[idx, 'pango_lineage']
            dict_row['pango_lineage'] = lineage

            # variant classication (VOI, VOC, VHC)
            dict_row['category'] = variant_category(lineage)

            # assign epiweek
            if len(dict_row['date']) > 0 and 'X' not in dict_row['date']:
                dict_row['epiweek'] = get_epiweeks(collection_date)
            else:
                dict_row['epiweek'] = ''

            # record sequence and metadata as found
            found.append(strain)
            if id not in metadata_issues.keys():
                lab_label[id] = strain
                outputDF = outputDF.append(dict_row, ignore_index=True)


    # process metadata from TSV
    dfN = dfN[dfN['strain'].isin(sequences.keys())]
    for idx, row in dfN.iterrows():
        strain = dfN.loc[idx, 'strain'].replace('hCoV-19/', '')
        if strain in sequences:
            if strain in outputDF['strain'].to_list():
                continue
            dict_row = {}
            date = ''
            for col in list_columns:
                if col == 'date':
                    date = dfN.loc[idx, col]
                dict_row[col] = ''
                if col in row:
                    dict_row[col] = dfN.loc[idx, col]

            # fix exposure
            columns_exposure = ['country_exposure', 'division_exposure']
            for level_exposure in columns_exposure:
                level = level_exposure.split('_')[0]
                if dict_row[level_exposure] in ['', None]:
                    dict_row[level_exposure] = dict_row[level]

            dict_row['iso'] = get_iso(dict_row['country'])
            dict_row['epiweek'] = get_epiweeks(date)
            found.append(strain)

            outputDF = outputDF.append(dict_row, ignore_index=True)

    # write new metadata files
    outputDF = outputDF.drop(columns=['region'])
    outputDF.to_csv(output1, sep='\t', index=False)

    # write sequence file
    exported = []
    with open(output2, 'w') as outfile2:
        # export new metadata lines
        for id, sequence in sequences.items():
            if id in lab_label and id not in metadata_issues.keys(): # export lab generated sequences
                if lab_label[id] not in exported:
                    entry = '>' + lab_label[id] + '\n' + sequence + '\n'
                    outfile2.write(entry)
                    print('* Exporting newly sequenced genome and metadata for ' + id)
                    exported.append(lab_label[id])
            else:  # export publicly available sequences
                if id not in exported and id in outputDF['strain'].tolist():
                    entry = '>' + id + '\n' + sequence + '\n'
                    outfile2.write(entry)
                    exported.append(id)

    if len(metadata_issues) > 0:
        print('\n\n### WARNINGS!\n')
        print('\nPlease check for metadata issues related to these samples and column (which will be otherwise ignored)\n')
        for id, columns in metadata_issues.items():
            print('\t- ' + id + ' (issues found at: ' + ', '.join(columns) + ')')

    print('\nMetadata file successfully reformatted and exported!\n')
