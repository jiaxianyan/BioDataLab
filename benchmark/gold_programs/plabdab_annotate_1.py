import numpy as np
import pandas as pd
from anarci import anarci


SCHEME_BORDERS = {
               # Start coordinates
               # CDR1, FR2, CDR2, FR3, CDR3, FR4
         'imgt': [27,  39,  56,   66,  105,  118, 129],
      'kabat_H': [31,  36,  50,   66,  95,   103, 114],
      'kabat_K': [24,  35,  50,   57,  89,    98, 108],
      'kabat_L': [24,  35,  50,   57,  89,    98, 108],
    'chothia_H': [26,  33,  52,   57,  95,   103, 114],
    'chothia_K': [24,  35,  50,   57,  89,    98, 108],
    'chothia_L': [24,  35,  50,   57,  89,    98, 108],
      'north_H': [23,  36,  50,   59,  93,   103, 114],
      'north_K': [24,  35,  49,   57,  89,    98, 108],
      'north_L': [24,  35,  49,   57,  89,    98, 108],
}

# { scheme -> { region -> list of position numbers } }
SCHEME_REGIONS = {
    scheme: {
        'FR1': list(range(1, borders[0])),
        'CDR1': list(range(borders[0], borders[1])),
        'FR2': list(range(borders[1], borders[2])),
        'CDR2': list(range(borders[2], borders[3])),
        'FR3': list(range(borders[3], borders[4])),
        'CDR3': list(range(borders[4], borders[5])),
        'FR4': list(range(borders[5], borders[6])),
    } for scheme, borders in SCHEME_BORDERS.items()
}

def check_cdr_integrity(
    record,
    scheme: str = 'imgt'
):
    seq = str(record['GBSeq_sequence'])
    seqid = str(record['GBSeq_accession-version'])
    results = anarci([(seqid, seq.upper())], scheme=scheme, output=False)

    # Unpack the results. We get three lists
    numbering, alignment_details, hit_tables = results
    if numbering[0] is None:
        return False
    
    domain_numbering, start_index, end_index = numbering[0][0]
    chain_type = alignment_details[0][0]['chain_type']
    pos_ins_list, aa_list = list(zip(*domain_numbering))
    if record['numbered'] != ''.join(aa_list): 
        # print(f'row {record.name} {seqid} numbered sequences:')
        # print(record['numbered'])
        # print(''.join(aa_list))
        return False
    else:
        aa_seq = ''.join(aa_list)
        cdr1_seq = ''.join([aa_seq[i] for i in SCHEME_REGIONS['imgt']['CDR1'] if i <= (len(aa_seq) - 1)]).replace('-', '')
        cdr2_seq = ''.join([aa_seq[i] for i in SCHEME_REGIONS['imgt']['CDR2'] if i <= (len(aa_seq) - 1)]).replace('-', '')
        cdr3_seq = ''.join([aa_seq[i] for i in SCHEME_REGIONS['imgt']['CDR3'] if i <= (len(aa_seq) - 1)]).replace('-', '')
        if f'{len(cdr1_seq) - 1}_{len(cdr2_seq)}_{len(cdr3_seq)}' == record['cdr_lengths']:
            if chain_type in ['L', 'K']:
                chain_type = 'L'
            assert chain_type == record['chain'], f'{chain_type} != {record["chain"]}?'
            # print(f'row {record.name} {seqid} cdr length:')
            # print(f'{len(cdr1_seq)}_{len(cdr2_seq)}_{len(cdr3_seq)}', record['cdr_lengths'])
            return f'{len(cdr1_seq)}_{len(cdr2_seq)}_{len(cdr3_seq)}'
    
    return False
    
def main(
    ref_csv_file: str = 'benchmark/dataset/PLAbDab/unpaired_sequences.csv',
    export_fasta_to: str = 'benchmark/dataset/PLAbDab/plabdab_annotate_1_input.fasta',
    ref_csv_to: str = 'benchmark/gold_results/plabdab_annotate_1.csv',
    random_sample_num: int = 500,
):
    print('Read PLAbDab dataset unpaired sequence CSV file ....')
    df = pd.read_csv(ref_csv_file)
    print(f'Random sample {random_sample_num} sequences and try to find strict data for benchmark ....')
    sub_df = df.iloc[np.random.permutation(range(len(df)))[:random_sample_num]].copy()
    
    sub_df['fn_ret'] = sub_df.apply(check_cdr_integrity, axis = 1)
    
    sub_df_v2 = sub_df.loc[(sub_df['fn_ret'] != False)].copy()
    sub_df_v2['Raw_Sequence'] = sub_df_v2['GBSeq_sequence'].apply(lambda x: x.upper())
    sub_df_v2 = sub_df_v2[['GBSeq_accession-version', 'Raw_Sequence', 'numbered', 'fn_ret', 'chain']]
    sub_df_v2 = sub_df_v2.rename(
        {
            'GBSeq_accession-version': 'Accession',
            'numbered': 'Numbered_Sequence',
            'fn_ret': 'CDR_Lengths',
            'chain': 'Chain_Type'
        },
        axis = 1,
    )
    sub_df_v2 = sub_df_v2.drop_duplicates(subset='Accession', keep="first")
    print(f'Found {len(sub_df_v2)}/{random_sample_num} unique sequence for benchmark ....')
    
    fasta_content = []
    for rid, row in sub_df_v2.iterrows():
        accid = row['Accession']
        raw_sequence = row['Raw_Sequence']
        fasta_content.append(f'>{accid}\n{raw_sequence}')
    
    export_fasta_to = Path(export_fasta_to)
    export_fasta_to.parent.mkdir(exist_ok=True, parents=True)
    export_fasta_to.write_text('\n'.join(fasta_content))
    print(f'Export input fasta file to {export_fasta_to}. Done!')
    
    ref_csv_to = Path(ref_csv_to)
    ref_csv_to.parent.mkdir(exist_ok=True, parents=True)
    sub_df_v2.to_csv(ref_csv_to, index = False)
    print(f'Export reference csv file to {ref_csv_to}. Done!')
    print('All done!')

if __name__ == "__main__":
    main()