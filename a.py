import pandas as pd
import numpy as np
import sys, shutil,os

# /persist_disk/mocks_and_isolates/eu_20220817-120350_mocks_and_isolates/isolates/220816-210858/220816-210858
#python a.py 220816-210858.k2.report 0.01 220816-210858.k2.out  ../190320035227.2.1.1/190320035227.2.1.1.R1.filtered.fa ../190320035227.2.1.1/190320035227.2.1.1.R2.filtered.fa 190320035227.2.1.1

def report_parser(x, y):
    """
    Parse reports and construct a list for mtb tabel backbone
    """
    # ~ self.LOGGER.info(f"Parsing {self.report_file} file")
    l = []
    with open(x) as fin:
        for line in fin:
            if int(line.split('\t')[2]) >= 3 and int(line.split('\t')[4]) not in [0,1]:
                taxid = int(line.split('\t')[4])
                counts = int(line.split('\t')[2])
                out = f'{taxid}\t{counts}'
                l.append(out.split('\t'))
                    # ~ self.info_list.append(out.split('\t'))
    
    df = pd.DataFrame(l, columns=['taxid', 'counts'])
    df['taxid'] = df['taxid'].astype(int)
    df['counts'] = df['counts'].astype(int)
    df['perc']= (df['counts']/df['counts'].sum())*100
    print(df)
    taxid_list = df['taxid'].loc[df['perc'] <= float(y)].to_list()
    return taxid_list

tax_list_bellow_perc_cutoff = report_parser(sys.argv[1], sys.argv[2])

def get_read_id(z):
    tax_list_bellow_perc_cutoff = report_parser(sys.argv[1], sys.argv[2])
    c = 1
    info_list =[]
    with open(z) as fin:
        for line in fin:
            _id = line.split('\t')[1]
            tax = line.split('\t')[2]
           
            for i in   tax_list_bellow_perc_cutoff:
                if int(i) == int(tax):
                    out =f'{tax}\t{_id}\toligotype_{c}_{tax}'
                    c+=1
                    info_list.append(out.split('\t'))
                     
    df = pd.DataFrame(info_list, columns=['original_taxid', 'read_id', 'oligo_head'])
    df_outname = sys.argv[4].replace("R1.filtered.fa", 'id_equivalence.tsv')
    df.to_csv(df_outname, sep='\t', index=False)
    return df


def get_fasta_from_id(r1_fasta,r2_fasta):
        r1_dict, r2_dict = {}, {}
        with open(r1_fasta) as fin:
            for line in fin:
                if line.startswith('>'):
                    read_head = line.split(' ')[0].strip('>')
                    #index_read = line.split(' ')[1].strip('\n')
                else:
                    seq = line
                    r1_dict[read_head] = seq

        with open(r2_fasta) as fin:
            for line in fin:
                if line.startswith('>'):
                    read_head = line.split(' ')[0].strip('>')
                    #index_read = line.split(' ')[1].strip('\n')
                else:
                    seq = line
                    r2_dict[read_head] = seq
        return r1_dict, r2_dict
            
def generate_oligo_fasta():
    
    tax = get_read_id(sys.argv[3])
    
    out_r1_fname = sys.argv[4].replace('filtered','cutoff')
    out_r2_fname = sys.argv[5].replace('filtered','cutoff')
    
    out_r1_mapname = sys.argv[4].replace('filtered.fa','cutoff.map.tsv')
    out_r2_mapname = sys.argv[5].replace('filtered.fa','cutoff.map.tsv')
    
    out_r1 = open(out_r1_fname, 'w')
    out_r2 = open(out_r2_fname, 'w')
    
    out_r1m = open(out_r1_mapname, 'w')
    out_r2m = open(out_r2_mapname, 'w')
    
    r1,r2 = get_fasta_from_id(sys.argv[4],sys.argv[5])
    for index, row in tax.iterrows():
        _r1_head = '>' + row['oligo_head'].replace('oligotype', 'oligotypeR1')
        _r1_seq = r1.get(row['read_id'])
        _r1 = f'{_r1_head}\n{_r1_seq}'
        map_r1 = f'{_r1_head.strip(">")}\t{sys.argv[6]}_{_r1_head.split(">oligotype")[1]}\n'
        
        _r2_head = '>' + row['oligo_head'].replace('oligotype', 'oligotypeR2')
        _r2_seq = r2.get(row['read_id'])
        _r2 = f'{_r2_head}\n{_r2_seq}'
        map_r2 = f'{_r2_head.strip(">")}\t{sys.argv[6]}_{_r2_head.split(">oligotype")[1]}\n'
        
        out_r1.writelines(_r1)
        out_r2.writelines(_r2)
        
        
        out_r1m.writelines(map_r1)
        out_r2m.writelines(map_r2)
        



def concat_files():
    in_r1_fname = sys.argv[4].replace('filtered','cutoff')
    in_r2_fname = sys.argv[5].replace('filtered','cutoff')
    out_fname  = in_r1_fname.replace('R1.cutoff.fa', 'cutoff.fa')
    
    
    in_r1_mapname = sys.argv[4].replace('filtered.fa','cutoff.map.tsv')
    in_r2_mapname = sys.argv[5].replace('filtered.fa','cutoff.map.tsv')
    out_mapname  = in_r1_mapname.replace('R1.cutoff.map.tsv', 'cutoff.map.tsv')

    with open(out_fname,'wb') as wfd:
        for f in [in_r1_fname, in_r2_fname]:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    
    with open(out_mapname,'wb') as wfd:
        for f in [in_r1_mapname, in_r2_mapname]:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    
    os.remove(in_r1_fname)
    os.remove(in_r2_fname)
    os.remove(in_r1_mapname)
    os.remove(in_r2_mapname)


generate_oligo_fasta()

concat_files()
