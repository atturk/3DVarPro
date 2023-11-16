from Bio import Entrez
from email_here import email
import Bio.PDB as BP
import pandas as pd
import re
from Bio.SeqUtils import IUPACData as aa


#cerca in clinvar le mutazioni per il gene e la condizione specificati
#restituisce dataframe con location mutazione e sostizuione AA in formato aawt, pos e aamut
def cvfetch(gene, pdb):
    record = cvrecord(gene)
    if type(record) == Warning:
        return Warning('**Error**: No mutations were found for the specified gene: "' + gene + '"')
    elif type(record) == list:
        df = record_to_df(record, pdb)
        if type(df) == Warning:
            return Warning("**Error**: The PDB file contains more than one chain/model. Currently, 3DVarPro supports single-chain PDB files only.")
        else:
            return df
     

#prende il record di cvfetch e lo trasforma in dataframe
def record_to_df(record, pdb):
    #parsing e controllo
    parser = BP.PDBParser()
    structure = parser.get_structure('proteina', pdb)
    model = structure[0]
    #Controlla se ha più di una catena
    if len(model) > 1:
        return Warning("**Error**: The PDB file contains more than one chain. Currently, 3DVarPro supports single-chain PDB files only.")
    else:
        chain = model['A']

        #setta il pattern che le mutazioni devono rispettare
        pattern = re.compile(r'^[A-Za-z]\d+[A-Za-z]$')
        #estrae dati da dizionario clinvar (record)
        df_data = []
        #per ogni entry dato che un record con una sola loc può avere più varianti
        for entry in record:
            variation_name = entry['variation_set'][0]['variation_name']
            trait_name = [element['trait_name'] for element in entry['trait_set'] if element['trait_name'] != 'not provided' and element['trait_name'] != 'not specified']
            mutations = entry['protein_change'].split(', ')
            clinical_significance = entry['clinical_significance']['description']
            #divide le mutazioni dal formato X123Y in [X, 123, Y]
            for mutation in mutations:
                #prima controlla che la mutazione rispetti il pattern
                if re.match(pattern, mutation):
                    aawt = mutation[0]
                    pos = ''.join(mutation[1:-1])
                    aamut = mutation[-1]
                    for residue in chain:
                        if residue.get_id()[1] == int(pos):
                            #ottiene bfactor e coordinate atomo CA per ogni residuo nella posizione della mutazione
                            bfactor = residue['CA'].get_bfactor()
                            coord=residue['CA'].get_coord()
                            df_data.append([variation_name, aawt, int(pos), aamut ,trait_name, clinical_significance, coord, bfactor])


        #crea dataframe
        df = pd.DataFrame(df_data, columns=['Location', 'WT AA','POS','MUT AA', 'Associated trait','Clinical Significance', 'Coordinates', 'Bfactor'])
        return df

#cerca in clinvar le mutazioni per il gene e la condizione specificati
#restituisce record (lista di dizionari) che poi va messo in record_to_df
def cvrecord(gene):

    # entrez email
    if email == "":
        raise ValueError("Please provide a valid email address in the email_here.py file.")
    else:
        Entrez.email = "italoruttico@gmail.com"

    #variabili
    database = "clinvar"
    max = 1000
    ricerca = gene + '''[gene] AND "mol cons missense"[filter]'''
    print(ricerca)

    #cerca
    handle = Entrez.esearch(db=database, 
                            term=ricerca,
                            retmax=max)
    record = Entrez.read(handle)
    handle.close()

    idlist = record["IdList"]

    if len(idlist) != 0:
        #ottieni i risultati
        handle = Entrez.esummary(db="clinvar", id=idlist,  rettype="text", retmode="tabular")
        record = Entrez.read(handle)
        handle.close()

        data = record.get('DocumentSummarySet').get('DocumentSummary')

        return data
    else:
        return Warning('**Error**: No mutations were found for the specified gene: "' + gene + '"')

def read_save_pdb(pdb_file):
    #legge file pdb caricato
    file_pdb = pdb_file.getvalue().decode("utf-8")

    #salva in proteina.pdb
    with open('proteina.pdb', 'w') as pdb:
        pdb.write(file_pdb)
    pdb_filename = 'proteina.pdb'
    return [pdb_filename, file_pdb]

def create_pdb(uploaded_file):
    pdb_decode = uploaded_file.getvalue().decode("utf-8")
        #salva in proteina.pdb
    with open('proteina.pdb', 'w') as f:
        f.write(pdb_decode)
    pdb_file = 'proteina.pdb'
    return pdb_file

def get_residues(pdb_file):
    parser = BP.PDBParser()
    data = parser.get_structure('mod',pdb_file)
    model = data.get_models()
    models = list(model)
    chains = list(models[0].get_chains())
    residues = list(chains[0].get_residues())
    return residues

def check_df(df, pdb_file):
    df2 = []
    residues = get_residues(pdb_file)
    for residue in residues:
        resn = aa.protein_letters_3to1_extended[residue.get_resname().title()]
        if [residue.get_id()[1], resn] in df[['POS', 'WT AA']].values.tolist():
            #seleziona riga del df che ha posizione e AA corrispondenti al residuo
            riga = df.loc[(df['POS'] == residue.get_id()[1]) & (df['WT AA'] == resn)]
            #crea un dataframe solo con le righe selezionate
            df2.append(riga)

    if len(df2) == 0:
        return Warning("**Error**: The filtered dataset is empty. Please check your filter settings.")
    else:
        df2 = pd.concat(df2)
    return df2