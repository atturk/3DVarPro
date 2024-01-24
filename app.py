from utils import *
import requests
import json
import time
import re
import streamlit as st
import pandas as pd
import Bio.PDB as BP
import py3Dmol as p3d
from stmol import *
from Bio.SeqUtils import IUPACData as aa
from unipressed import IdMappingClient

def disable():
    st.session_state["disabled"] = True

def showexample():
    st.session_state.ex = True

def sessionstatedfu_and_fileup():
    sessionstatedfu()
    fileuploaded()
    
def sessionstatedfu():
    st.session_state.dfu = None

def sessionstatedff():
    st.session_state.dff = None

def idonlysearch():
    st.session_state.idonly = True

def fileuploaded():
    st.session_state.fileup = True

def get_pdb():
        st.markdown('<h4>Upload the PDB file of the wild-type protein</h4>', unsafe_allow_html=True)
        uploaded_file = st.file_uploader(label='pdb uploader',type=['pdb'], label_visibility='collapsed', on_change=sessionstatedfu_and_fileup, disabled=st.session_state.disabled)
        if uploaded_file is not None:
            pdb_file = create_pdb(uploaded_file)
            st.write('File uploaded successfully')
            return pdb_file

def get_id():
    st.markdown('<h4>Insert the identifier </h4>', unsafe_allow_html=True)
    st.info("You can use the UniProtKB AC/ID or the Gene Name", icon='ℹ️')
    id = st.text_input('gene', label_visibility='collapsed', on_change=sessionstatedfu, disabled=st.session_state.disabled).upper()
    return id

def id_type(identifier):
    regex = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    if re.match(regex, identifier):
        return 'uniprot'
    else:
        return 'gene'

def id_to_genename_and_uniprotid(id):
    if id !='':
        identifier=id_type(id)

        if identifier == 'gene':
            request = IdMappingClient.submit(
                source="GeneCards", dest="UniProtKB", ids={id})
            time.sleep(2.0)
            res = list(request.each_result())
            if res == []:
                return Warning
            else:
                res=res[0].get('to')
                uniprotid=res
                gene=id
                return uniprotid, gene

        else:
            request = IdMappingClient.submit(
                source="UniProtKB_AC-ID", dest="GeneCards", ids={id})
            time.sleep(2.0)
            res = list(request.each_result())
            if res == []:
                return Warning
                st.warning("Identifier not valid")
                st.stop()
            else:
                res=res[0].get('to')
                uniprotid=id
                gene=res
                return uniprotid, gene
            
def get_af_pdb(uniprotid):
    url = f"https://alphafold.com/api/prediction/{uniprotid}"
    try:
        response = requests.get(url)
        pdburl = json.loads(response.text)[0].get('pdbUrl')
        response = requests.get(pdburl)
        with open('proteina.pdb', 'w') as f:
            f.write(response.text)  
        pdb_file = 'proteina.pdb'
        return pdb_file
    except:
        return Warning
    
def fetch_data(pdb_file, gene):

        #metto df in session state dfu (df unfiltered) così da non usare cvfetch ogni volta
        if st.session_state.dfu is None:
            st.session_state.dfu = cvfetch(gene, pdb_file)

        #restituisco warning
        if type(st.session_state.dfu) == Warning:
            st.error(st.session_state.dfu)

        #errore imprevisto se non è dataframe
        elif type(st.session_state.dfu) != pd.DataFrame and type(st.session_state.dfu) != Warning:
            st.error("**Warning**: An error occurred while fetching results from clinvar and setting up the dataframe.")
        
        #altrimenti stampa dataframe (in session state)
        else:
            st.info("The following table shows the mutations associated with the selected gene. The mutations are taken from the NCBI ClinVar database. You can download the dataframe and use the magnifying glass icon on the upper right corner to search for specific strings. It is also possible to sort columns by clicking on the column header.", icon='ℹ️')
            st.dataframe(st.session_state.dfu, hide_index=True)

def df_filter(pdb_file):
    with st.form('Filtro'):
        #titolo form
        st.markdown('<h5>Filter and viewer settings</h4>', unsafe_allow_html=True)

        #filtri Bfactor e trait
        bfactor = st.slider('Bfactor', min_value=0.0, max_value=100.0, value=(0.0, 100.0), step=0.1, help="In experimental structures, the B-factor column of the PDB file describes the displacement of the atomic positions from an average (mean) value (mean-square displacement). In computational predictions, the B-factor is replaced by pLDDT a measure of the confidence in the prediction. In experimental structures, the lower the B-factor the higher the quality, while in computational structures, the higher the pLDDT the higher the quality.")
        tratti_unici = list(set([item for sublist in st.session_state.dfu['Associated trait'].tolist() for item in sublist]))
        trait = st.multiselect('Associated trait', tratti_unici, help='Here, you can filter by the traits associated with the mutations in the dataset.')
        mutpos = st.number_input("Mutation position", min_value=0, max_value=len(get_residues(pdb_file)), value=0, step=1, help='Here, you can filter by the position of the mutation in the protein.')

        #stili
        c1, c2 = st.columns(2)
        with c1:
            colorepro = st.color_picker('Protein color', value='#ffffff', help='Here, you can select the color of the wild-type residues.')
            stilepro = st.selectbox('Protein style', ['cartoon', 'line', 'stick', 'sphere'], help='Here, you can select the style of the wild-type residues.')
        with c2:
            colorevar = st.color_picker('Mutated residues color', value='#ffffff', help='Here, you can select the color of the mutated residues.')
            stilevar = st.selectbox('Mutated residues styles', ['cartoon', 'line', 'stick', 'sphere'], help='Here, you can select the style of the mutated residues.')


        #colori default, color by Bfactor ed etichette
        colordef = st.toggle('Use default colors', help='Select this to use the default color scheme.', value=False)
        colorbyb = st.toggle('Color by B-factor', help='Select this to color the chain by B-factor. Lower B-factor values are represented in red, while higher B-factor values are represented in green.')
        visetich = st.toggle('Show labels', help='Select this to show the labels of the mutated residues.')

        #zoom su residuo
        zoomsuresiduo = st.number_input("Zoom to residue number:", min_value=0, max_value=len(get_residues(pdb_file)), value=0, step=1, help='Here, you can zoom to a specific residue in the structure.')
        
        #submit button (ogni volta che si submitta df filtered torna a None nel session state)
        submitted = st.form_submit_button('Apply filters and show structure', on_click=sessionstatedff)

        return bfactor, trait, colorepro, stilepro, colorevar, colordef, stilevar, colorbyb, visetich, zoomsuresiduo, mutpos, submitted
    
def init_structure(pdb_file):
    view = p3d.view(width=800, height=800)
    view.addModel(open(pdb_file, 'r').read(), 'pdb')
    view.setBackgroundColor('white')
    return view

def color_structure(view, colorepro, colorevar, colordef, colorbyb, stilevar, stilepro, pdb_file, residues):
    if colordef is True:
        if colorbyb or colorevar != '#ffffff' or colorepro != '#ffffff':
            st.warning("**Warning**: You've selected both 'Use default colors' and a custom color. The default color will be used.")
        view.setStyle({'chain': 'A'}, {stilepro: {'color': 'spectrum'}})
        return view
    elif colorevar != '#ffffff' or colorepro != '#ffffff':
        if colorbyb:
            st.warning("**Warning**: You've selected both 'Color by B-factor' and a custom color. 'Color by B-factor' will be ignored.")
        for residue in residues:
            if residue.get_id()[1] in st.session_state.dff['POS'].tolist():
            #qui si potrebbe usare view.addcolor ma vabbe
                view.setStyle({'resi': residue.get_id()[1]}, {stilevar: {'color': colorevar}})
            else:
                view.setStyle({'resi': residue.get_id()[1]}, {stilepro: {'color': colorepro}})
        return view
    elif colorbyb:
        residues = get_residues(pdb_file)
        for residue in residues:
            atoms = residue.get_atoms() 
            for atom in atoms:
                if atom.get_name() == 'CA':
                    bfactor = atom.get_bfactor()
                    normbfact = bfactor/100
                    r=(1-(normbfact))*255
                    g=(normbfact)*255
                    b=0
                    color = 'rgb('+str(r)+','+str(g)+','+str(b)+')'
                    view.setStyle({'resi': residue.get_id()[1]}, {stilepro: {'color': color}})
        return view
    else:
        st.warning("**Warning**: No color selected. The rainbow spectrum will be used for visualization.")
        view.setStyle({'chain': 'A'}, {stilepro: {'color': 'spectrum'}})
    return view

def style_structure(view, stilepro, stilevar, residues):
    for residue in residues:
        if residue.get_id()[1] in st.session_state.dff['POS'].tolist():
            view.addStyle({'resi': residue.get_id()[1]}, {stilevar: {}})
        else:
            view.addStyle({'resi': residue.get_id()[1]}, {stilepro: {}})
    return view

def add_labels(view, df, residues):
    # aggiunge stile e label per ogni mutazione
    coord_list = st.session_state.dff['Coordinates'].tolist()
    posizioni = st.session_state.dff['POS'].tolist()
    i=0
    for coord in coord_list:
        x, y, z = coord
        label = ''.join([st.session_state.dff['WT AA'].tolist()[i], str(posizioni[i])])
        view.addLabel(label, {'position': {'x': float(x), 'y': float(y), 'z': float(z)}})
        i+=1
    return view

def main():
    st.title('3DVarPro')
    st.markdown('''
                
                3DVarPro is a web app that lets you visualize clinical trait-associated point mutations on the 3D structure of the wild-type protein. 
                This app was developed by me, Attilio Turco, during my internship at [CASSMedChem](http://www.cassmedchem.unito.it). 
                A public release of the code is available [here](https://github.com/atturk/3DVarPro).

                Changelog:
                - version 1.1: Added "Show example" feature
                - version 1.2: Added filter by mutation position 
                - version 1.2.1:  You can now search using both the Gene Name and the UniprotKB AC-ID. 
                **It is also possible to search using the identifier only** (beta feature)
                
                If you encounter any issue, please contact me at attilio.turco@edu.unito.it

                Feel free to contact me if you'd like to suggest any feature!
                
                ''', unsafe_allow_html=True)

    if 'dfu' not in st.session_state:
        st.session_state.dfu = None

    if 'fileup' not in st.session_state:
        st.session_state.fileup = False

    if "disabled" not in st.session_state:
        st.session_state.disabled  = False

    if 'idonly' not in st.session_state:
        st.session_state.idonly = False

    if 'ex' not in st.session_state:
        st.session_state.ex = False


    id = ""
    pdb_file = None

    gene = ""
    uniprotid = ""

    if pdb_file is None and id == '':
        st.button('Show example', on_click=showexample, disabled=st.session_state.disabled)
        if st.session_state.ex:
            pdb_file = 'AF-Q96Q42-F1-model_v4.pdb'
            id = 'ALS2'
            st.markdown('PDB file: [AF-Q96Q42-F1-model_v4.pdb](https://alphafold.ebi.ac.uk/entry/Q96Q42) <br> Gene: [ALS2](https://www.ncbi.nlm.nih.gov/clinvar/?term=ALS2%5Bgene%5D+AND+%22mol+cons+missense%22%5Bfilter%5D)', unsafe_allow_html=True)

    if st.session_state.disabled == False:
        st.warning("**Warning**: Currently, 3DVarPro does not support PDB files containing more than one model, more than one chain, or discontinuous chains. Please ensure your input file meets these criteria for optimal visualization.")
    else:
        st.warning("**Warning**: in order to avoid bugs, you can't run the example or change the Gene/PDB file after running a search (or the example itself). Please refresh the page if you want to perform another search.", icon="⚠️")

    #carica id
    if id == "":
        id = get_id()
    
    #carica pdb
    if pdb_file is None:
        pdb_file = get_pdb()

    if id != "":
        unid_and_gene = id_to_genename_and_uniprotid(id)
        if type(unid_and_gene) == tuple:
            uniprotid = unid_and_gene[0]
            gene = unid_and_gene[1]
        else:
            st.warning("No results were found for the given ID at https://www.uniprot.org/id-mapping. Make sure the ID is correct.")

    if not st.session_state.fileup and not st.session_state.ex:
        st.button("Or Search using the identifier only (beta)", on_click=idonlysearch)
        if st.session_state.idonly:
            if id != '':
                afpdb = get_af_pdb(uniprotid)
                if type(afpdb) == str:
                    pdb_file = afpdb
                    afurl = f'https://alphafold.ebi.ac.uk/entry/{uniprotid}'
                    st.info(f'The following Alphafold model was used: {afurl}', icon='ℹ️')
                else:
                    st.warning("Error while retrieving the pdb model. Try uploading the PDB file.")
                    st.stop()
            else: 
                st.warning("Please insert the identifier first.")

    #dati da clinvar
    if pdb_file is not None and gene != '':
        if st.session_state.disabled == False:
            disable()
            st.rerun()
        fetch_data(pdb_file, gene)

        
    submitted = False
    #form per filtrare df
    if pdb_file is not None and gene != '' and type(st.session_state.dfu) == pd.DataFrame:

        if 'dff' not in st.session_state:
            st.session_state.dff = None

        bfactor, trait, colorepro, stilepro, colorevar, colordef, stilevar, colorbyb, visetich, zoomsuresiduo, mutpos, submitted = df_filter(pdb_file)

    if submitted:

        #filtra per bfactor
        st.session_state.dff = st.session_state.dfu[(st.session_state.dfu['Bfactor'] >= bfactor[0]) & (st.session_state.dfu['Bfactor'] <= bfactor[1])]
        #filtra per trait
        if trait != []:
            st.session_state.dff = st.session_state.dff[st.session_state.dfu['Associated trait'].apply(lambda lista: any(item in lista for item in trait))]
        st.info('Below is the filtered dataframe and the generated structure.', icon='ℹ️')
        st.dataframe(st.session_state.dff, hide_index=True)
        #filtra per posizione
        if mutpos != 0:
            st.session_state.dff = st.session_state.dff[st.session_state.dff['POS'] == mutpos]


        #VARIANTI NELLA STRUTTURA
        #secondo filtro dataframe per tenere solo AA che effettivamente esistono nella proteina
        if type(check_df(st.session_state.dff, pdb_file)) == Warning:
            st.warning("**Error**: The filtered dataset is empty. Please check your filter settings.")
            st.stop()
        else:
            st.session_state.dff = check_df(st.session_state.dff, pdb_file)

        #lista residui che servirà dopo
        residues = get_residues(pdb_file)

        #crea visualizzatore struttura
        st.markdown('<h4>3. Structure viewer</h4>', unsafe_allow_html=True)

        view = init_structure(pdb_file)
        
        #colore
        view = color_structure(view, colorepro, colorevar, colordef, colorbyb, stilevar, stilepro, pdb_file, residues)

        #stile
        if stilepro != stilevar:
            view = style_structure(view, stilepro, stilevar, residues)
        #etichette
        if visetich:
            view = add_labels(view, st.session_state.dff, residues)

        #zoom su residuo
        if zoomsuresiduo != 0:
            view.zoomTo({'resi': zoomsuresiduo})
        else:
            view.zoomTo({'chain': 'A'})

        #mostra struttura
        showmol(view, width=800, height=800)


if __name__ == '__main__':
    main()