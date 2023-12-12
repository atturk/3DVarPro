from utils import *
import streamlit as st
import pandas as pd
import Bio.PDB as BP
import py3Dmol as p3d
from stmol import *
from Bio.SeqUtils import IUPACData as aa

def showexample():
    st.session_state.ex = True
    
def sessionstatedfu():
    st.session_state.dfu = None

def sessionstatedff():
    st.session_state.dff = None

def get_pdb():
        st.markdown('<h4>Upload the PDB file of the wild-type protein.</h4>', unsafe_allow_html=True)
        uploaded_file = st.file_uploader(label='pdb uploader',type=['pdb'], label_visibility='collapsed', on_change=sessionstatedfu)
        if uploaded_file is not None:
            pdb_file = create_pdb(uploaded_file)
            st.write('File uploaded successfully')
            return pdb_file

def get_gene():
    st.markdown('<h4>Insert the gene name</h4>', unsafe_allow_html=True)
    gene = st.text_input('gene', label_visibility='collapsed', on_change=sessionstatedfu).upper()
    return gene

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
    st.markdown('3DVarPro is a web app that lets you visualize clinical trait-associated point mutations on the 3D structure of the wild-type protein. This app was developed by me, Attilio Turco, during my internship at [CASSMedChem](http://www.cassmedchem.unito.it). A public release of the code is available [here](https://github.com/atturk/3DVarPro).')
    st.warning("**Warning**: Currently, 3DVarPro does not support PDB files containing more than one model, more than one chain, or discontinuous chains. Please ensure your input file meets these criteria for optimal visualization.")

    if 'dfu' not in st.session_state:
        st.session_state.dfu = None

    gene = None
    pdb_file = None

    #carica pdb
    pdb_file = get_pdb()

    #carica gene
    gene = get_gene()

    if 'ex' not in st.session_state:
        st.session_state.ex = False
    
    if pdb_file is None and gene == '':
        st.button('Show example', on_click=showexample)
        if st.session_state.ex:
            pdb_file = 'AF-Q96Q42-F1-model_v4.pdb'
            gene = 'ALS2'
            st.markdown('PDB file: [AF-Q96Q42-F1-model_v4.pdb](https://alphafold.ebi.ac.uk/entry/Q96Q42) <br> Gene: [ALS2](https://www.ncbi.nlm.nih.gov/clinvar/?term=ALS2%5Bgene%5D+AND+%22mol+cons+missense%22%5Bfilter%5D)', unsafe_allow_html=True)

    #dati da clinvar
    if pdb_file is not None and gene != '':
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
        #filtra per posizione
        if mutpos != 0:
            st.session_state.dff = st.session_state.dff[st.session_state.dff['POS'] == mutpos]
        st.info('Below is the filtered dataframe and the generated structure.', icon='ℹ️')
        st.dataframe(st.session_state.dff, hide_index=True)
        print(colordef)


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
