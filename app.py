from utils import *
import streamlit as st
import pandas as pd
import Bio.PDB as BP
import py3Dmol as p3d
from stmol import *
from Bio.SeqUtils import IUPACData as aa

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

        #se non ci sono varianti warning
        if type(st.session_state.dfu) == Warning:
            st.warning(st.session_state.dfu)

        #errore imprevisto se non è dataframe
        elif type(st.session_state.dfu) != pd.DataFrame:
            st.error("**Warning**: An error occurred while fetching results from clinvar and setting up the dataframe.")
        
        #altrimenti stampa dataframe (in session state)
        else:
            st.dataframe(st.session_state.dfu, hide_index=True)

def df_filter(pdb_file):
    with st.form('Filtro'):
        #titolo form
        st.markdown('<h5>Filter and viewer settings</h4>', unsafe_allow_html=True)

        #filtri Bfactor e trait
        bfactor = st.slider('Bfactor', min_value=0.0, max_value=100.0, value=(0.0, 100.0), step=0.1)
        tratti_unici = list(set([item for sublist in st.session_state.dfu['Associated trait'].tolist() for item in sublist]))
        trait = st.multiselect('Associated trait', tratti_unici)

        #stili
        c1, c2 = st.columns(2)
        with c1:
            colorepro = st.color_picker('Protein color', value='#ffffff')
            stilepro = st.selectbox('Protein style', ['cartoon', 'line', 'stick', 'sphere'])
        with c2:
            colorevar = st.color_picker('Mutated residues color', value='#ffffff')
            stilevar = st.selectbox('Mutated residues styles', ['cartoon', 'line', 'stick', 'sphere'])

        #color by Bfactor ed etichette
        colorbyb = st.toggle('Color by Bfactor')
        visetich = st.toggle('Show labels')

        #zoom su residuo
        zoomsuresiduo = st.number_input("Zoom to residue number:", min_value=0, max_value=len(get_residues(pdb_file)), value=0, step=1)
        
        #submit button (ogni volta che si submitta df filtered torna a None nel session state)
        submitted = st.form_submit_button('Apply filters and show structure', on_click=sessionstatedff)

        return bfactor, trait, colorepro, stilepro, colorevar, stilevar, colorbyb, visetich, zoomsuresiduo, submitted
    
def init_structure(pdb_file):
    view = p3d.view(width=800, height=800)
    view.addModel(open(pdb_file, 'r').read(), 'pdb')
    view.setBackgroundColor('white')
    return view

def color_structure(view, colorepro, colorevar, colorbyb, stilevar, stilepro, pdb_file, residues):
    if colorevar != '#ffffff' or colorepro != '#ffffff':
        if colorbyb:
            st.warning("**Warning**: You've selected both 'Color by B-factor' and a custom color. 'Color by B-factor' will be ignored.")
        for residue in residues:
            if residue.get_id()[1] in st.session_state.dff['POS'].tolist():
            #qui si potrebbe usare view.addcolor ma vabbe
                view.setStyle({'resi': residue.get_id()[1]}, {stilevar: {'color': colorevar}})
            else:
                view.setStyle({'resi': residue.get_id()[1]}, {stilepro: {'color': colorepro}})
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
    st.write('3DVarPro is a web app that lets you visualize clinical trait-associated point mutations on the 3D structure of the wild-type protein.')
    st.warning("**Warning**: Currently, 3DVarPro does not support PDB files containing more than one model, more than one chain, or discontinuous chains. Please ensure your input file meets these criteria for optimal visualization.")

    if 'dfu' not in st.session_state:
        st.session_state.dfu = None

    gene = None
    pdb_file = None

    #carica pdb
    pdb_file = get_pdb()

    #carica gene
    gene = get_gene()

    #dati da clinvar
    if pdb_file is not None and gene != '':
        fetch_data(pdb_file, gene)
        
    submitted = False
    #form per filtrare df
    if pdb_file is not None and gene != '' and type(st.session_state.dfu) == pd.DataFrame:

        if 'dff' not in st.session_state:
            st.session_state.dff = None

        bfactor, trait, colorepro, stilepro, colorevar, stilevar, colorbyb, visetich, zoomsuresiduo, submitted = df_filter(pdb_file)

    if submitted:

        #filtra per bfactor
        st.session_state.dff = st.session_state.dfu[(st.session_state.dfu['Bfactor'] >= bfactor[0]) & (st.session_state.dfu['Bfactor'] <= bfactor[1])]
        #filtra per trait
        if trait != []:
            st.session_state.dff = st.session_state.dff[st.session_state.dfu['Associated trait'].apply(lambda lista: any(item in lista for item in trait))]
        st.dataframe(st.session_state.dff, hide_index=True)


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
        view = color_structure(view, colorepro, colorevar, colorbyb, stilevar, stilepro, pdb_file, residues)

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