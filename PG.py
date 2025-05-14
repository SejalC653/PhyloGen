import streamlit as st
import subprocess
import os
from Bio import SeqIO, Phylo, AlignIO
import matplotlib.pyplot as plt
from io import StringIO
import tempfile
import pandas as pd
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, ParsimonyTreeConstructor, ParsimonyScorer, NNITreeSearcher
import base64
import requests
import time

# ✅ First Streamlit Command
st.set_page_config(
    page_title="PhyloGen - Home",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Custom CSS --- 
st.markdown("""
    <style>
        .content-box {
            padding: 2.5rem;
            border-radius: 12px;
            margin-top: 20px;
            box-shadow: 0px 0px 10px rgba(0,0,0,0.05);
        }

        .features-box {
            background-color: #cce4f6; /* Light blue background for Features section */
            padding: 2.5rem;
            border-radius: 12px;
            box-shadow: 0px 0px 10px rgba(0,0,0,0.05);
        }

        h1 {
            color: #114b5f;
            margin-bottom: 0 !important;  /* Removed space below title */
        }

        .stTabs [data-baseweb="tab-list"] {
            gap: 20px;
            margin-top: 0 !important;  /* Removed space above tabs */
        }

        .stTabs [data-baseweb="tab"] {
            font-size: 18px;
            font-weight: 600;
            padding: 12px 24px;
            border-radius: 10px 10px 0 0;
            background-color: #cce4f6;
            color: #114b5f;
            transition: background-color 0.3s ease;
            margin-top: 0;  /* No margin to prevent any gap */
        }

        .stTabs [aria-selected="true"] {
            background-color: #114b5f;
            color: white;
        }

    </style>
""", unsafe_allow_html=True)

# --- Tabs ---



with st.sidebar:
    st.title("PhyloGen V1.0")
    st.header("Navigation")
    page = st.sidebar.radio("Go to", ["🏠️ Home", "📘 User Guide", "🧪 PhyloGen Tool", "ℹ️ About"])
    if page == "🏠️ Home":
        st.divider()
        st.header("Workflow")
        st.markdown("""
        **FASTA to Phylogenetic Tree Tool**

        - Upload FASTA sequences  
        - Perform multiple sequence alignment  
        - Build phylogenetic trees  
        - Export results
        """)
        st.divider()

    elif page == "🧪 PhyloGen Tool":
        st.header("Workflow")
        st.markdown("""
        **FASTA to Phylogenetic Tree Tool**

        - Upload FASTA sequences  
        - Perform multiple sequence alignment  
        - Build phylogenetic trees  
        - Export results
        """)
        st.divider()
        st.markdown("Choose a sequence type")
        sequence_type = st.radio("Select Sequence Type", ["dna", "protein"])
        st.divider()



if page == "🏠️ Home":
    st.markdown('<div class="content-box">', unsafe_allow_html=True)

    st.markdown("## 🧬 Welcome to PhyloGen")
    st.markdown(""" 
    **PhyloGen** is an intuitive web application for building phylogenetic trees from DNA or protein sequences in FASTA format.  
    Designed for students, researchers, and life scientists, PhyloGen streamlines alignment and tree construction using trusted tools—without needing to write code or use the command line.
    """)

    # Removed the light blue box – replaced with regular content-box
    st.markdown("### 🛠️ Features")
    st.markdown("""  
    - 📁 **FASTA File Upload** – Upload .fasta or .fa files with 3+ sequences  
    - 🧬 **Multiple Sequence Alignment** – Choose from **MAFFT**, **MUSCLE**, or **ClustalW**  
    - 🌳 **Tree Construction Methods**:
        - Distance-based: **Neighbor Joining**, **UPGMA**
        - Character-based: **Maximum Parsimony**, **Maximum Likelihood (FastTree)**
    - ⚙️ **Customizable Parameters** – Gap penalties, substitution models for DNA or protein  
    - 🔢 **Substitution Model** - Choose a substitution model for construction of matrix
        - DNA: **Blossom62**, **Blossom45**, **Blossom50**, **Blastp**, **Dayhoff**, **Pam250**, 
               **Pam70**, **Pam30**, **Identity**  
        - Protein: **Megablast**, **Blastn**, **Trans**, **Identity**
    - 📊 **Interactive Preview** – See uploaded sequence details before processing  
    - 📤 **Export Results**:
        - Aligned sequences (FASTA)  
        - Tree files (Newick format)  
        - Tree images (PNG, ready for publication)
    """)

    st.markdown("</div>", unsafe_allow_html=True)


def get_image_base64(image_path):
    with open(image_path, "rb") as image_file:
        encoded = base64.b64encode(image_file.read()).decode()
    return f"data:image/png;base64,{encoded}"

if page == "📘 User Guide":
    st.markdown('<div class="content-box">', unsafe_allow_html=True)

    st.markdown("### 📘 User Guide")

    def image_border(img_path, caption):
        img_base64 = get_image_base64(img_path)
        st.markdown(
            f'''
            <div style="border: 2px solid black; display: inline-block; padding: 5px; margin-bottom: 10px;">
                <img src="{img_base64}" width="500" />
                <p style="text-align: center;">{caption}</p>
            </div>
            ''',
            unsafe_allow_html=True
        )

    st.markdown("1. **Upload FASTA File** – Ensure the file has 3 or more sequences.")
    image_border("1choosefile.png", "Choose a multi-fasta file")

    st.markdown("2. **Preview Sequences** - Sequences IDs and Description can be displayed.")
    image_border("2prevseq.png", "Display sequence information in table")

    st.markdown("3. **Select Sequence type** - DNA or Protein.")
    image_border("3seqtype.png", "Choose the type of data")

    st.markdown("""
    4. **Select Alignment Tool**  
    - Choose from **MAFFT**, **MUSCLE**, **ClustalW** or **TCoffee**  
    - Set gap penalties and substitution models  
    - Choose tree construction method:  
      ~ Distance-based: **Neighbor Joining**, **UPGMA**  
      ~ Character-based: **Maximum Parsimony**, **Maximum Likelihood**  
    """)
    image_border("4alnsettng.png", "Set alignment options and tree method")

    st.markdown("5. **View & Download** – Preview and download results.")
    image_border("5restree.png", "Download phylogenetic tree as PNG")
    image_border("6resaln.png", "Download aligned sequences as FASTA")
    image_border("7resdata.png", "Download tree data as Newick format")

    st.markdown("</div>", unsafe_allow_html=True)

    # Path to the existing PDF file
    pdf_path = "PGManual.pdf"  # Replace with your actual PDF file path

    # Provide a download button to download the PDF file
    with open(pdf_path, "rb") as pdf_file:
        pdf_bytes = pdf_file.read()

    st.download_button(
        label="Download User Manual",  # Button label
        data=pdf_bytes,                   # PDF data in binary
        file_name="PGManual.pdf",      # Filename to be given to the downloaded PDF
        mime="application/pdf"            # MIME type for PDF
      )

   



if page == "🧪 PhyloGen Tool":
    st.markdown("""
    <style>
        .main {background-color: #bd73eb;}
        .stButton>button {border-radius: 5px;}
        .stSelectbox, .stSlider {padding: 5px;}
        .st-eb {background-color: #e8f4f8;}
        h1 {color: #2a5a72;}
        h2 {color: #3d7b96;}
    </style>
    """, unsafe_allow_html=True)
    
    # Functions
    def validate_fasta(uploaded_file):
        try:
            with tempfile.NamedTemporaryFile(delete=False) as tmp:
                tmp.write(uploaded_file.getvalue())
                tmp.seek(0)
                records = list(SeqIO.parse(tmp.name, "fasta"))
                if len(records) < 3:
                    return False, "At least 3 sequences required"
                return True, ""
        except Exception as e:
            return False, f"Invalid FASTA format: {str(e)}"

    def run_ml(input_file, output_file):
        cmd = "FastTree -gtr -nt alignment.fasta > tree.nwk"
        subprocess.run(cmd, shell=True, check=True)

    # Main App
    st.title("PhyloGen Tool")
    st.markdown("Upload your FASTA file and generate publication-quality phylogenetic trees!")

    # File Upload Section
    uploaded_file = st.file_uploader("Choose a file", type=["fasta", "fa"])

    if uploaded_file is not None:
        st.success(f"Uploaded file: {uploaded_file.name}")
        file_content = uploaded_file.read().decode("utf-8")
        sequences = list(SeqIO.parse(StringIO(file_content), "fasta"))

        if st.checkbox("Preview Sequences"):
            seq_data = []
            for seq in sequences:
                seq_data.append({
                    "ID": seq.id,
                    "Description": seq.description,
                    "Length": len(seq.seq),
                    "Type": "dna" if "U" not in seq.seq else "protein"
                })
            st.dataframe(pd.DataFrame(seq_data), height=200)


    # Processing Section
    email = "your@email.com"
    is_valid, error_msg = validate_fasta(uploaded_file)
    if uploaded_file and is_valid:
        with st.expander("⚙️ STEP 2: Processing Options", expanded=True):
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("Alignment Settings")
                align_tool = st.selectbox("Alignment Algorithm", ["mafft", "muscle", "clustalo", "tcoffee"])
                gap_open = st.number_input("Gap Open Penalty (max=3.0)", value=2.7)
                gap_extend = st.number_input("Gap Extension Penalty (max=1.0)", value=1.8)
            with col2:
                st.subheader("Tree Construction")
                tree_method = st.selectbox("Construction Method", ["Distance based", "Character based"])
                if tree_method == "Distance based":
                    dmat = st.selectbox("Distance based", ["Neighbour Joining", "UPGMA"])
                else:
                    dmat = st.selectbox("Character based", ["Maximum Parsimony", "Maximum Likelihood"])
                substitution_model = st.selectbox(
                    "Substitution model",
                    ["blosum62", "blosum45", "blosum50", "blastp", "dayhoff", "pam250", "pam70", "pam30", "identity"]
                    if sequence_type == "protein" else
                    ["megablast", "blastn", "trans", "identity"]
                )

        if st.button("Run Full Analysis", type="primary"):
            with st.status("Running Analysis...", expanded=True) as status:
                temp_fasta = "temp.fasta"
                temp_aln = "alignment.fasta"
                tree = "output_file"
                with open(temp_fasta, "w") as f:
                    f.write(file_content)
                    with st.spinner("Submitting job..."):

                        # Submit job to the selected tool
                        run_url = f"https://www.ebi.ac.uk/Tools/services/rest/{align_tool}/run"
                        data = {
                            "email": email,
                            "sequence": file_content,
                            "stype": sequence_type,
                            "gapopen": str(gap_open),
                            "gapextend": str(gap_extend)
                        }

                        try:
                            response = requests.post(run_url, data=data)
                            response.raise_for_status()
                            job_id = response.text.strip()
                            st.success(f"Job submitted! ID: {job_id}")
                        except Exception as e:
                            st.error(f"Submission failed: {e}")
                            st.stop()

                    # Poll for completion
                    status_url = f"https://www.ebi.ac.uk/Tools/services/rest/{align_tool}/status/{job_id}"
                    start_time = time.time()
                    timeout = 600  # 10 minutes

                    with st.spinner("Waiting for job to finish..."):
                        while True:
                            job_status = requests.get(status_url).text.strip()
                            if job_status == "FINISHED":
                                break
                            elif job_status in ["ERROR", "FAILURE", "NOT_FOUND"]:
                                st.error(f"Job failed with status: {job_status}")
                                st.stop()
                            elif time.time() - start_time > timeout:
                                st.error("Job timed out.")
                                st.stop()
                            time.sleep(5)
                            st.write(f"Status: {job_status}...")

                    # Fetch the result
                    result_url = f"https://www.ebi.ac.uk/Tools/services/rest/{align_tool}/result/{job_id}/fa"
                    result_response = requests.get(result_url)
                    if result_response.status_code == 200:
                        temp_aln = result_response.text
                st.write("Building phylogenetic tree...")
                
                def distmat(alignment_file):
                    calculator = DistanceCalculator(substitution_model)
                    return calculator.get_distance(alignment)
                with open("alignment.fasta", "w") as f:
                    f.write(temp_aln)

                if tree_method == "Distance based":
                    alignment = AlignIO.read("alignment.fasta", "fasta")
                    constructor = DistanceTreeConstructor()
                    tree = constructor.nj(distmat(alignment)) if dmat == "Neighbour Joining" else constructor.upgma(distmat(temp_aln))
                    Phylo.write(tree, "tree.nwk", "newick")

                elif tree_method == "Character based":
                    alignment = AlignIO.read("alignment.fasta", "fasta")
                    if dmat == "Maximum Parsimony":
                        calculator = DistanceCalculator('identity')
                        distance_matrix = calculator.get_distance(alignment)
                        constructor = DistanceTreeConstructor()
                        starting_tree = constructor.nj(distance_matrix)
                        scorer = ParsimonyScorer()
                        searcher = NNITreeSearcher(scorer)
                        pars_constructor = ParsimonyTreeConstructor(searcher, starting_tree)
                        tree = pars_constructor.build_tree(alignment)
                        Phylo.write(tree, "tree.nwk", "newick")
                    else:
                        run_ml(temp_aln, tree)
                        tree = Phylo.read("tree.nwk", "newick")

                num_taxa = len(tree.get_terminals())
                fig = plt.figure(figsize=(8, max(4, num_taxa * 0.4)))
                ax = fig.add_subplot(1, 1, 1)
                Phylo.draw(tree, axes=ax, label_func=lambda clade: clade.name if clade.is_terminal() else None,
                           branch_labels=lambda clade: f"{round(clade.branch_length, 3)}" if clade.branch_length else '',
                           do_show=False)
                ax.set_xlabel("Branch Length")
                ax.set_ylabel("Taxa")
                ax.set_title(f"{dmat}", fontsize=14)
                plt.tight_layout()
                tree_image_path = "tree.png"
                plt.savefig(tree_image_path)
                plt.close()

                status.update(label="Analysis Complete!", state="complete", expanded=False)

            # Results
            st.subheader("Results")
            t1, t2, t3 = st.tabs(["Tree Visualization", "Alignment", "Tree Data"])

            with t1:
                if os.path.exists(tree_image_path):
                    with open(tree_image_path, "rb") as img_file:
                        st.download_button("Download Tree Image", data=img_file, file_name="phylogenetic_tree.png", mime="image/png")
                    st.pyplot(fig)
                else:
                    st.error("Error: The phylogenetic tree image could not be generated.")

            with t2:
                with open("alignment.fasta") as f:
                    st.download_button("Download Alignment", f.read(), "alignment.fasta")
                st.code(open("alignment.fasta").read(), language='fasta')

            with t3:
                tree_nwk_path = "tree.nwk"
                if os.path.exists(tree_nwk_path):
                    with open(tree_nwk_path, "r") as tree_file:
                        st.download_button("Download Newick", tree_file.read(), "tree.nwk")
                    st.code(open(tree_nwk_path).read(), language='newick')

    # Cleanup
    if uploaded_file:
        for f in ["temp.fasta", "alignment.fasta", "tree.txt", "tree.nwk"]:
            if os.path.exists(f):
                os.remove(f)
if page == "ℹ️ About":
    st.title("📌 About The Developer")

    st.markdown("""
    ---
     👨‍💻 Author 
    ---
    Sejal Chaudhari
    ---
    Master's Student in Bioinformatics  
    DES Pune University, Maharashtra, India
   
    ---
    ### 🎯 Purpose
    - To simplify sequence alignment and tree building for educational and academic purposes.
    - To provide an interactive visualization platform for understanding evolutionary relationships.
    - To demonstrate integration of core bioinformatics algorithms in a Python web framework.

    ---
    ### 🌟 Key Features
    - 🧬 FASTA Upload & Sequence Preview  
    - 🧰 Multiple Sequence Alignment (MAFFT, MUSCLE, ClustalW)  
    - 🌳 Phylogenetic Tree Construction (NJ, UPGMA, Maximum Parsimony, ML)  
    - 📐 Adjustable Gap Penalties & Substitution Models  
    - 📤 Downloadable Outputs (Newick, FASTA, Tree Image)  
    - 📊 Tree Visualization with Branch Lengths  
    - 🧠 DNA and Protein Modes  

    ---
    ### 🛠️ Tools & Technologies Used
    - Python, Streamlit for Web Interface  
    - Biopython for Sequence Parsing, Alignment, and Phylogenetics  
    - Matplotlib for Tree Rendering  
    - Subprocess calls for CLI-based tools like MAFFT, MUSCLE, and FastTree  

    ---
    ### 💡 Benefits  
    - Beginner-friendly layout with expandable configuration options  
    - Clear visualization and downloadable results  
    - Modular architecture — easy to enhance with new tools  
    - Easy, usable interface
    - Publication quality trees just a click away!

    ---
    ### 🚀 Planned Enhancements
    - Enhanced tree editing tools (advanced visulaisation features)
    - Session saving and project export  
    - Support for multi-gene and whole-genome datasets  
    - Integration with online databases (e.g., GenBank fetch)  

    ---
    ### 👨‍🏫 Mentorship & Acknowledgements
    **Special thanks to Dr. Kushagra Kashyap**  
    Assistant Professor (Bioinformatics), Department of Life Sciences  
    DES Pune University, for his academic guidance and mentorship that shaped the design and development of this tool.
    🔗 [Connect on LinkedIn](https://www.linkedin.com/in/dr-kushagra-kashyap-b230a3bb)

    ---
    ### 🔗 Connect & Contact
    📧 Email: [sejalc2468@gmail.com](mailto:sejalc2468@gmail.com)  

    *Thank you for exploring PhyloGen Tool! Your feedback is always welcome.*
    """)
