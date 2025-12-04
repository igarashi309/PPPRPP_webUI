import os
from pathlib import Path
from app import run_unified_analysis

# User provided sequences
PROTEIN_SEQ = """>PPR791(Rf1a)_IR24
MARRAASRAVGALRSDGSIQGRGGRAGGSGAEDARHVFDELLRRGRGASIYGLNRALADVARDSPAAAVSRYNRMARAGADEVTPDLCTYGILIGCCCRAGRLDLGFAALGNVIKKGFRVDAIAFTPLLKGLCADKRTSDAMDIVLRRMTELGCIPNVFSYNILLKGLCDENRSQEALELLHMMADDRGGGSPPDVVSYTTVINGFFKEGDSDKAYSTYHEMLDRGILPDVVTYNSIIAALCKAQAMDKAMEVLNTMVKNGVMPDCMTYNSILHGYCSSGQPKEAIGFLKKMRSDGVEPDVVTYSLLMDYLCKNGRCMEARKIFDSMTKRGLKPEITTYGTLLQGYATKGALVEMHGLLDLMVRNGIHPDHYVFSILICAYAKQGKVDQAMLVFSKMRQQGLNPNAVTYGAVIGILCKSGRVEDAMLYFEQMIDEGLSPGNIVYNSLIHGLCTCNKWERAEELILEMLDRGICLNTIFFNSIIDSHCKEGRVIESEKLFELMVRIGVKPNVITYNTLINGYCLAGKMDEAMKLLSGMVSVGLKPNTVTYSTLINGYCKISRMEDALVLFKEMESSGVSPDIITYNIILQGLFQTRRTAAAKELYVRITESGTQIELSTYNIILHGLCKNKLTDDALQMFQNLCLMDLKLEARTFNIMIDALLKVGRNDEAKDLFVAFSSNGLVPNYWTYRLMAENIIGQGLLEELDQLFLSMEDNGCTVDSGMLNFIVRELLQRGEITRAGTYLSMIDEKHFSLEASTASLFIDLLSGGKYQEYYRFLPEKYKSFIESLSC"""

RNA_SEQ = """>atp6-orf79_BT
ggcgttggtttttccaacgaaaaacgaattcgaacatgaaatgaaaaaacaaatgctagttgtcattttttgtaacgcaatacacttccgcgtttgtcccgccccctaaaatcttcggctcctcgtttttattcaattatgaaattactcattttgatggagatttatagcatcattcaag"""

def test_run():
    print("Starting Test Run...")
    
    # Call the function directly
    # inputs=[prot_text, prot_file, rna_text, rna_file, tool_select, analysis_select, job_name_input]
    
    result = run_unified_analysis(
        protein_text=PROTEIN_SEQ,
        protein_file=None,
        rna_text=RNA_SEQ,
        rna_file=None,
        tool_choice="PPRFinder",
        analysis_type="kobayashi_gap0_ps",
        job_name="test_run_rf1a"
    )
    
    # outputs=[status_output, df_output, csv_output, plot_browser, plot_detail, df_state, motif_table, motif_plot]
    status, df, csv, fig_browser, fig_detail, df_state, motif_df, motif_fig = result
    
    print(f"Status: {status}")
    if "Error" in status:
        print("TEST FAILED")
    else:
        print("TEST PASSED")
        print(f"CSV Path: {csv}")
        print(f"Motifs Found: {len(motif_df)}")
        print(f"Analysis DF Shape: {df.shape}")
        
        if not df.empty:
            print("Testing show_manual_detail...")
            # Pick a position from the DF
            test_pos = df.iloc[0]['Genome_Position']
            print(f"Testing position: {test_pos}")
            
            from app import show_manual_detail
            # Test Forward
            fig_fwd, seq_fwd = show_manual_detail(test_pos, "Forward", df)
            print(f"Forward Figure Title: {fig_fwd.layout.title.text}")
            print(f"Forward Sequence: {seq_fwd}")
            print(f"Number of traces: {len(fig_fwd.data)}")
            if len(fig_fwd.data) > 1:
                print(f"Second trace name: {fig_fwd.data[1].name}")
            
            # Test Reverse (might fail if no site there, but should run)
            fig_rev, seq_rev = show_manual_detail(test_pos, "Reverse", df)
            print(f"Reverse Figure Title: {fig_rev.layout.title.text}")
            
            # Test NCBI Fetch and Analysis with AP012527
            print("\nTesting NCBI Fetch and Analysis with AP012527...")
            from ppr_analysis import ncbi
            from app import run_ncbi_analysis # Import the new function wrapper if needed, or just use run_unified_analysis directly?
            # Actually test_run calls run_unified_analysis directly for the first part, 
            # but for the NCBI part we should probably call run_ncbi_analysis to test the full flow including file generation?
            # But run_ncbi_analysis is in app.py.
            from app import run_ncbi_analysis
            
            test_accession = "AP012527"
            fasta_path = ncbi.fetch_fasta(test_accession, output_dir=".")
            print(f"Fetched FASTA: {fasta_path}")
            
            if os.path.exists(fasta_path):
                print("FASTA file exists.")
                
                # Run analysis with fetched FASTA
                with open(fasta_path, "r") as f:
                    rna_text_ncbi = f.read()
                    
                print("Running analysis on fetched sequence...")
                # Call run_ncbi_analysis to test the full flow
                # inputs=[ncbi_accession, prot_text, prot_file, tool_choice, analysis_type, job_name]
                # Test PS_Scan API
                print("Testing with PS_Scan API...")
                generator = run_ncbi_analysis(
                    accession_id=test_accession,
                    protein_text=PROTEIN_SEQ,
                    protein_file=None,
                    tool_choice="PS_Scan", # Changed to PS_Scan
                    analysis_type="kobayashi_gap0_ps",
                    job_name="test_ncbi_ap012527_ps_scan"
                )
                
                # Iterate to get the final result
                result_ncbi = None
                print("Iterating generator...")
                for update in generator:
                    result_ncbi = update
                    # print(f"Update: {update[0]}") # Print status message
                
                if result_ncbi:
                    status_n, fasta_path_n, bedgraph_path_n, gff_path_n, df_n, csv_n, fig_browser_n, log_n = result_ncbi
                    print(f"NCBI Analysis Status: {status_n}")
                    print(f"NCBI Analysis DF Shape: {df_n.shape if df_n is not None else 'None'}")
                    print(f"FASTA Path for Download: {fasta_path_n}")
                    print(f"BedGraph Path for Download: {bedgraph_path_n}")
                    print(f"GFF3 Path for Download: {gff_path_n}")
                    
                    if df_n is not None and not df_n.empty:
                        if os.path.exists(fasta_path_n) and os.path.exists(bedgraph_path_n):
                            print("Downloadable files (FASTA, BedGraph) generated successfully.")
                        else:
                            print("Error: Downloadable files missing.")
                            
                        if gff_path_n and os.path.exists(gff_path_n):
                             print("GFF3 file generated successfully.")
                        else:
                             print("Warning: GFF3 file missing or not fetched.")
                    
                    # Clean up FASTA
                    if os.path.exists(fasta_path):
                        os.remove(fasta_path)
                else:
                    print("Generator yielded no results.")
            else:
                print("FASTA file NOT found.")
        else:
            print("Analysis DF is empty. Cannot test show_manual_detail.")

if __name__ == "__main__":
    test_run()
