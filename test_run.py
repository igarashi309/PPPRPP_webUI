import os
from pathlib import Path
from app import run_unified_analysis

# User provided sequences
PROTEIN_SEQ = """>PPR791(Rf1a)_IR24
MARRAASRAVGALRSDGSIQGRGGRAGGSGAEDARHVFDELLRRGRGASIYGLNRALADVARDSPAAAVSRYNRMARAGADEVTPDLCTYGILIGCCCRAGRLDLGFAALGNVIKKGFRVDAIAFTPLLKGLCADKRTSDAMDIVLRRMTELGCIPNVFSYNILLKGLCDENRSQEALELLHMMADDRGGGSPPDVVSYTTVINGFFKEGDSDKAYSTYHEMLDRGILPDVVTYNSIIAALCKAQAMDKAMEVLNTMVKNGVMPDCMTYNSILHGYCSSGQPKEAIGFLKKMRSDGVEPDVVTYSLLMDYLCKNGRCMEARKIFDSMTKRGLKPEITTYGTLLQGYATKGALVEMHGLLDLMVRNGIHPDHYVFSILICAYAKQGKVDQAMLVFSKMRQQGLNPNAVTYGAVIGILCKSGRVEDAMLYFEQMIDEGLSPGNIVYNSLIHGLCTCNKWERAEELILEMLDRGICLNTIFFNSIIDSHCKEGRVIESEKLFELMVRIGVKPNVITYNTLINGYCLAGKMDEAMKLLSGMVSVGLKPNTVTYSTLINGYCKISRMEDALVLFKEMESSGVSPDIITYNIILQGLFQTRRTAAAKELYVRITESGTQIELSTYNIILHGLCKNKLTDDALQMFQNLCLMDLKLEARTFNIMIDALLKVGRNDEAKDLFVAFSSNGLVPNYWTYRLMAENIIGQGLLEELDQLFLSMEDNGCTVDSGMLNFIVRELLQRGEITRAGTYLSMIDEKHFSLEASTASLFIDLLSGGKYQEYYRFLPEKYKSFIESLSC"""

RNA_SEQ = """>atp6-orf79_BT
ttacttaggaaagactacacgaatagaggtgccccaccacgaatagtcaacctctatagggtgcccggcttttatattatgatctatgattttttccattaaagccgggtcgtacagaggggcctgatactttagcaaagcaaacaacagggcgcctcctacaacgacaccgaagatgaaatatggaagaccgttagtccctcgggtagtggagaagagccatcggaccagatttgccatagcacaataaaataatatttttaagctccgctctccttatattggagagacttgtttggggtcggcgttggtttttccaacgaaaaacgaattcgaacatgaaatgaaaaaacaaatgctagttgtcattttttgtaacgcaatacacttccgcgtttgtcccgccccctaaaatcttcggctcctcgtttttattcaattatgaaattactcattttgatggagatttatagcatcattcaagtaaatacaaattgagatcgtagaaacatgagcttgtaatatagctacacctaattccagaccggttaatgctagaactataaataagggaccaagatctcctatgaaatagaaaatattattcagaaatagcatagtccaagcgaacccacttaaaatctttactgaactatgaccggccatcatattagcaaataaacgtattcctgagcttaatgcacgaaaacaatgagagattagctcaaggagtactaaaaaaggtgctaatggcagtgggactcccgctggtaataagaagctaaaaaaatgaagcccatgtctttgaaatccaacgatcgtaatgcctataaaaatggaaaatgaaagagccaaagtaatgagaaaatgacttgtcactgtgaagctaaagggtatcataccctggggattacgaaataacgaaaaagtaaaagtgaccgagatgcgagggaaaaacttttgtttaacatttccacctatttgttcgtttaccaggttcagcacgaaatcataaataagctctaccaaggattgccatgcatttggcactgactttcccccgccctttttcgtaacaactccaaataaagatggcaccaaaacgacagttaccgccatagacaaggattcatttgtgaatgataaataaaagttacccatatcaagaccaaataatggatcaatgaaaaattgatccaaggggcttgaaattactggtccgccgcctccagggataactatatcgccgttcctccattcagccccatctaacaactctattccgttctgccagttattgaaacggatacgctgcgccctctcaaatggttcaaaaaaagcgtctatttctgcttgtctccttctcttcaacgagttcacgttataatcattgaggagtttccagatcgagtctctctgattcaaacccattattaccacatgattgtgatcgaaattcat"""

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
        print(motif_df.head())

if __name__ == "__main__":
    test_run()
