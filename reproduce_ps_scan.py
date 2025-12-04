import os
import sys
from ppr_analysis import motif_generation

# Dummy data for testing
PROTEIN_SEQ = """>PPR791(Rf1a)_IR24
MARRAASRAVGALRSDGSIQGRGGRAGGSGAEDARHVFDELLRRGRGASIYGLNRALADVARDSPAAAVSRYNRMARAGADEVTPDLCTYGILIGCCCRAGRLDLGFAALGNVIKKGFRVDAIAFTPLLKGLCADKRTSDAMDIVLRRMTELGCIPNVFSYNILLKGLCDENRSQEALELLHMMADDRGGGSPPDVVSYTTVINGFFKEGDSDKAYSTYHEMLDRGILPDVVTYNSIIAALCKAQAMDKAMEVLNTMVKNGVMPDCMTYNSILHGYCSSGQPKEAIGFLKKMRSDGVEPDVVTYSLLMDYLCKNGRCMEARKIFDSMTKRGLKPEITTYGTLLQGYATKGALVEMHGLLDLMVRNGIHPDHYVFSILICAYAKQGKVDQAMLVFSKMRQQGLNPNAVTYGAVIGILCKSGRVEDAMLYFEQMIDEGLSPGNIVYNSLIHGLCTCNKWERAEELILEMLDRGICLNTIFFNSIIDSHCKEGRVIESEKLFELMVRIGVKPNVITYNTLINGYCLAGKMDEAMKLLSGMVSVGLKPNTVTYSTLINGYCKISRMEDALVLFKEMESSGVSPDIITYNIILQGLFQTRRTAAAKELYVRITESGTQIELSTYNIILHGLCKNKLTDDALQMFQNLCLMDLKLEARTFNIMIDALLKVGRNDEAKDLFVAFSSNGLVPNYWTYRLMAENIIGQGLLEELDQLFLSMEDNGCTVDSGMLNFIVRELLQRGEITRAGTYLSMIDEKHFSLEASTASLFIDLLSGGKYQEYYRFLPEKYKSFIESLSC"""

def test_ps_scan():
    print("Testing PS_Scan API...")
    try:
        # Create temp file for protein seq
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fasta") as tmp:
            tmp.write(PROTEIN_SEQ)
            tmp_path = tmp.name
            
        # Run PS_Scan API
        if not os.path.exists("test_ps_scan_output"):
            os.makedirs("test_ps_scan_output")
            
        pkl_path = motif_generation.run_ps_scan_api(tmp_path, "test_ps_scan_output")
        print(f"PS_Scan API Result PKL: {pkl_path}")
        
        # Clean up
        os.remove(tmp_path)
    except Exception as e:
        print(f"PS_Scan API Failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_ps_scan()
