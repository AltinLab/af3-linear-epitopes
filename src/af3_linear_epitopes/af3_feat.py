from mdaf3.AF3OutputParser import AF3Output


def extract_residue_pLDDT(row, inf_path):
    af3 = AF3Output(inf_path / row["job_name"])
    row["pLDDT"] = [
        float(atm_pLDDT.mean())
        for atm_pLDDT in af3.get_mda_universe().residues.tempfactors
    ]
    return row
