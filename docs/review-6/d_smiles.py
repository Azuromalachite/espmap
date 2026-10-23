"""Does the fallback computed_smiles of the forced swapped-ethanol record parse?"""
from rdkit import Chem
for s in ["[CH][OH2]C", "CCO"]:
    m = Chem.MolFromSmiles(s)
    print(repr(s), "parses" if m is not None else "DOES NOT PARSE (RDKit, sanitised)")
    m2 = Chem.MolFromSmiles(s, sanitize=False)
    if m2 is not None:
        m2.UpdatePropertyCache(strict=False)
        print("   unsanitised: formal charge", Chem.GetFormalCharge(m2),
              "radicals", sum(a.GetNumRadicalElectrons() for a in m2.GetAtoms()),
              "O explicit valence", [a.GetExplicitValence() for a in m2.GetAtoms() if a.GetSymbol()=="O"])
