import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


import os

# Print the current working directory
print("Current Path:", os.getcwd())
files = os.listdir()
print("Files in current directory:", files)

print("started with the predictor")


# Load the pretrained model and encoder
xgb_model = joblib.load("~/srv/scripter/xgb_model.pkl")
onehot_encoder = joblib.load("~/srv/scripter/onehot_encoder.pkl")

def generate_ecfp(molecule, radius=2, bits=1024):
    """Generate ECFP fingerprint for a given molecule."""
    if molecule is None:
        return None
    return list(AllChem.GetMorganFingerprintAsBitVect(molecule, radius, nBits=bits))

def predict_single(smiles: str, protein_name: str):
    """Predicts binding probability for a single molecule and protein pair."""
    
    # Convert SMILES to RDKit molecule
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid SMILES string")
    
    # Generate ECFP fingerprint
    ecfp = generate_ecfp(molecule)
    if ecfp is None:
        raise ValueError("Could not generate ECFP fingerprint")

    # One-hot encode the protein_name
    protein_onehot = onehot_encoder.transform([[protein_name]])[0]

    # Combine features into a single NumPy array
    X_test = np.array(ecfp + list(protein_onehot)).reshape(1, -1)

    # Predict probability
    probability = xgb_model.predict_proba(X_test)[0, 1]  # Probability of binding

    return probability

# Example usage
if __name__ == "__main__":
    smiles_input = "O=C(C[C@H](Nc1nc(NCc2nccs2)nc(Nc2ccc3c(c2)CNC3=O)n1)c1cccc(Cl)c1Cl)N[Dy]"
    protein_input = "BRD4"
    
    predicted_prob = predict_single(smiles_input, protein_input)
    print(f"🔬 Predicted Binding Probability: {predicted_prob:.4f}")