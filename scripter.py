import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator
import requests


token = token()
inputs = data.get('inputs')
base_url = data.get('base')
trigger_compound = inputs.get('trigger_compound')
compound_name = trigger_compound.get('name')
compound_id = trigger_compound.get('id')
compound_smiles = trigger_compound.get('structure')

print(base_url)
print(f"Compound name: {compound_name}, Compound ID: {compound_id}")
print(f"SMILES: {compound_smiles}")

smiles_input = compound_smiles
protein_input = "BRD4"

# Load the pretrained model and encoder
xgb_model = joblib.load("/srv/scripter/xgb_model.pkl")
onehot_encoder = joblib.load("/srv/scripter/onehot_encoder.pkl")

def generate_ecfp(molecule, radius=2, bits=1024):
    """Generate ECFP fingerprint for a given molecule."""
    if molecule is None:
        return None
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=bits)
    return list(morgan_gen.GetFingerprint(molecule))

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

    return float(probability)
    
def update_bp_value(base_url, token, compound_id, value):
    """Update compound binding probability via PUT request."""
    url = f"{base_url}/api/v1/compounds/{compound_id}"
    
    # Create the payload according to the specified format
    payload = {
        "token": token,
        "item": {
            "BRD4 Binding Probability": value
        }
    }
    
    # Set headers
    headers = {
        "Content-Type": "application/json"
    }
    
    try:
        # Make the PUT request
        response = requests.put(url, json=payload, headers=headers)
        
        # Check if the request was successful
        response.raise_for_status()
        
        print(f"Successfully updated compound {compound_id}. Status code: {response.status_code}")
        return response.json() if response.text else None
        
    except requests.exceptions.RequestException as e:
        print(f"Error updating compound {compound_id}: {str(e)}")
        return None

predicted_prob = predict_single(smiles_input, protein_input)
print(f"🔬 Predicted Binding Probability: {predicted_prob:.4f}")

# Call the update function
result = update_bp_value(base_url, token, compound_id, predicted_prob)