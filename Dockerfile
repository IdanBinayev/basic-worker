# Use an official Python runtime as a parent image
FROM python:3.9

# Set the working directory in the container
WORKDIR /app

# Copy the model and the script into the container
COPY xgb_model.pkl onehot_encoder.pkl predictor.py ./

# Install dependencies
RUN pip install rdkit duckdb xgboost scikit-learn pandas numpy joblib

# Set the command to run the predictor script
CMD ["python", "predictor.py"]
