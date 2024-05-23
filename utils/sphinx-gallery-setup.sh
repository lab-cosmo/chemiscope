#!/bin/bash

# Setup npm packages and build js part
echo "Running npm install and npm run build..."
npm install
npm run build

# Check if npm build was successful
if [ $? -ne 0 ]; then
  echo "npm build failed. Exiting script."
  exit 1
fi

# Install Python requirements
echo "Navigating to docs folder..."
cd ../docs || { echo "docs folder not found. Exiting script."; exit 1; }

echo "Installing Python requirements..."
pip install -r requirements.txt

# Check if pip install was successful
if [ $? -ne 0 ]; then
  echo "Failed to install Python requirements. Exiting script."
  exit 1
fi

# Build Sphinx Gallery html files
echo "Running make html..."
make html

# Check if make html was successful
if [ $? -ne 0 ]; then
  echo "Failed to generate HTML documentation. Exiting script."
  exit 1
fi

echo "Setup executed successfully."
