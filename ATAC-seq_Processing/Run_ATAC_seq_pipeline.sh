#!/bin/bash

# This script runs the looper command with the specified configuration file.
# Usage: ./run_looper.sh
# Make sure that 'looper' is installed and accessible in your PATH.

# Check if looper is installed
if ! command -v looper &> /dev/null
then
    echo "Error: looper is not installed. Please install it to proceed."
    exit 1
fi

# Define the configuration file
CONFIG_FILE="Metastasis.yaml"

# Run the looper command
echo "Running looper with configuration file: $CONFIG_FILE"
looper run $CONFIG_FILE

# Check if the looper command was successful
if [ $? -eq 0 ]; then
    echo "Looper command executed successfully."
else
    echo "Error: Looper command failed."
    exit 1
fi





