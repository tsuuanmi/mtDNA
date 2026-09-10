#!/bin/bash

# Creates improved directory structure for mtDNA project
BASE_DIR="$(pwd)"

# Create main directories
mkdir -p "${DATA_DIR}/raw"
mkdir -p "${DATA_DIR}/metadata"

mkdir -p "${RESULTS_DIR}/automate_pipeline/blastn"
mkdir -p "${RESULTS_DIR}/automate_pipeline/tracy"
mkdir -p "${RESULTS_DIR}/automate_pipeline/comparison"
mkdir -p "${RESULTS_DIR}/manual_pipeline/regenerate"

mkdir -p "${RESULTS_DIR}/manual_pipeline/manual"
mkdir -p "${RESULTS_DIR}/manual_pipeline/comparison"
mkdir -p "${RESULTS_DIR}/manual_pipeline/regenerate"

mkdir -p "${RESULTS_DIR}/fasta"
mkdir -p "${RESULTS_DIR}/reports"

mkdir -p "${BASE_DIR}/logs"

# Set appropriate permissions
chmod -R 755 "${BASE_DIR}"

# Check and download tools
TOOLS_DIR="${BASE_DIR}/tools"
mkdir -p "${TOOLS_DIR}"

# Check for UGENE
if [ ! -d "${TOOLS_DIR}/ugene-50.0" ]; then
    echo "Downloading UGENE..."
    wget https://github.com/ugeneunipro/ugene/releases/download/50.0/ugene-50.0-linux-x86-64.tar.gz -P "${TOOLS_DIR}/"
    cd "${TOOLS_DIR}"
    tar -xzf ugene-50.0-linux-x86-64.tar.gz
    rm ugene-50.0-linux-x86-64.tar.gz
    cd "${BASE_DIR}"
else
    echo "UGENE already installed"
fi

# Check for Tracy
if [ ! -f "${TOOLS_DIR}/tracy" ]; then
    echo "Downloading Tracy..."
    wget https://github.com/gear-genomics/tracy/releases/download/v0.7.8/tracy_v0.7.8_linux_x86_64bit -P "${TOOLS_DIR}/"
    chmod +x "${TOOLS_DIR}/tracy_v0.7.8_linux_x86_64bit"
    mv "${TOOLS_DIR}/tracy_v0.7.8_linux_x86_64bit" "${TOOLS_DIR}/tracy"
else
    echo "Tracy already installed"
fi

# Create virtual environment
curl -LsSf https://astral.sh/uv/install.sh | sudo sh

uv venv --python 3.12

uv sync

# Check if .env exists, if not create it
if [ ! -f ".env" ]; then
    cp .env.example .env
fi