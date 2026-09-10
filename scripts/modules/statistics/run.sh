#!/bin/bash

set -e  # Exit on error
set -u  # Exit on undefined variable

# Load environment variables
source .env

# Step 2: Merge data and results
python -m src.modules.statistics.merge
