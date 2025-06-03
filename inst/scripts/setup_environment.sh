#!/bin/bash
# setup_environment.sh

# Exit on error
set -e

# Define variables
ENV_PREFIX="/share/maize/frodrig4/conda/env/bzeaseq"
YAML_FILE="bzeaseq.yml"
BASE_DIR="/rsstu/users/r/rrellan/BZea/bzeaseq"

# Print section header
print_header() {
    echo "============================================"
    echo "  $1"
    echo "============================================"
}


print_header "Setting up conda environment for bzeaseq"

# Create conda environment from YAML file
print_header "Creating conda environment from YAML"
conda env create -p $ENV_PREFIX -f $YAML_FILE -y
print_header "Activating environment"
eval "$(conda shell.bash hook)"
conda activate $ENV_PREFIX


# Create environment activation script
print_header "Creating activation script"
cat > activate_bzeaseq.sh << 'EOF'
#!/bin/bash

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /share/maize/frodrig4/conda/env/bzeaseq


# Set environment variables for pipeline
# Update this
export B73_REFERENCE="/rsstu/users/r/rrellan/BZea/bzeaseq/Zm-B73-REFERENCE-NAM-5.0.fa" 

# Print environment info
echo "Maize genomics environment activated"
echo "Python: $(which python)"
echo "R: $(which R)"
echo "GATK: $(which gatk)"
echo "BWA: $(which bwa)"
EOF

chmod +x activate_bzeaseq.sh

print_header "Testing the environment"
source activate_bzeaseq.sh

# Test key tools
echo "Testing key tools..."
python --version
R --version
bwa || true
samtools --version
bcftools --version
gatk --help || true
bedtools --version

print_header "Setup complete!"
echo "To activate the environment, run:"
echo "source activate_bzeaseq.sh"