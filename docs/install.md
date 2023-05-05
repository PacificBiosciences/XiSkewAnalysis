# Installing XiSkewAnalysis
## Instructions
The following instructions assume familiarity with command line and that both `git` and `conda` are installed.

1. Clone this repository into a local folder: `git clone https://github.com/PacificBiosciences/XiSkewAnalysis.git`.
2. Create and activate a conda environment for the XiSkewAnalysis.
4. Test the script file by running it with the help option (`-h`).
5. Visit the [User guide](./user_guide.md) for details on running XiSkewAnalysis.

## Example
```bash
# clone the repository
git clone https://github.com/PacificBiosciences/XiSkewAnalysis.git
cd XiSkewAnalysis
# create and activate the conda environment
conda env create --file ./conda/xci_skew.yaml --name xci_skew
conda activate xci_skew
# execute help instructions
python3 ./scripts/SampleXCI.py -h
```