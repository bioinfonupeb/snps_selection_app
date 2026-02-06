# SNPs_QTL_Selection
A desktop application for filtering Single Nucleotide Polymorphisms (SNPs) from VCF files for Quantitative Trait Locus (QTL) analysis. This tool provides an intuitive graphical interface for loading, processing, and filtering genetic variant data.

# Features

- **User-Friendly GUI**: Built with Tkinter for easy interaction
- **Multiple Filtering Algorithms**: Five different filters for SNP selection
- **VCF File Support**: Standard Variant Call Format processing
- **Sample Configuration**: Flexible column mapping for different sample types
- **Real-time Logging**: Track processing steps and results
- **Export Results**: Save filtered variants to new VCF files


# Installation

## Prerequisites

- Python 3.6 or higher
- Tkinter (usually included with Python)

## Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/snps_qtl_selection.git
   cd snps_qtl_selection
    ```
2. No additional packages required! The application uses only Python standard libraries.

# Usage

## Running the Application

   ```bash
   python snps_qtl_selection.py
   ```

# Workflow

1. **Load VCF File**
    - Click "Browse" to select your VCF file
    - Click "Load VCF" to parse the file
    
2. **Configure Sample Columns**
    - The application auto-detects sample columns
    - Assign each of the 4 required sample types:
      - Parental Superior
      - Parental Inferior
      - Pool Superior
      - Random Pool
    - Click "Save Configuration" to store the mapping
    
3. **Select Filter**
   - Choose from 5 available filters:
     - **Minimum Allele Count**: Filter by minimum read count
     - **Percent Threshold**: Filter by allele frequency percentage
     - **Dominant Reference Allele**: Filter for dominant alleles
     - **Difference from Maximum Frequency**: Filter by frequency differences
     - **Random Pool Mean**: Filter using random pool statistics
     
4. **Configure Filter Parameters**
   - Enter required parameters for the selected filter
   - View help information for each filter type

5. **Apply Filter**
   - Click "Apply Filter" to process variants
   - View results in the log panel

6. **Save Results**
   - Click "Save Result" to export filtered variants
   - Choose location and filename