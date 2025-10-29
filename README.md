# Advanced Enzyme Kinetics Database Analyzer

This project is a Python-based graphical user interface (GUI) application designed for the comprehensive analysis of enzyme kinetics data. It allows users to load, filter, visualize, and analyze data related to enzyme mutations, kinetic parameters (kcat, Km), and associated protein structures.

## 🔬 Features

* **Interactive GUI:** Built with Tkinter, providing a user-friendly interface for data interaction.
* **Data Analysis:** Utilizes `pandas` for powerful data manipulation and analysis of enzyme kinetics data from CSV files.
* **Rich Visualization:** Generates various plots (scatter plots, bar charts, statistical analyses) using `matplotlib` and `seaborn`.
* **3D Protein Visualization:** Integrates `py3Dmol` and `BioPython` to fetch, display, and analyze 3D protein structures from the PDB.
* **Data Export:** Allows users to export filtered data and analysis results to Excel (`.xlsx`) files with custom formatting.
* **Statistical Analysis:** Includes statistical tools from `scipy` to compare datasets.

## 📋 Requirements

This application requires Python 3 and the following libraries:

* `numpy`
* `matplotlib`
* `requests`
* `biopython`
* `py3Dmol`
* `Pillow (PIL)`
* `pandas`
* `openpyxl`
* `seaborn`
* `scipy`

## ⚙️ Installation

1.  **Clone the repository:**
    ```bash
    git clone <your-repository-url>
    cd <your-repository-name>
    ```

2.  **Create a virtual environment (recommended):**
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
    ```

3.  **Install the required libraries:**
    You can create a `requirements.txt` file with the libraries listed above and install them all at once:
    ```bash
    pip install -r requirements.txt
    ```
    *(Alternatively, install them manually: `pip install numpy matplotlib requests biopython py3Dmol pillow pandas openpyxl seaborn scipy`)*

## 🚀 Usage

1.  **Ensure your data file is present:**
    This program is designed to work with the `all_enzyme_kcat_km_ddG_all_clean_pH.csv` file. Make sure this file is in the same directory as the Python script.

2.  **Run the application:**
    ```bash
    python "file data edited.py.py"
    ```

3.  **Using the Application:**
    * The main window will open.
    * Use the various tabs and buttons to load data, apply filters, generate plots, and view protein structures.
    * Analysis results and plots can be viewed within the app or exported.

## 📁 Data

The core dataset for this application is `all_enzyme_kcat_km_ddG_all_clean_pH.csv`. This file contains detailed information on enzyme mutations, organisms, kinetic parameters, and corresponding protein data.

**Note:** If this CSV file is very large (e.g., >100MB), you should use **Git LFS (Large File Storage)** or host the file elsewhere (like Google Drive or Zenodo) and provide a download link in this README. GitHub has a strict file size limit.
