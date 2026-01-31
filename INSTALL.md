# Installation Guide

## Requirements

- **Python**: 3.8 or higher
- **Operating System**: Linux, macOS, or Windows
- **Disk Space**: ~50 MB for installation

## Standard Installation

### Option 1: Using pip (Recommended)

```bash
# Navigate to project directory
cd /path/to/sanger_quality

# Install package
pip install .

# Verify installation
sangerqc --help
```

### Option 2: Development Mode

For development or if you plan to modify the code:

```bash
# Install in editable mode
pip install -e .

# Install with development dependencies
pip install -e ".[dev]"
```

### Option 3: Using pipx (Isolated Environment)

```bash
# Install pipx if not already installed
python3 -m pip install --user pipx
python3 -m pipx ensurepath

# Install sanger-qc-trim
pipx install /path/to/sanger_quality
```

## Verification

After installation, verify everything works:

```bash
# Check command is available
sangerqc --help

# Run verification script
python verify_algorithms.py

# Run test suite (if dev dependencies installed)
pytest tests/
```

## Dependencies

The following packages will be automatically installed:

### Required
- **biopython** ≥1.79 - Biological sequence parsing
- **pandas** ≥1.3.0 - Data manipulation
- **numpy** ≥1.21.0 - Numerical computing
- **pyarrow** ≥6.0.0 - Parquet file format
- **typer** ≥0.9.0 - CLI framework
- **tqdm** ≥4.62.0 - Progress bars

### Optional (Development)
- **pytest** ≥7.0.0 - Testing framework
- **pytest-cov** ≥3.0.0 - Test coverage
- **black** ≥22.0.0 - Code formatting
- **ruff** ≥0.1.0 - Linting

## Platform-Specific Notes

### Linux (Tested ✅)
No additional steps required.

### macOS
Should work without issues. If you encounter problems with dependencies:
```bash
# Install Xcode command line tools if needed
xcode-select --install

# Or use Homebrew Python
brew install python@3.10
```

### Windows
Recommended to use:
- **WSL2** (Windows Subsystem for Linux) - Recommended
- **Anaconda/Miniconda**
- **Native Windows** with Python from python.org

## Conda Installation (Alternative)

If you prefer conda:

```bash
# Create new environment
conda create -n sangerqc python=3.10

# Activate environment
conda activate sangerqc

# Install dependencies
conda install -c conda-forge biopython pandas numpy pyarrow tqdm

# Install package
pip install .
```

## Troubleshooting

### "Command not found: sangerqc"

1. Check if installation succeeded:
   ```bash
   pip list | grep sanger
   ```

2. Ensure pip bin directory is in PATH:
   ```bash
   echo $PATH
   ```

3. Reinstall in user mode:
   ```bash
   pip install --user .
   ```

### Permission Errors

Use `--user` flag:
```bash
pip install --user .
```

Or use a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install .
```

### Dependency Conflicts

Create a fresh virtual environment:
```bash
python -m venv fresh_env
source fresh_env/bin/activate
pip install --upgrade pip
pip install .
```

## Uninstallation

```bash
pip uninstall sanger-qc-trim
```

## Upgrading

```bash
# If installed normally
pip install --upgrade .

# If installed in development mode
cd /path/to/sanger_quality
git pull  # if using git
pip install --upgrade -e .
```

## Docker Installation (Advanced)

A Dockerfile is not currently provided, but you can create one:

```dockerfile
FROM python:3.10-slim

WORKDIR /app
COPY . /app

RUN pip install --no-cache-dir .

ENTRYPOINT ["sangerqc"]
```

Build and run:
```bash
docker build -t sangerqc .
docker run -v $(pwd)/data:/data sangerqc all /data -o /data/results
```

## Running Without Installation

If you prefer not to install:

```bash
# Run directly with Python
python -m sanger_qc_trim.cli --help

# Note: You'll need to install dependencies first
pip install biopython pandas numpy pyarrow typer tqdm
```

## Verification Checklist

After installation, verify:

- [ ] Command available: `sangerqc --help`
- [ ] Subcommands work: `sangerqc qc --help`
- [ ] Test data runs: `sangerqc all seq_data/Y19/ -o test_output`
- [ ] Tests pass: `pytest tests/` (if dev dependencies installed)
- [ ] Algorithms verified: `python verify_algorithms.py`

## Getting Help

If installation fails:

1. Check Python version: `python --version` (must be ≥3.8)
2. Update pip: `pip install --upgrade pip`
3. Check error messages in logs
4. Try creating a fresh virtual environment
5. Open an issue with error details

## System Requirements Summary

| Component | Requirement |
|-----------|-------------|
| Python | ≥3.8 |
| RAM | ~100 MB |
| Disk Space | ~50 MB |
| OS | Linux, macOS, or Windows |
| Internet | Required for installation |

## Quick Test

After installation:

```bash
# Quick test
sangerqc all seq_data/Y19/ -o quick_test
ls quick_test/qc/
cat quick_test/qc/summary.json
```

Success! You're ready to use `sangerqc`.
