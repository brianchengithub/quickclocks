# GitHub Repository Setup Guide

This guide explains how to set up a GitHub repository for the Epigenetic Clock Calculator using GitHub Desktop.

## Repository Structure

Your repository should have this structure:

```
EpigeneticClockCalculator/
├── .gitignore              # Files to exclude from version control
├── DESCRIPTION             # R package metadata
├── LICENSE                 # MIT License
├── NAMESPACE               # R package exports
├── README.md               # Main documentation
├── install_dependencies.R  # Dependency installation script
├── diagnose_packages.R     # Troubleshooting script
├── R/                      # R source code
│   ├── calculate_clocks.R      # Main function
│   ├── clock_availability.R
│   ├── clock_calculator.R
│   ├── clock_implementations.R
│   ├── dunedin_pace.R
│   ├── epitoc2.R
│   ├── imputation.R
│   ├── input_handler.R
│   ├── methylcipher.R
│   ├── output_handler.R
│   ├── pc_clocks.R
│   ├── preprocessing.R
│   ├── sesame_functions.R
│   ├── utils.R
│   └── zzz.R
├── inst/
│   └── extdata/
│       └── reference_betas.rds  # Reference data for imputation (~8MB)
└── man/                    # Documentation (auto-generated)
```

## Setting Up with GitHub Desktop

### Step 1: Install GitHub Desktop

Download and install from: https://desktop.github.com/

### Step 2: Create a New Repository

1. Open GitHub Desktop
2. Click **File → New Repository** (or Ctrl+N / Cmd+N)
3. Fill in the details:
   - **Name**: `EpigeneticClockCalculator` (or your preferred name)
   - **Description**: "Unified calculator for 40+ epigenetic clocks"
   - **Local Path**: Choose where to save the repository
   - **Initialize with README**: ❌ Uncheck (we have our own)
   - **Git Ignore**: Select "R" from dropdown
   - **License**: MIT License
4. Click **Create Repository**

### Step 3: Copy Your Files

1. Open the repository folder (click "Show in Explorer/Finder")
2. Copy all the package files into this folder:
   - Copy `R/` folder
   - Copy `inst/` folder
   - Copy `DESCRIPTION`, `NAMESPACE`, `LICENSE`
   - Copy `README.md`
   - Copy `install_dependencies.R`, `diagnose_packages.R`
   - Replace the auto-generated `.gitignore` with ours

### Step 4: Commit Your Files

1. Return to GitHub Desktop
2. You should see all files listed under "Changes"
3. Review the changes (make sure no large data files are included)
4. Write a commit message:
   - **Summary**: "Initial commit: Epigenetic Clock Calculator v2.0.0"
   - **Description**: (optional) "40+ epigenetic clocks with EpiDISH cell deconvolution"
5. Click **Commit to main**

### Step 5: Publish to GitHub

1. Click **Publish repository** (top right)
2. Choose settings:
   - **Name**: Keep as is or modify
   - **Description**: Keep as is or modify
   - **Keep this code private**: ✅ Check if you want private repo
3. Click **Publish Repository**

Your repository is now on GitHub!

## Important: Check File Sizes

Before committing, make sure these large files are NOT included:

| File | Size | Should Commit? |
|------|------|----------------|
| `inst/extdata/reference_betas.rds` | ~8 MB | ✅ Yes (needed for imputation) |
| `PCClocks_data.qs2` | ~2 GB | ❌ No (downloaded automatically) |
| `*.idat` files | Variable | ❌ No (user data) |
| Manifest cache files | ~15 MB each | ❌ No (downloaded automatically) |

The `.gitignore` file should prevent large files from being committed, but always verify!

## Updating the Repository

### Making Changes

1. Edit files in your local repository folder
2. Open GitHub Desktop - changes appear automatically
3. Review changes in the "Changes" tab
4. Write a descriptive commit message
5. Click **Commit to main**
6. Click **Push origin** to upload changes

### Pulling Updates

If you work on multiple computers:
1. Open GitHub Desktop
2. Click **Fetch origin** to check for updates
3. Click **Pull origin** if updates are available

## Recommended GitHub Settings

### Repository Settings (on github.com)

1. Go to your repository on github.com
2. Click **Settings**

**Recommended settings:**

- **General → Features**: Enable "Issues" for bug tracking
- **General → Features**: Enable "Discussions" for Q&A
- **Branches → Branch protection**: Consider protecting `main` branch

### Add Topics

On your repository page, click "Add topics" and add:
- `epigenetics`
- `dna-methylation`
- `r-package`
- `bioinformatics`
- `epigenetic-clock`
- `aging`

### Add a Release

After initial setup:
1. Click **Releases** (right sidebar)
2. Click **Create a new release**
3. Tag: `v2.0.0`
4. Title: "v2.0.0 - Initial Release"
5. Description: Copy key features from README
6. Click **Publish release**

## Installation Instructions for Users

After publishing, users can install with:

```r
# From GitHub
devtools::install_github("yourusername/EpigeneticClockCalculator")

# Or if private repository
devtools::install_github("yourusername/EpigeneticClockCalculator", 
                         auth_token = "your_personal_access_token")
```

## Troubleshooting

### "File too large" Error

If you accidentally try to commit a large file:
1. The `.gitignore` should prevent this
2. If it happens, remove the file from the commit
3. In GitHub Desktop: right-click the file → "Discard changes"

### "Repository not found" Error

- Check if the repository is private and you're logged in
- Verify the repository name/URL is correct

### Syncing Issues

1. Click **Repository → Repository Settings**
2. Verify the remote URL is correct
3. Try **Fetch origin** then **Pull origin**

## Keeping Dependencies Updated

Periodically update your package dependencies:

```r
# Check for updates
old.packages()

# Update Bioconductor packages
BiocManager::install(ask = FALSE)

# Update GitHub packages
devtools::install_github("danbelsky/DunedinPACE")
devtools::install_github("HigginsChenLab/methylCIPHER")
```

---

For more help with GitHub Desktop, see: https://docs.github.com/en/desktop
