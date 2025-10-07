# Setting Up Jupyter Book with GitHub Pages for xQTL manuscript

Data and code to reproduce each figure in the xQTL manuscript

## 1. Installation and Environment Setup

Install Jupyter Book and related tools

```bash
pip install jupyter-book ghp-import
```

## 2. Initialize A Jupyter Book

Create a new Jupyter Book in the root directory of this project:

```bash
jupyter-book create website
```

## 3. Set Up Book Structure

The created template will have:
- `_config.yml`: Configuration for the book
- `_toc.yml`: Table of contents structure
- Sample content files

Let's customize these for this project:

### Custom _config.yml

Create or modify the `_config.yml` file:

```yaml
# Book settings
title: FunGen-xQTL Manuscript Repository
author: FunGen-AD Consortium
copyright: "2025"  # Current year

# Force re-execution of notebooks on each build
execute:
  execute_notebooks: force

# Define the name of the latex output file for PDF builds
latex:
  latex_documents:
    targetname: book.tex

# Add a bibtex file for citations
bibtex_bibfiles:
  - references.bib

# Information about where the book exists on the web
repository:
  url: https://github.com/StatFunGen/xqtl-paper
  branch: main

# Add GitHub buttons to your book
html:
  use_issues_button: true
  use_repository_button: true
```

### Custom _toc.yml

Create or modify `_toc.yml` file to include all figure notebooks and analysis sections:

```yaml
# Table of contents
# Learn more at https://jupyterbook.org/customize/toc.html

format: jb-book
root: README
chapters:
  - file: main_text/1_project_overview/index
  - file: main_text/2_single_context_cis/index
  - file: main_text/3_single_context_multigene_trans/index
  - file: main_text/4_multi_context/index
  - file: main_text/5_AD_xQTL_integration/index
  - file: main_text/6_AD_xQTL_genes/index
```

## 4. Create Content Files

The repository is now structured with the following organization:

### Main Text Structure
```
main_text/
├── 1_project_overview/
│   ├── index.md (FunGen-xQTL atlas and companion projects)
│   └── figure_data/README.md
├── 2_single_context_cis/
│   ├── index.md (Single context cis-xQTL fine-mapping)
│   └── figure_data/README.md
├── 3_single_context_multigene_trans/
│   ├── index.md (Multi-gene and trans-xQTL fine-mapping)
│   └── figure_data/README.md
├── 4_multi_context/
│   ├── index.md (Propagation of multi-context xQTL effects)
│   └── figure_data/README.md
├── 5_AD_xQTL_integration/
│   ├── index.md (Alzheimer's disease (AD) loci integration)
│   └── figure_data/README.md
└── 6_AD_xQTL_genes/
    ├── index.md (Insights into AD risk genes)
    └── figure_data/README.md
```

### Supplementary Materials
```
supplementary_tables/
└── README.md (Placeholder linking to Google Drive document)
```

### Root Files
- `README.md`: Landing page for the Jupyter Book
- `references.bib`: Bibliography file for citations

Each `figure_data` folder contains RDS format data files with the minimum information necessary to reproduce plot figures, allowing users to check numerical values.

## 5. Build the Book

```bash
jupyter-book build . --config website/_config.yml --toc website/_toc.yml
```

This will generate the HTML files in the `_build/html` directory.

## 6. Deploy to GitHub Pages

Use `ghp-import` to publish your book:

```bash
ghp-import -n -p -f _build/html
```

This will:
- Copy the contents of `_build/html` to a branch called `gh-pages`
- Push this branch to GitHub
- Configure it as a GitHub Pages source

## 7. Set Up GitHub Pages in Repository Settings

1. Go to your repository on GitHub
2. Navigate to Settings > Pages
3. Ensure the source is set to the `gh-pages` branch and the folder is `/ (root)`
4. Save the settings

The site should be live at: https://statfungen.github.io/xqtl-paper/

## 8. Automatic Deployment with GitHub Actions

The repository includes a GitHub Actions workflow (`.github/workflows/deploy-book.yml`) that automatically rebuilds and publishes your book whenever you push changes to the main branch.

The workflow:
- Sets up Python 3.9
- Installs jupyter-book and ghp-import
- Builds the book using your custom configuration
- Deploys to GitHub Pages

This ensures your website is always up-to-date with the latest changes in your repository.

## 9. Adding Content

To add new analyses or figures:

1. Create notebook files (.ipynb) or markdown files (.md)
2. Add them to the `_toc.yml` file in the appropriate section
3. Push to the main branch - the site will automatically rebuild

## Troubleshooting

- Ensure all files referenced in `_toc.yml` exist
- Check that notebook dependencies are properly specified
- Verify GitHub Pages is enabled in repository settings
- Check GitHub Actions logs for build errors