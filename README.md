# PharmVIP(Pharmacogenomic Variant Analysis and Interpretation Platform) - Guideline Package

### Install a local project in "editable" mode for develop package.

```shell
conda env create --file environment.yml
conda env export --no-builds | grep --invert-match "^prefix: " > environment.yml
conda env remove --name ENVIRONMENT_NAME
pip install -e .
```