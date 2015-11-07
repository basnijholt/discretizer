git config core.attributesfile .gitattributes
git config filter.dropoutput_ipynb.clean "python scripts/ipynb_output_filter.py"
git config filter.dropoutput_ipynb.smudge cat
