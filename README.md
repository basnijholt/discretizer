# Left to be done
## Main stuff
* fixing recursive discretization to work separately on every summand
* shortening and reading hoppings (super important)
* generalazing it to work on matrices (important)
* transforming it into functions (check how Anton did it)
* writing tests, lot of tests (maybe wait till packaging)
* making package
* wrapper of string input (nice, but least on priorty)

# Development guide


## Ipymd
Please use ``ipymd`` to work on markdown ``md`` instead of ``ipynb`` files.
If you use your own laptop please do:
```
pip install git+https://github.com/rossant/ipymd
```
and add to ~/.jupyter/jupyter
```
c.NotebookApp.contents_manager_class = 'ipymd.IPymdContentsManager'
```

In case it doesn't work ask Rafal. On hpc1 should work by default (dev2 tested)


## If we decide to use ``ipynb`` notebooks we may need output filters.
Let's activate them now, just in case, with:
```
sh activate_filters.sh
```
