# Development guide


## Ipymd
Please use ``ipymd`` to read markdown ``md`` instead of ``ipynb`` files.
If you don't use your own laptop please do:
```
pip install git+https://github.com/rossant/ipymd
```
and add to ~/.jupyter/jupyter
```
c.NotebookApp.contents_manager_class = 'ipymd.IPymdContentsManager'
```

In case it doesn't work ask Rafal.


## If we decide to use ``ipynb`` notebooks we may need output filters.
Let's activate them now, just in case, with:
```
sh activate_filters.sh
```
