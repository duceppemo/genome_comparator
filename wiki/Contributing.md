# Contributing

Bug reports and pull requests are welcome on [GitHub](https://github.com/duceppemo/genome_comparator/issues).

## Development setup
```
conda env create -f environment.yml
conda activate genome_comparator
pip install --no-deps -e .
pytest
```
Tests run automatically on GitHub Actions for Python 3.10, 3.12 and 3.14.

## Documentation
This wiki is maintained in the [`wiki/`](https://github.com/duceppemo/genome_comparator/tree/master/wiki) folder of
the main repository. It is published to the GitHub wiki automatically when changes are pushed to `master`.
**Do not edit the wiki on GitHub directly**: changes would be overwritten. Edit the files in `wiki/` and open a pull
request instead.

* Page file names become page titles: `Output-files.md` → "Output files".
* Link to other pages without the extension: `[Usage](Usage)`.
* `_Sidebar.md` is the navigation shown on every page.
