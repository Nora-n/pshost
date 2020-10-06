# ps-host
## Requirements
This module uses `pymage` and `astrobject`. To install, `pip install pymage` and
`pip install astrobject`

## Installation
Do wherever you want to save the folder and then
```bash
git clone https://github.com/Nora-n/pshost.git
cd pshot
python setup.py install
```

## Imports
For now, only `massmeasure` is defined. To import it:
```python
from pshost import massmeasure
```

## Usage
From a `ra`, `dec`, and `z`, simply instantiate the class with:
```python
mm = massmeasure.MassMeasure(ra, dec, z)
```
Then, to get the host galaxy parameters:
```python
mm.get_hostgalaxy_params()
```
which returns only the mass estimation with asymmetric errors, but stores
information as properties of `mm`.

## To be done
- [ ] Write properties
- [ ] Detail documentation
