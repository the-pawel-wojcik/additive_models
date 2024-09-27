# How to use

1. Prepare input files 
```bash
cfour_parser -j output.c4 | jq > output.json
processors/print_roots.py -c output.json | jq > basis+cclevel.json
```
Do it for the CCSD and CCSDT versions.

2. Use the `pprint_dT.py` to print the ΔT correction.
