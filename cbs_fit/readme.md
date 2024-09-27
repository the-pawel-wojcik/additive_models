# Fitting EOMEE energies

The `find_cbs.py` program takes a JSON file as an input and runs a complete
basis set (CBS) extrapolation. The program returns a JSON with the CBS values.
It also prints various characteristics on the way and displays fitting graphs.

## Prepare input for the fitting
Go to the directory with your CFOUR outputs and run:
```bash
cfour_parser -j output.c4 > output.json
processors/print_roots.py -c output.json | jq > cbs_input_basis.json
```
If you are using the CFOUR's `basis = SPECIAL`, you will need to manually
adjust the value in the output.

Merge all the `cbs_input_basis.json` files into a single input with 
```bash
find . -name 'cbs_input_*' -exec cat '{}' + | jq -s '.' > cbs_input.json
```
I use the convention that the single input file name marks the CC level and the
basis set family, i.e. in the last command `cbs_input = ccsd+pwCVnZ.json`.

## Running CBS extrapolation 
Run your input through the `find_cbs.py` script. Save the output as the same
file name with the "+cbs" suffix, e.g., `ccsd+pwCVnZ+cbs.json`.

# Specific to EOMEE

## See what you got
Use the `pprint_final_energies.py` to pretty print the final energies.
```bash
./pprint_final_energies.py -s cbs_input{,+cbs}.json
```

## Use the CBS energies in your xsim simulation 
Use the `turn_cbs_into_xsim_input.py` script.
