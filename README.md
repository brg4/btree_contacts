# btree_contacts

## Description

**A fortran program to calculate contacts for nested theta structures of replicating chromosomes.**

Gilbert et al., *Frontiers in Cell and Developmental Biology*, 2023, DOI:[10.3389/fcell.2023.1214962](https://doi.org/10.3389/fcell.2023.1214962)

## Repository directory structure

 - `src` - *source files and Makefile for program*
 - `example` - *example to demonstrate program*

## Installation

Program installation requires a fortran compiler. The program was tested with gfortran-v12.1.0.

Run `make contact_calc` in the `src` directory. The executable (**contact_calc**) will be in `src`.

## Usage

Prepare an input file.

Run with: `./contact_calc --i_f=<input file> --o_d=<output directory> --o_l=<output label> --l=<log file> --n_t=<num threads>`

Output:

- `rval_(label).bin` - *binary file storing 1D array of elements ($A_{ij}$, double) of symmetric matrix*
- `rval_mapped_(label).bin` - *binary file storing 1D array of elements ($A_{ij}$, double) of mapped symmetric matrix*
- `pairs_(label).bin` - *binary file storing 2D array of integer pairs encoding indices ($(i,j)$, int32) of symmetric matrix elements*
- `pairs_mapped_(label).bin` - *binary file storing 2D array of integer pairs encoding indices ($(i,j)$, int32) of symmetric matrix elements*

See `example` directory for a demonstration.

## Support
brg4@illinois.edu

## Authors and acknowledgment
Benjamin R. Gilbert - brg4@illinois.edu

## Project status
This project is under development.
