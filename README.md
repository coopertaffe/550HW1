# 550HW1
MATLAB code to solve the Stokes equations in 2D on a staggered grid.
## Contents
- `550HW1.pdf` - Compiled Homework PDF
- `HW1.tex` - LaTeX file
- `figures/` - jpegs of figures used in Homework report
- `MATLAB code/` - code used to create figures

## Usage
#### View Homework
- Open `550HW1.pdf` to read the homework.

To compile the homework from its LaTeX source, run:
```bash
latexmk HW1.tex
```

#### Run Code
1. Clone or download this repository.
2. Open the MATLAB scripts in MATLAB.
3. Run the desired script to reproduce the numerical results.
   
`HW1_part1.m` solves stokes and produces plot of the magnitude of **u**. User can vary grid size by changing variable N.
`HW1_ERROR.m` gives a plot of relative error based on difference grid sizes. User can vary vector of grid sizes by changing variable Nv.
In `HW1CharliesRule.m` : section 1 times backslash and GMRES and section 2 changes the tolerance of GMRES. 

