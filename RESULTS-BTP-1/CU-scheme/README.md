# 1-Dimensional Central Upwind Scheme

To run the code correctly, make sure the following files are present in the directory

- main.cpp
- 1D_CentralUpwind.h
- CUNumericalFlux.h
- ExtendCells.h
- InitializeRiemannProblem.h
- SSPRK.h
- Constants.h
- OutputResult.h
- PrimitiveVariables.h
- Utility.h
- plot.py

--- 

## Method-1 Running in a Linux or related environment (any bash terminal)

- Use the *run.sh* bash file to run the code in a single command

---

## Method-2 Running on non-unix environment

- STEP-1 Create the following directories in the root with the same spelling (if present make sure these directories are empty)
    - env
    - plots
    - result
    - result1
    - result2

- STEP-2 Run the following commands in order
    - g++ main.cpp -o prg.x
    - run the executable ***prg.x***

- STEP-3 Plot the results with the following command
    - python plot.py ( or python3 plot.py )

---