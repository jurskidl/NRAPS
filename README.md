# Nuclear Reactor Analysis Particle Simulator (NRAPS)
This project started as a class project to simulate reactor physics. This entailed building both a Monte Carlo (MC) solver and a deterministic solver. Although the class has ended, I want to improve the vode base as a whole in order to continue the learning and improvement process and I will therefore continue making commits until I am fully satisfied with the state of the project.
The project consists of several peices that work together to function. These parts are:
  1. Input Processing
  2. Mesh Generation
  3. Physics Solver
     - Monte Carlo
     - Deterministic
  4. Plotting Results

## Input Processing
A lot of work went into the input pocessing to ensure it ran quickly at the time the project was submitted. Currently, the code processes the input file using a memmap that enables it to read the file extremely quickly. 
TODO: add assumptions made on the input file formatting.
Multithreading was tested and was relatively inconsequential as evidenced in the run time figure below. The time it takes to spin up multiple threads takes more time than the additional threads can reduce the processing time.
<img width="600" height="371" alt="Multithread_input" src="https://github.com/user-attachments/assets/683fe816-d47d-4d66-93f5-cebd6bae0f89" />
Input parsing has a negligible impact on runtime compared to file-system overhead. Total execution time scales near-constantly until input length exceeds 110,000 times the base size, after which parsing overhead becomes statistically significant. This is evidenced in the image below.
<img width="600" height="371" alt="Input_lines" src="https://github.com/user-attachments/assets/221513b1-043e-4a05-8cab-e08103439e83" />
The table below indicates the intended improvements to this portion of the project:
|                            Task description                              |        Status      |
| :----------------------------------------------------------------------- | :----------------: |
| Implement array to save data instead of hashing since data set is known. | :heavy_check_mark: |
| Fully define input file specification                                    |         :x:        |
| Implement Chunking to contol programs memory usage                       |         :x:        |

## Mesh Generation
This portion sets up the fidelety with which the problem will be solved and can be set separately for water and fuel regions. This doesn't have specific improvements, but generally speaking exploring SIMD implementations and other optimizations are of interest. If specifics are determined, they will be added into a table below.

## Physics Solver
### Monte Carlo
The primary motivation to improve the monte carlo code is to remove the dependence on the rand crate. Since the goal is statistical quality and not cryptography, the two methods considered for pseudo-random number generation are [Xoshiro128++](https://prng.di.unimi.it/) and [PCG-Random](https://www.pcg-random.org/). PCG-Random was chosen as the base of the random number generator due to its better statistical quality at a relatively small computational cost.
Other items to change include refactoring to improve testability and ultimately exploring the costs/benefits of a compute shader for the particle simulation. These are summarized int the table below.
|                            Task description                              |        Status      |
| :----------------------------------------------------------------------- | :----------------: |
| Develop self-maintained psuedo-random number generator                   | :heavy_check_mark: |
| Refactor code to improve readability and maintainability                 |         :x:        |
| Create unit tests that ensure accuracy and validity of segments of code  |         :x:        |
| Create GPU compute shader to perform particle calulcations               |         :x:        |
To implement the GPU compute shader, numerous steps will need to be completed to ensure the pipeline runs efficiently. This requires some sort of [event based algorithm](https://www.sciencedirect.com/science/article/pii/S1738573317302966#sec4) to vectorize the solutions to prevent branch divergence. There will also need to be implementation (Woodcock tracking](https://www.yiningkarlli.com/projects/specdecomptracking/references/Woodcock1965.pdf)
### Deterministic Solver
The Deterministic solver currently relies on the NAlgebra crate, whereas the goal is to develop an internal matrix solver. This internal solution method should be able to solve the problem in 4 different ways: 1) direct matrix inversion, 2) Gaussian Elimination, 3) Jacobi Method, and 4) successive over-relaxation (SOR).
|                            Task description                              |        Status      |
| :----------------------------------------------------------------------- | :----------------: |
| Develop self-maintained matrix solver                                    | :heavy_check_mark: |
| Refactor code to improve readability and maintainability                 |         :x:        |
| Create unit tests that ensure accuracy and validity of segments of code  |         :x:        |
| Develop all four solution methods                                        |         :x:        |
## Plotting Results
The current plotting method saves the data to a file and executes python code via the terminal to plot using the matplotlib library in python. This should at the very least be transferred over to rust and run via the PyO3 crate.
