
| File                      | Content                                                                                               |
| --------------------------| ------------------------------------------------------------------------------------------------------|
| src/header.h              | A header file containing all function declarations and constant definitions                           |
| src/time_integration.cpp  | A source file containing the time integration function                                                |
| src/sys_physics.cpp       | A source file containing the functions that initialize the lattice or calculate physical quantities   |
| src/io.cpp                | A source file containing all the functions needed for reading input data and writing output data      |
| src/main.cpp              | The main program source file that calls the functions and performs the simulation                     |

Compilation: (Use makefile in the project directory)
```
make moldy
```
This will create an executable file called moldy in the run directory.

To erase all executables and .o -files, type make clean into the command line

