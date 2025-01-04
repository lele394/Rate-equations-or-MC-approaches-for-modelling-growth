# Rate equations or MC approaches for modelling growth

The following is a project submission for the Numerical Methods 2 unit of Master 2 CompuPhys. Documentation on the code can be found on `O10-Technical-Report`. Exloitations and results can be found in notebooks `001` to `004`. Finally, a benchmark can be found in notebook `012`.

## License

This work is provided under the [CC-BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/deed.en) licence.


# Notes on the project

## Notebooks

Jupyter notebooks can be found in the `Numerical Work` folder. They use the simulation implementation files to execute simulations and answer questions.

>#### Task notebooks
>
>- 001-notebook-task1
>- 002-notebook-task2
>- 003-notebook-task3
>- 004-notebook-task4

>#### Simulation related
>- 012-notebook-benchmark

## Markdown

Markdown reports are available in the numerical work. They are used when no computations are required.

>- 010-Technical-Report
>- 011-writeup_better_implementation


## Reports

Reports can be found the PDF formatin the `PDF reports` folder. Additionally, HTML versions can be found in `Numerical Work/HTML reports`. Notebooks can also be found and ran directly in `Numerical Work`.

Markdown files can be found in `HTML reports`.

`graphs` folder are identical and contain image ressources used in mardown inside of the notebook of task 2.


## Python files

Multiple python files can be found in `Numerical Work`.

#### Simulation related

- task2and3 : Single-threaded version of the simulation used in task 2 and 3.
- task2_parallel : Parallel version of the simulation used in task 2 and 3.
- task4_parallel : Parallel version of the simulation used in task 4.

#### Other files

- cluster : Sandbox file, code version used on the cluster. was modified multiple time.
- video_converter : can be used to convert hdf5 data from simulation to mp4 videos.
- test : self explanatory

## Saved data format.

Saved data are stored in hdf5 files using the h5 extension.
