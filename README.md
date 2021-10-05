# GCOM-2021-0117

## **On mild solutions of the Laplacian fractional Langevin equations with anti-periodic type boundary conditions**

The python codes here are to perform the numerical experiments:

+ Example **4.1**, **(a)** and **(b)** 
+ Example **4.2**
+ Example **4.3** 
  + **(a)** without noise to input
  + **(b)** with random noise to input 

To do, we have 5 tasks. Three first tasks should be done at first and only once. We will focus only on two last ones. 

## **1. Install Anaconda for Python 3** 

  a) Windows: 
  
  + https://repo.anaconda.com/archive/Anaconda3-2021.05-Windows-x86_64.exe
  + Guide: 
      + https://docs.anaconda.com/anaconda/install/windows/
      + https://www.datacamp.com/community/tutorials/installing-anaconda-windows

  b) Linux: 
  + https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
  + Guide:     
    + https://docs.anaconda.com/anaconda/install/linux/

## **2. Install Python libraries**

  a) Open a Command Prompt (Windows) or Terminal (Linux)
  
  b) Run the following commands
  
    conda install numpy
    conda install scipy
    conda install matplotlib

## **3. Obtain this package**

  Download https://codeload.github.com/tranqv/GCOM-2021-0117/zip/refs/heads/main
  
  and extract the zip file. One should have **GCOM-2021-0117-main** folder.
  
  Inside the subdirectory **code** of the folder **GCOM-2021-0117-main**, three examples are prepared with their main programs as follows:
     
  + **main_ex01a.py**:  Example 1a,  
  + **main_ex01b.py**:  Example 1b, 
  + **main_ex02a.py**:  Example 2,  
  + **main_ex03a.py**:  Example 3a, without noise,
  + **main_ex03b.py**:  Example 3b, with noise to the input. 

  All the programs above perform **Procedure (P)** in Section 4 of the manuscript. Read comments included in and messages inside the **print()** function of the codes to gain information of the computation steps, and check with output to the screen when executing.
  
  Setting for each of the examples is given in file **inc_ex*.py**, i.e., **inc_ex1a.py**, **inc_ex1b.py**, **inc_ex2a.py**, **inc_ex3a.py**, **inc_ex3b.py**.
  
  Auxiliary procedures for the computation are defined in **inc_sub2.py**, while the ones for plotting are in **inc_post.py**.
 
  ## **4. Run the test cases**
  
  a) **Method 1:** Using Jupyter Notebook (installed together with Anaconda)
  
  Change directory to **GCOM-2021-0117-main**, and to **code**, then we click on one of the following files
  
  + **main_ex01a.ipynb**:  Example 1a,  
  + **main_ex01b.ipynb**:  Example 1b, 
  + **main_ex02a.ipynb**:  Example 2,  
  + **main_ex03a.ipynb**:  Example 3a, without noise,
  + **main_ex03b.ipynb**:  Example 3b, with noise to the input. 

  They are copied from the corresponding file **main_ex*.py** to run in Jupyter Notebook.

  Click to choose a file, move to the new tab, then click on _Kernel_ and _Restart & Run All_ to execute the program.

  b) **Method 2:** Using Anaconda Prompt (Windows) or Terminal (Linux)
  
  + Change directory to **GCOM-2021-0117-main**, and to **code**

  + Inside **code**, run the examples using the following commands 
  
  ```
      python3 main_ex01a.py
      python3 main_ex01b.py
      python3 main_ex02a.py
      python3 main_ex03a.py
      python3 main_ex03b.py
  ```
   where the first command is for Example 1a, and so on.

### **Notes:**
  + At the first time we run one of the above commands, two folders, i.e., **figs** and **xout**, will be created.
    + **figs** saves figures from all the main programs
    + **xout** saves the results (numerical mild soluions, error estimates) in ascii format, i.e., .txt and .csv.
  + The **main_ex03a.py** (or **main_ex03a.ipynb**) must be excecuted in advance of **main_ex03b.py** (or **main_ex03b.ipynb**) since the latter program needs the output from the preceeding one, while the others (for Examples 1a, 1b, 2, and 3a) can run independently.
  + Check output to your screen with the **log_*.txt** files.
      
## **5. Postprocessing**
   
  After we finish running the examples, we can run 
  
    main_post.ipynb
  
  with Jupyter Notebook (read **4. Method 2**), or invoke the command 
  
      python3 main_post.py
      
  in Anaconda Prompt (Windows) or Terminal (Linux) to plot the results of Examples 1a, 1b, 2, 3a. 
  
  Note that, in **main_post.py** (or **main_post.ipynb**), we must set the values of **fprefix** and **lmax** accordingly.
  While all figures for Example 3b were already made when we ran **main_ex03b.py** (or **main_ex03b.ipynb**).
    

  
## ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) **Hint:** Do this in advance to check if the codes work properly.

+ For Examples 1a, 1b, 2, and 3a, we should run **main.ipynb** or **main.py**.

    After opening this file (i.e. with Jupyter Notebook for **main.ipynb** or with a text editor for **main.py**)
    
    + set the value of **lmax** moderately, e.g. **lmax = 3**, 
    + assign the value of **fprefix** as it was prepared inside the code, e.g.
    ```python
    #fprefix = "ex01a"
    #fprefix = "ex01b"
    fprefix = "ex02a"       # remove the comment "#" to run Example 2
    #fprefix = "ex03a"
    ```
    + run the code as in 4a) Method 1 or 4b) Method 2.

+ For Example  3b, we should run **main_ex3b.py** or **main_ex3b.ipynb**.

    + set the value of **luref** to be small, e.g. **luref = 3**.

  Make sure that the folder **xout** does already exists and includes the results of Example 3a, i.e. files with the prefix **ex03a**. 
  
## Notes on Examples 1a and 1b:

  At the end of the computation process, some warnings about "divide by zero encountered" apprear. The reason is U_N = Uexact = 0 for all N, therefore all related error estimates are also zero. These messages should not crash the program. Please ignore them. 

## **References:**

+ [1] https://github.com/khinsen/mittag-leffler
+ [2] https://github.com/tranqv/Mittag-Leffler-function-and-its-derivative
+ [3] https://se.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function

Python libraries: 
+ [4] https://numpy.org/doc/stable/reference/
+ [5] https://docs.scipy.org/doc/scipy/reference/interpolate.html 
+ [6] https://docs.scipy.org/doc/scipy/reference/special.html
+ [7] https://matplotlib.org/
