# GCOM-2021-0117

#### **On mild solutions of the Laplacian fractional Langevin equations with anti-periodic type boundary conditions**

The python codes here are to perform numerical experiments of 

+ Example **4.1**, **(a)** and **(b)** 
+ Example **4.2**
+ Example **4.3** 
  + **(a)** without noise to input
  + **(b)** with random noise to input 


**1) Install Anaconda for Python 3** 

  a) Windows: 
  
  + https://repo.anaconda.com/archive/Anaconda3-2021.05-Windows-x86_64.exe
  + Guide: 
      + https://docs.anaconda.com/anaconda/install/windows/
      + https://www.datacamp.com/community/tutorials/installing-anaconda-windows

  b) Linux: 
  + https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
  + Guide:     
    + https://docs.anaconda.com/anaconda/install/linux/

**2) Install Python libraries**

  a) Open a Command Prompt (Windows) or Terminal (Linux)
  
  b) Run the following commands
  
    conda install numpy
    conda install scipy
    conda install matplotlib

**3) Obtain this package**

  Download https://codeload.github.com/tranqv/GCOM-2021-0117/zip/refs/heads/main
  
  and extract the zip file. One should have the folder **GCOM-2021-0117-main**.
  
  Inside the subdirectory **code** of the folder **GCOM-2021-0117-main**, three examples are prepared with their main programs as follows:
     
    + main_ex01a.py  Example 1a,  
    + main_ex01b.py  Example 1b, 
    + main_ex02a.py  Example 2,  
    + main_ex03a.py  Example 3a, without noise,
    + main_ex03b.py  Example 3b, with noise to the input. 

  All programs perform **Procedure (P)** in Section 4 of the manuscript. 
  
  Read comments inside the codes to gain informations and check with output to the screen when executing.
 
  **4) Run the test cases**
  
  a) Method 1: Using Jupyter Notebook
  
  Change directory to **GCOM-2021-0117-main**, and **code**, then we can see the following files:
  
      main_ex01a.ipynb
      main_ex01b.ipynb
      main_ex02a.ipynb
      main_ex03a.ipynb
      main_ex03b.ipynb
  
  Click to choose a file. Move to the new tab, click _Kernel_ and _Restart & Run All_ to execute the program.
  
  ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+)
  **Hint:** We should run **main.ipynb** firstly to check if the codes works properly.
  + set the value of **lmax** moderately, e.g. **lmax = 5**, and 
  + assig the value of **fprefix** as it was prepared inside the code, then 
  + _Kernel -> Restart & Run All_ 

  b) Method 2: Using Command Prompt (in Windows) or Terminal (in Linux)
  
  + Change directory to **GCOM-2021-0117-main**, and then the **code** inside 

  + Inside the folder **code**, run the examples as the following commands 
  
  ```
      python3 main_ex01a.py
      python3 main_ex01b.py
      python3 main_ex02a.py
      python3 main_ex03a.py
      python3 main_ex03b.py
  ```
   where the first command is for Example 1a, and so on.
  
  ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+)
  **Hint:** We should run the program **main.py** firstly to check if the codes works properly or not. 
  + set the value of **lmax** moderately, e.g. **lmax = 5**, and 
  + assign the value of **fprefix** as it was prepared inside the code, then 
  + perform this command as in Step 3d.
``` 
  python3 main.py
```

  **Notes:**
  + At the first time we run one of the above commands, two folders, i.e., **figs** and **xout**, will be created.
    + **figs** saves figures from all the main programs
    + **xout** saves the results (numerical mild soluions, error estimates) in ascii format, i.e., .txt and .csv.
  + The **main_ex03a.py** must be excecuted in advance of **main_ex03b.py** since the latter program needs the output 
  from the preceeding one, while the others (for Examples 1a, 1b, 2, and 3a) can run independently.
      
**5) Postprocessing**
   
  After running the examples, we can run 
  
    main_post.ipynb
  
  with Jupyter Notebook, or invoke the command 
  
      python3 main_post.py
      
  in a Command Prompt (in Windows) or Terminal (in Linux) to plot the results of Examples 1a, 1b, 2, 3a. 
  
  Note that, within **main_post.py**, we must set the values of **fprefix** and **lmax** accordingly.
  While figures for Example 3b were made when we run the command **python3 main_ex03b.py** as above.
    

**References:**

+ [1] https://github.com/khinsen/mittag-leffler
+ [2] https://github.com/tranqv/Mittag-Leffler-function-and-its-derivative
+ [3] https://se.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function
