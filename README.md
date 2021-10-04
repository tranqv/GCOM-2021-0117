# GCOM-2021-0117
**On mild solutions of the Laplacian fractional Langevin equations with anti-periodic type boundary conditions**

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

**3) Run the test cases**

  a) Download this package 
  
  https://codeload.github.com/tranqv/GCOM-2021-0117/zip/refs/heads/main
  
  and extract it. One should have the folder **GCOM-2021-0117-main**.
  
  Inside the subdirectory **code** of the folder **GCOM-2021-0117-main**, three examples are prepared with their main programs as follows:
     
    + main_ex01a.py  Example 1a,  
    + main_ex01b.py  Example 1b, 
    + main_ex02a.py  Example 2,  
    + main_ex03a.py  Example 3a, without noise,
    + main_ex03b.py  Example 3b, with noise to the input. 

  All programs perform Procedure (P) of the manuscript. 
  
  Please read comments inside the codes to gain informations and check with output to the screen.
  
  We can use **main.py** to check if the codes works properly or not by: 
  + setting the value of **lmax** moderately, e.g. **lmax = 5**, and 
  + assigning the value of **fprefix** as it was prepared inside the code.
  
  b) Open a Command Prompt (Windows) or Terminal (Linux)
  
  c) Change directory to **GCOM-2021-0117-main**, and then the **code** inside 
  
  d) Inside the folder **code**, run the examples as the following commands 
  
      python3 main_ex01a.py
      python3 main_ex01b.py
      python3 main_ex02a.py
      python3 main_ex03a.py
      python3 main_ex03b.py
    
  where the first command is for Example 1a, and so on.
    
  Notes:
  + At the first time we run one of the above commands, two folders, i.e., **figs** and **xout**, will be created.
    + **figs** saves figure output from all the main program
    + **xout** saves the results in ascii format.
  + The **main_ex03a.py** must be excecuted in advance of **main_ex03b.py** since the latter program needs the output 
  from the preceeding one, while the others (for Examples 1a, 1b, 2, and 3a) can run independently.
      
**4) Postprocessing**
   
  After running the examples, we can use the command 
  
      python3 main_post.py
      
  to plot the results of Examples 1a, 1b, 2, 3a. 
  
  Note that, within **main_post.py** we must set the values of **fprefix** and **lmax** accordingly.
  While figures for Example 3b were made when we run the command **python3 main_ex03b.py** as above.
    

**References:**

+ [1] https://github.com/khinsen/mittag-leffler
+ [2] https://github.com/tranqv/Mittag-Leffler-function-and-its-derivative
+ [3] https://se.mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function
