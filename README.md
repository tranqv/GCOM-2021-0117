# GCOM-2021-0117
On mild solutions of the Laplacian fractional Langevin equations with anti-periodic type boundary conditions

1) Install Anaconda for Python 3 

  a) Windows: 
    https://repo.anaconda.com/archive/Anaconda3-2021.05-Windows-x86_64.exe
    Guide: 
      https://docs.anaconda.com/anaconda/install/windows/
      https://www.datacamp.com/community/tutorials/installing-anaconda-windows

  b) Linux: 
    https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
    Guide:     
      https://docs.anaconda.com/anaconda/install/linux/

2) Install Python libraries 

  a) Open a Command Prompt (Windows) or Terminal (Linux)
  
  b) Run the following commands
  
    conda install numpy
    conda install scipy
    conda install matplotlib

3) Run the test cases

  a) Download this package ans unzip it. One should have the folder **GCOM-2021-0117**
  
     Inside the subdirectory **code** of the folder **GCOM-2021-0117**, three examples are prepared with their main 
     programs as follows:
     
    + main_ex01a.py  Example 1a,  
    + main_ex01b.py  Example 1b, 
    + main_ex02a.py  Example 2,  
    + main_ex03a.py  Example 3a, 
    + main_ex03b.py  Example 3b. 
    
    One can use main.py to check if the codes works properly or not by 
      + setting the value of "lmax" not so large, e.g. lmax=5, and 
      + assigning the value of "fprefix" as it was prepared inside the code.
  
  b) Open a Command Prompt (Windows) or Terminal (Linux)
  
  c) Change directory to "GCOM-2021-0117" 
  
  d) Run the examples as the following commands 
  
      python3 main_ex01a.py
      python3 main_ex01b.py
      python3 main_ex02a.py
      python3 main_ex03a.py
      python3 main_ex03b.py
    
    where the first command is for Example 1a and so on.
    
    Note that the "main_ex03a.py" must be run in advance of "main_ex03b.py" since the latter program needs the output 
    from the preceeding one, while the others can run independently.
      
4) Postprocessing 
   
    After running the examples, we can use the command "python3 main_post.py" to plot the results of Examples 1a, 1b, 2, 3a. 
    Note that, we must set the values of "fprefix" and "lmax" accordingly.
    While figures for Example 3b were made when we run the command "python3 main_ex03b.py" as above.
    
    
    
