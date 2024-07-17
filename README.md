This program, w2kernel, converts weight function in frequency to     
kernel function in time domain.                                      
w2kernel works interactvely, just start the linux binary file in the bin directory, then        
1. you need to give the number of beads in the PIMD simulation       
2. the temperature of the simulation must be entered                 
3. the program asks for the timestep between two snapshots           
The weight function is expected in the file of "gx_function.dat".    
The weight function should be given from 0 to 11999.95 with 0.05     
step size. The weight function can be                                
 a) determined from QTB simulation in CP2K                           
 b) determined from a fortran code:                                  
    https://github.com/madaraszadam/PI_weight_functions              
    https://doi.org/10.26434/chemrxiv-2024-x36vm-v2                  
 c) downloaded from this site:                                       
    https://zenodo.org/records/10702413                              
The resulting kernel function is written in the file of              
"kernel.dat"                                                         
The algorithm is based on Fourier transformation, more details can be         
found in the following paper:                                        
Dénes Berta, Dávid Ferenc, Imre Bakó and Ádám Madarász               
"Nuclear Quantum Effects from the Analysis of Smoothed Trajectories: 
Pilot Study for Water"                                               
https://doi.org/10.1021/acs.jctc.9b00703

In the example directory you can run the program with this command:
"chmod +x w2kernel; ./w2kernel < input.txt"

If you want to compile the source code, then you need the install fftw3 from https://www.fftw.org/
If the fftw3 is installed in the default directorioes, then you can compile the source code with gfortran with the next command:
"gfortran w2kernel.f -I /usr/local/include -L /usr/local/lib/ -lfftw3"
                             
email: madarasz.adam@ttk.hu                                          
