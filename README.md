# Hodgkin-Huxley-MatLab-Repo
In this repo we endeavor to simulate Hodgkin-Huxley electrical dynamics on neurons.

## Software Requirements
* Most recent version of MatLab, [MatLab](https://www.mathworks.com/products/matlab.html)
* Recommended (for Windows Users): gitforwindows, [gitwindows](https://gitforwindows.org/)
* MacOS users using the bash is fine 

## Example Usage
The following steps have been tested on Windows 10 (will update when tested on Mac and Linux)
To run the code:
1. First open a  bash/terminal (whatever!) on your desktop (or wherever!) and
execute the following <code>git clone https://github.com/jarosado0911/Hodgkin-Huxley-MatLab-Repo.git</code> or you can download the ZIP of the repo. 
Below is a picture of what this may look like 

![gitclone](images/gitclone.PNG)

2. Next, open MatLab and navigate into the <code>Hodgkin-Huxley-MatLab-Repo</code>. It should look like this:

![matlab](images/matlab.PNG)

3. Next, execute in the MatLab command window <code>addpath('simulation_core');</code> this will tell MatLab to look inside the <code>simulation_core</code> folder for the functions we wish to call! Nothing will happen when you do this!
4. Now, let us try to run a simulation. 
  - Execute <code>neuronSimTest('sample_geometry',[0],'C:/path_to_your_Desktop/Desktop/sampleTest',1);</code>
  - This will tell MatLab to look inside the <code>sample_geometry</code> folder for the <code>.swc</code> geometry file, 
  - it will run the simulation on 0 refinement, 
  - the output will be saved to the desktop in <code>sampleTest</code>, 
  - and <code>1</code> will save the voltage state of the entire cell at every time step, set it to <code>0</code> if you don't want to save voltage data.
As the simulation is running you will see this:

<p float="left">
  <img src="images/running.PNG" alt="while running" />
  <img src="images/finished.PNG" alt="when complete" /> 
</p>

The following output will be in the output folder you specified:
 - Inside <code>sampleTest</code> there is a subfolder called <code>0Ref2.0000000</code>, this folder contains all the data from the simulation
 - Inside <code>0Ref2.0000000</code> there is a subfolder called <code>data</code>, this contains all the voltage data at every time step, if you set the simulation to <code>0</code> then this folder will be empty.
 - Inside <code>0Ref2.0000000</code> there are <code>.dat</code> files, these are the voltage, m,n,h state data for only the BRANCH points of the cell. These files will always be saved!
 - Inside <code>0Ref2.0000000</code> there are two <code>.mat</code> files corresponding to the voltage at the soma and the time values.

## Making a Movie

5. Picking up from the example above execute in MatLab

<code>makematlabmov('C:/path_to_your_Desktop/Desktop/sampleTest/0Ref2.0000000','sample_geometry/228-13MG.CNG_segLength=8_1d_ref_0.swc','C:/path_to_your_Desktop/Desktop/sample');</code>

  - where <code>path_to_your_Desktop</code> is the file path to where your Desktop is.
  - <code>'C:/path_to_your_Desktop/Desktop/sampleTest/0Ref2.0000000'</code> this is the folder with all the data from your simulation
  - <code>'sample_geometry/228-13MG.CNG_segLength=8_1d_ref_0.swc'</code> this is the precise <code>.swc</code> geometry file that you used in your simulation run.
  - <code>'C:/path_to_your_Desktop/Desktop/sample'</code> this is the path to where you want to save the video and the name, in this case it will be named <code>sample.mp4</code> and will be saved on the desktop.

6. As it is running the MatLab output will show <code>frame = <some_numbers></code> and a window will open showing the plot, like below

![](images/sample.gif)
