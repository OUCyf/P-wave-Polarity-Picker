# P-wave-Polarity-Picker
A GUI for picking first motion of P-wave and map into beachball



# How to use P-wave Polarity Picker

1. Install Dependencies

	- `Install anaconda or miniconda` 
	- `conda install numpy`
	- `conda install matplotlib`
	- `conda install obspy -c conda-forge`


2. Generate Velocity Model List via taup Method

	- v_model.nd is the velocity model for this earthquake event, and you can also set your own velocity model for your own events, according to https://docs.obspy.org/packages/obspy.taup.html

	- In your terminal: `python build_v_model.py`
	
    - This command will generated the `v_model.npz` file, which will be used in Picker GUI.


3. Pick Now

	- In your terminal: `python picker_gui.py`
	- Load data
	- Set the strike, dip, rake, and v_model value. For this earthquake, the strike is 291, dip is 85, rake is -1, and the v_model is v_model.npz which is generated in step-2. I also provide 2 general v_model which do not need v_model.npz file, just set the v_model to prem or iasp91.
	- Click plot button
	- Now you can pick the p-wave polarity, the mouse left click is up (red color), right click is down (black color). You can also use scroll down/up to zoom in/out figure. Delete button will delete the latest pick point. XY_init button will set the figure to the original size.
	- After picking, you need save those infos into the file. Click the save button, you will save to `./picks.txt` file


4. Now try your case.



# License

MIT
	 