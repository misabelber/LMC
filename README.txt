This repository stores the codes, tools and models for the analysis of Limits on Dark Matter in the Large Magellanic Cloud with CTA.

==========CONFIGURATION===========================================================================

--Pipelines in Python/Ctools:

	    Activate ctools environment:
	    	     - source activate ctools	    

	    To run the pipelines in python/ctools which are mainly for simulate data and create 
	    models it's necessary to configure the right paths in the file config.sh . 

--Pipelines in C++/Root:
	    
	    Pipelines for C++/Root requires the repository "Math": https://github.com/misabelber/Math 
	    
	    To run the pipelines in C++/Root which are for likelihood fitting and upper limits calculation it's necessary to configure the right paths in the file config.sh. 
	    
	    	    
==========DEPENDENCIES============================================================================
	Ctools 1.5
	Python 2.7
	Root 5.34
	"Math" repository	    
	
=========FOLDERS DESCRIPTION==================================================

-- models:

	This folder contains .xml models for all the simulated sources, both
	baryonic and Dark Matter.

-- pipelines:
	In this folder the principal codes for the analysis are stored.
	- pipes_in_py: Contain pipelines in python/ctools that produce
	the data for the analysis. They simulate observations of the
	sources, bin the results and build models of predicted counts.

	. pipes_in_C: Here there are the codes in C++/ROOT to perform the
	likelihood analysis for the upper limits calculation and the study
	of source correlations.

-- results:

	This is the place where upper limits results (in form of textfiles)
	from pipelines/pipes_in_C/LMC_repository.C will be stored.

-- spectra:

	In this folder textfiles with spectra information of the sources
	are stored.

-- tools:

	In this folder there are codes for diverse utilities, from producing
	spectra or .xml models to plotting macros.

	- tools_in_py: This folder contains codes in python/ctools for
	production of spectra, .xml models and model edition. Here
	the Cirelli et al. DM spectra can be managed.

	- tools_in_C: Here there are some C++/ROOT scripts for plotting
	upper limits results.  	      	   	 


============USAGE INSTRUCTIONS===========================================

-- Preparing the data:

   - Use /tools/tools_in_py/producetxtfileW.py to extract DM spectra from
   Cirelly et al. "AtProduction_gammas.dat" file. Edit wanted DM masses
   or final state.

   - Use /tools/tools_in_py/PSmodel_maker.py to producle .xml model files
   for point sources.

   - Use /tools/tools_in_py/DMmodel_maker.py to produce .xml model files
   for DM realizations.

   - Use and modify /tools/tools_in_py/model_editor.py to make changes
   in .xml model files.

-- Simulating the data and model producing:

   - Run pipelines/pipes_in_py/observations_KSP.py to simulate
   observations of wanted sources following the KSP pointing pattern.
   It also builds a model of predicted counts obtained from the .xml
   model + simulation result.

   - Run pipelines/pipes_in_py/producedataDM_KSP.py to simulate
   observations of Dark Matter of different masses at once following the
   KSP pointing pattern.It also builds a model of predicted counts.

   - Run pipelines/pipes_in_py/producedataPS.py to simulate observations
   of the baryonic point sources all at once following the KSP pointing
   pattern. It also builds a model of predicted counts.


-- Likelihood Analysis:

   - Before using, remove any existing *_C* file in /pipelines/pipes_in_C/

   - Inside a ROOT interpreter run .x /pipelines/pipes_in_C/load.C to
     load the functions for analysis contained in
     /pipelines/pipes_in_C/LMCrepository.C.

   - Edit the pertinent global variables from /pipelines/pipes_in_C/
     LMCrepository.C such as N (number of Baryonic Background), dm_mass,
     names of .fits files, etc. 
	
   - From /pipelines/pipes_in_C/LMCrepository.C run "Check_Correlations"
     function to study the correlation factors between Dark Matter and the
     other Baryonic Components and to be sure that the maximum Likelihood
     is properly found.
     Arguments are: firstbin,lastbin,dm_mass (in TeV). Usually firstbin=0,
     lastbin=20.
      
   - From /pipelines/pipes_in_C/LMCrepository.C run "Bands()" to calculate
   upper limits on dark matter for given dark matter mass. It runs the
   calculation over different realizations of the data (default is 10,
   ideal is > 100). so in the end a mean and rms can be computed for the
   upper limit with containment bands.
   Arguments are: dm_mass(in TeV), bestcase (Boolean,
   True for a supposed very well known background so the baryonic
   components will be forced to fit the model with a given tolerance.
   False to let every source free in the fitting so the correlations
   will affect more the DM upper limit result), tol (tolerance level
   for best case, default is tol=0.01 which mean that baryonic sources
   normalization won't differ more than 0.01 from 1).
   Resulting upper limits are stored in results/ as a textfile.

-- Results visualization:

   - Run /tools/tools_in_C/plotbands.C as ROOT macro to plot the DM
   upper limit mean and containment bands (ay 95% and 84% CL) obtained
   from /pipelines/pipes_in_C/LMCrepository.C "Bands()".
     
