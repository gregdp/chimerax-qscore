Implementation of the map-model Q-Score algorithm (https://www.nature.com/articles/s41592-020-0731-1) in ChimeraX.

To install
* download latest version .whl file from dist/ folder
* in ChimeraX command line: toolshed install [path to .whl file]

Example usage 1:
* From Tools -> Validation -> Model-map Q-score
* Select model, map, press Calculate
* In Plot display, use mouse wheel to zoom in/out; click on point to display a residue and the map around it
* Select Ribbon: Residue-Average to color ribbon with Q-scores

Example usage 2:
* To generate a file with per-residue Q-scores:
* In command line: "qscore #1 toVolume #2 useGui False mapResolution 3.1"
* #1 refers to residues to calculate Q-scores for
* #2 refers to volume
* mapResolution is used to generate a 'expected Q-score' baseline to compare resulting Q-scores to

Example usage 3:
* To show residues/atoms and densities:
* In command line: "q #1:100-110 showResInMap #6 distNear 4"
* #1:100-110 refers to residues in model #1 to show
* showResInMap #6 is where to take densities from
* (optional) distNear 4 takes otehr residues with atoms within 4Å

Example maps and models:
* Streptavidin: 1.7Å -- emd:31083, pdb:7efc -- "open 31083 from emdb; open 7efc from pdb"
* Sars2 Orf3a: 2.9Å -- emd: 22136 pdb: 6xdc -- 
* Ribozyme: 3.1Å -- emd:31385 pdb:7ez0 --

