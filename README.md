# ElectromechanicalPipeline
This repository will house the electromechanical data analysis pipeline developed during my PhD. The purpose is simply to showcase data analysis techniques.

The files in this repo are divided into three types:
1) Pig Data files, including subsections of the raw data for analysis purposes. These may be electrophysiological data (either .BDF or .MAT) and impedance planimetry data (.CSV).
2) Correlation Functions, which perform numerous analysis tasks.
3) Analysis scripts, which read in appropriate data and call appropriate correlation functions in turn.

These three data files are all required for meaninful data outputs. For licensing reasons, the raw contractile and electrophysiological data has been ignored from this repository.
Importantly, electrophysiological data has been prepared in the Gastrointestinal Electrical Mapping Suite (GEMS), to which I do not have the rights to share, and is formatted uniquely. Additionally, some of the correlation functions invoke GEMS functions, which will once again be unavailable in this repository.
