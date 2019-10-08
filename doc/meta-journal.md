# Some points to to into account while starting this new project

* Keep all data in data sub directory  
  I could even separate data into:
  * incoming raw data
  * incoming external data
  * generated data (including intermediate data)
  * final results data
  The idea being that the 'analyses' directory remains small and can be put under version control in its entirety

* Keep stringent project-journal. Generate html and pdf versions of journal at every git commit

* Create separate top 'code' directory which will contain all scripts, etc. This will make it easy to separate version control for code on one hand and analyses on the other hand.

* Given that data is separate from 'analyses' directories scrips always have to specify output directory

* Generated data always labeled with date tag, and digest hash

* Generate html and pdf docs of every "knitr doc"
