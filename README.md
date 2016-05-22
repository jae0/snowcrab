
NOTE to all: 
  
The year in the file "current.assessment.year.r" must be changed every year.
It must be kept separate from "initialize.local.environment" as running in parallel mode 
requires occasionally overrding some parameters in "p". This override cannot be completed as  
"initialize.local.environment.r" is sourced with every initialisation of a a new CPU.

Begin in sequence due to data dependencies starting with 1.snowcrab.r, 2. ....
   


    

