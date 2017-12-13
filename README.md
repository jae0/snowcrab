# bio.snowcrab 

Utilities to help develop and/or use snowcrab assessment tools in the bio framework.

Installation:

```
install.packages( "devtools", ask=F, dependencies=TRUE )   
require ("devtools")
install_github( "jae0/bio.snowcrab" )
```

Setup environment as indicated in: https://github.com/jae0/aegis.env



NOTE to all:

Stock assessment of Canada's Maritimes Region snow crab (Chionoectes oplio) leveraging aegis*, bio*, and stm* packages. 

To install you need to bootstrap from https://github.com/jae0/aegis.env directly: 

```
  devtools::install_github( "jae0/aegis.env" )
```

Then, you need to have an Rprofile set up properly. An example can be seen in aegis.env/R/project.Rprofile.example.R, or use the following, being careful to define the required R-global variables:

```.
libPaths("~/R")
homedir = path.expand("~")
tmpdir = file.path( homedir, "tmp" )
work_root = file.path( homedir, "work" )    ### replace with correct path to work directory (local temporary storage)
code_root = file.path( homedir, "bio" )   ### replace with correct path to the parent directory of your git-projects
data_root = file.path( homedir, "bio.data" )   ### replace with correct path to your data

# store your passwords and login here and make sure they are secure
passwords = file.path( homedir, ".passwords" )
if (file.exists(passwords)) source( passwords )

require( aegis.env ) 
```


Thereafter, you can used the bootstrapped environment to install the other basic tools: 

```
  aegis.env::project.libraryInstall(DS="snowcrab")
```

If you have a local git clone of the required packages, you can install with:

```
  aegis.env::project.libraryInstall(DS="snowcrab", local=TRUE)  

```

For usage, examples can be found in https://github.com/jae0/aegis/inst/scripts and https://github.com/jae0/bio.snowcrab/inst/scripts. 

