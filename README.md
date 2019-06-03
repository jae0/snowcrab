Stock assessment of Canada's Maritimes Region snow crab (Chionoectes oplio) leveraging aegis*, bio*, and stm* packages.


Installation:


1. To install:

```
  install.packages( "devtools", ask=F, dependencies=TRUE ) # to inter-operate with github
  devtools::install_github( "jae0/aegis" ) # to bootstrap by installing directly from github
  aegis::project.libraryInstall(DS="snowcrab") # install bio.snowcrab and other required packages
```



2. Then, you need to have an Rprofile set up properly. Use the following, being careful to define the required R-global variables (see also: https://github.com/jae0/aegis/src/master/R/project.Rprofile.example.r):


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

require( aegis )
```


If you have a local git clone of the required packages, you can install with:

```
  aegis::project.libraryInstall(DS="snowcrab", local=TRUE)

```

For usage, examples can be found in aegis.*.
