
  # list of short code snippets to help in
  # move/reformat of file names and locations using R or bash

  filetofind = "ep364.txt"
  dirtolook = "data/snowcrab/ENS\ Snow Crab 2004 Survey/"
  cmd("find", dirtolook, "-iname", filetofind) 

  cmd("cp ENS\ Snow\ Crab\ 2004\ Survey/*/Minilog/*/*.T* minilog/")


