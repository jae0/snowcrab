
export.to.clipboard = function (results) {
  write.table(results, "clipboard", sep="\t")
  return ("Now you can paste the table into excel, or word, etc")
}


