fn = "/GitHub/str2phy/ex/ferritins_M20/tree_comparison/output_filenames_v1.txt"
lines = readLines(fn, skipNul=TRUE)

TF = startsWith(x=lines, prefix="#")
lines = lines[TF==FALSE]

TF = lines == ""
lines = lines[TF==FALSE]

length(lines)