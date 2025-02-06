.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    crayon::green("Welcome to RNAmeta! This package is currently under active development.")
  )
  packageStartupMessage(
    crayon::yellow("You can visit our interactive online website:\nhttp://rnainformatics.cn:3838/RNAmeta")
  )
  packageStartupMessage(
    crayon::cyan("or visit our GitHub page:\nhttps://github.com/Rwtkc/RNAmeta/tree/master")
  )
  packageStartupMessage(
    crayon::blue("If you have any questions, please contact the developer @Rwtkc: rwtkc@foxmail.com")
  )
}
