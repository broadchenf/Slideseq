
args = commandArgs(trailingOnly=TRUE)


sayHello <- function(){
  print('hello')
  print(args[[1]])
}

sayHello()