.First.lib = function(libname, pkgname)
{
  cat('| \n')
  cat('| Rgtsp - generalized top scoring pairs in R\n')
  cat('|\n')
  cat('| Loading...')

    s = search()
    library.dynam('Rgtsp', pkgname, libname, now=FALSE)
    
  cat('OK\n')
}

