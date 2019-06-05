
### shinyApps.io installing minfiShinyl ###

library(rsconnect)

bioc <- local({
  env <- new.env()
  on.exit(rm(env))
  evalq(source("http://bioconductor.org/biocLite.R", local=T), env)
  biocinstallRepos()
})

options(repos = BiocManager::repositories())


rsconnect::setAccountInfo(name='margarciavalero',
                          token='X',
                          secret='X')


rsconnect::deployApp('C:/Users/alumne/Desktop/TFM/data/shinyApp')


