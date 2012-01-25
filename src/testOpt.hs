import Phylo.Opt
import Phylo.NLOpt
testFunc [a,b] = (a*a)+(b*b)

main = do ans <- bobyqa [0.01,0.01] 1E-04 [4.9,-0.4] testFunc [Just (-5.0),Just (-5.0)] [Just 5.0,Just 5.0]
          print ans


