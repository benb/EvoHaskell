import System.Environment (getArgs)
import Alignment
import Tree
import Alignment.Dist
import Text.Printf
import System.Console.GetOpt
import Text.ParserCombinators.Parsec

main = do args <- getArgs
          let (dist,files) = if ((head args)=="-g") 
                then (homGapDist,(drop 1 args))
                else (homDist,args)
          (diff dist files) >>= putStrLn

diff dist (x:y:z:xs) = do a <- parseFastaFile x
                          b <- parseFastaFile y
                          treeStr <- readFile z
                          return (case (parse parseTree "" (treeStr)) of 
                                Left err -> show err
                                Right t -> (printf "%.6f" (dist a b)) ++ (show (enforceBi t)))
diff dist x = return "Usage: phydist <fasta1> <fasta2> <tree>"
                        
       
