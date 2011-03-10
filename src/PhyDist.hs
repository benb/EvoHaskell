import System.Environment (getArgs)
import Alignment
import Tree
import Alignment.Dist
import Text.Printf
import System.Console.GetOpt
import Text.ParserCombinators.Parsec

main = do args <- getArgs
          parseCommand args >>= putStrLn
 
parseCommand [String] -> String
parseCommand ("-g":xs) = diff homGapDist xs
parseCommand ("-t":xs) = diffTree homTreeDist xs
parseCommand xs = diff homDist xs

diffTree dist (x:y:z:xs) = do a <- parseFastaFile x
                              b <- parseFastaFile y
                              return (case (parse parseTree "" (treeStr)) of 
                                      Left err -> show err
                                      Right t -> goTree dist (compatible t a) (compatible t b) t a b -- (printf "%.6f" (dist a b)) ++ (show (enforceBi t)))

diffTree dist x = return "Usage: phydist <fasta1> <fasta2> <tree>"

diff dist (x:y:xs) = do a <- parseFastaFile x
                        b <- parseFastaFile y
                        printf "%.6f" (dist a b)
diff dist x = return "Usage: phydist <fasta1> <fasta2>"


goTree dist false x t a b = "Tree is incompatible with first alignment"
goTree dist true false t a b = "Tree is incompatible with second alignment"
goTree dist true true t a b = printf "%.6f" (dist t a b)
