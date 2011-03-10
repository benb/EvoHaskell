import System.Environment (getArgs)
import Alignment
import Tree
import Alignment.Dist
import Text.Printf
import System.Console.GetOpt
import Numeric

main = do args <- getArgs
          ans <- parseCommand args 
          putStrLn ans
          return ()
 
parseCommand :: [String] -> IO String
parseCommand ("-g":xs) = diff homGapDist xs
parseCommand ("-t":xs) = diffTree homTreeDist xs
parseCommand xs = diff homDist xs

diffTree dist (x:y:z:xs) = do a <- parseFastaFile x
                              b <- parseFastaFile y
                              treeStr <- readFile z
                              return (case (readBiNewickTree treeStr) of 
                                      Left err -> show err
                                      Right t -> goTree dist (compatible t a) (compatible t b) t a b) where

diffTree dist x = return "Usage: phydist <fasta1> <fasta2> <tree>"

diff dist (x:y:xs) = do a <- parseFastaFile x
                        b <- parseFastaFile y
                        return $ showEFloat Nothing (dist a b) ""
diff dist x = return "Usage: phydist <fasta1> <fasta2>"


goTree dist False x t a b = "Tree is incompatible with first alignment"
goTree dist True False t a b = "Tree is incompatible with second alignment"
goTree dist True True t a b = showEFloat Nothing (dist t a b) ""
