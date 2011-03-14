import System.Environment (getArgs)
import Alignment
import Tree
import Alignment.Dist
import Text.Printf
import System.Console.GetOpt
import Numeric
import Control.Monad

main = do args <- getArgs
          ans <- parseCommand args 
          putStrLn ans
          return ()
 
parseCommand :: [String] -> IO String
parseCommand ("-g":xs) = diff homGapDist xs
parseCommand ("-t":xs) = diffTree homTreeDist xs
parseCommand xs = diff homDist xs

diffTree dist (x:y:z:xs) = do a <- parseAlignmentFile parseFasta x
                              b <- parseAlignmentFile parseFasta y
                              treeStr <- readFile z
                              let t = readBiNewickTree treeStr
                              return $ goTree dist ((liftM2 compatible) t a) ((liftM2 compatible) t b) t a b

diffTree dist x = return "Usage: phydist <fasta1> <fasta2> <tree>"

diff :: (ListAlignment -> ListAlignment -> Double) -> [String] -> IO String
diff dist (x:y:xs) = do a <- parseAlignmentFile parseFasta x
                        b <- parseAlignmentFile parseFasta y
                        let ans = (liftM2 dist) a b
                        case ans of 
                                (Right output) -> return $ showEFloat Nothing output ""
                                (Left err) -> return err

diff dist x = return "Usage: phydist <fasta1> <fasta2>"


goTree dist (Right False) x t a b = "Tree is incompatible with first alignment"
goTree dist (Right True) (Right False) t a b = "Tree is incompatible with second alignment"
goTree dist (Right True) (Right True) (Right t) (Right a) (Right b) = showEFloat Nothing (dist t a b) ""
