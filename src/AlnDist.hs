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
parseCommand ("-g":xs) = Main.diff homGapDist xs
parseCommand ("-n":xs) = Main.diff hom0Dist xs
parseCommand ("-t":xs) = diffTree homTreeDist xs
parseCommand xs = Main.diff homDist xs

diffTree :: (Node->ListAlignment->ListAlignment->(Int,Int)) -> [String] -> IO String
diffTree dist (x:y:z:xs) = do rawa <- parseAlignmentFile parseFasta x
                              rawb <- parseAlignmentFile parseFasta y
                              let a = fmap sortAlignment rawa
                              let b = fmap sortAlignment rawb
                              treeStr <- readFile z
                              let t = readBiNewickTree treeStr
                              return $ goTree dist ((liftM2 compatible) t a) ((liftM2 compatible) t b) t a b

diffTree dist x = return usage

usage = "Usage: alndist <-t> <-g> <-n> <fasta1> <fasta2> <tree>"

diff :: (ListAlignment -> ListAlignment -> (Int,Int)) -> [String] -> IO String
diff dist (x:y:xs) = do rawa <- parseAlignmentFile parseFasta x
                        rawb <- parseAlignmentFile parseFasta y
                        let a = fmap sortAlignment rawa
                        let b = fmap sortAlignment rawb
                        let ans = (liftM2 dist) a b
                        case ans of 
                                (Right (numPairs,numDiff)) -> return $ (show numDiff) ++ " / " ++ (show numPairs) ++ " = " ++ (showEFloat Nothing ((fromIntegral numDiff)/(fromIntegral numPairs)) "")
                                (Left err) -> return err

diff dist x = return usage


goTree :: (Node->ListAlignment->ListAlignment->(Int,Int)) -> Either String Bool  -> Either String Bool -> Either String Node -> Either String ListAlignment -> Either String ListAlignment -> String
goTree dist (Right False) x t a b = "Tree is incompatible with first alignment"
goTree dist (Right True) (Right False) t a b = "Tree is incompatible with second alignment"
goTree dist (Right True) (Right True) (Right t) (Right a) (Right b) = (show numDiff) ++ " / " ++ (show numPairs) ++ " = " ++ (showEFloat Nothing ((fromIntegral numDiff)/(fromIntegral numPairs)) "") 
                                                                   where numPairs = snd ans
                                                                         numDiff = fst ans
                                                                         ans = dist t a b 

