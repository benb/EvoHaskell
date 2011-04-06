import System.Environment (getArgs)
import System.Console.GetOpt
import Alignment
import Data.List
import Control.Monad

data Flag = Sort
options = [ Option ['S'] ["sort"] (NoArg Sort) "Sort alignment" ]

main = do args <- getArgs
          aln <- case getOpt RequireOrder options args of 
                     ([],(file:xs),[]) -> parseAlignmentFile parseFasta file
                     ((Sort:xs),(file:ys),[]) -> do myaln <- parseAlignmentFile parseFasta file
                                                    return $ liftM sortAlignment myaln
                     (_,_,msgs) -> error $ concat msgs
          case aln of 
                     Just a  -> putStrLn $ intercalate "\n" $ map (\(x,y) -> x ++ "\t" ++ y) $ (names a) `zip` (sequences a)
                     Nothing -> putStrLn "error"
          return ()
 

