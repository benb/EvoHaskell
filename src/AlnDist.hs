import System.Environment (getArgs,getProgName)
import System.IO
import Alignment
import Tree
import Alignment.Dist
import Text.Printf
import System.Console.GetOpt
import Numeric
import Control.Monad

data Options = Options  { optVersion    :: Bool
                        , optFunc       :: ListAlignment -> ListAlignment -> Either String (Int,Int)
}

options :: [ OptDescr (Options -> IO Options) ]
options = [ Option ['g'] ["gap"] (NoArg (\opt-> return opt {optFunc = safeCompare homGapDist})) "Homolgy distance with labelled gaps",
            Option ['n','s'] ["ssp"] (NoArg (\opt -> return opt {optFunc = safeCompare hom0Dist})) "Symmetrised Sum-Of-Pairs",
            Option ['h'] ["hom"] (NoArg (\opt -> return opt {optFunc = safeCompare homDist})) "Homology distance (default)",
            Option ['t'] ["tree"] (ReqArg (\arg opt -> do treeIO <- (liftM readBiNewickTree) (readFile arg)
                                                          let tree = case treeIO of
                                                                       Right t -> t
                                                                       Left err -> error err
                                                          return $ opt {optFunc = (safeTreeCompare homTreeDist tree)}) "TREE" )
            "Homology distance with tree-labelled gaps"
          ]
safeCompare :: (ListAlignment -> ListAlignment -> (Int,Int)) -> ListAlignment -> ListAlignment -> Either String (Int,Int)
safeCompare dist aln1 aln2 = case compatibleAlignments aln1 aln2 of 
                                        False -> Left "Incompatible alignments"
                                        True -> Right $ dist aln1 aln2

safeTreeCompare dist tree aln1 aln2 = case (compatible tree aln1) of
                                        False -> Left "Tree is incompatible with first alignment"
                                        True -> case (compatible tree aln2) of 
                                                False -> Left "Tree is incompatible with second alignment"
                                                True -> safeCompare (dist tree) aln1 aln2

startOptions = Options {optVersion = False,
                        optFunc = safeCompare homDist }

main = do args <- getArgs
          let (actions, nonOptions, errors)=getOpt Permute options args
          opts <- foldl (>>=) (return startOptions) actions
          let Options {optFunc = f} = opts
          alignments <- mapM readAln nonOptions
          me <- getProgName
          case alignments of 
               (x:y:[]) -> case (f x y) of 
                                   Left err -> hPutStrLn stderr err
                                   Right (d,n) -> putStrLn $ (show n) ++ " / " ++ (show d) ++ " = " ++ (show ((fromIntegral n)/(fromIntegral d)))
               _        -> hPutStrLn stderr (usageInfo me options)



readAln :: String -> IO ListAlignment
readAln x = do rawa <- parseAlignmentFile parseFasta x
               return $ case rawa of 
                          Right aln -> aln
                          Left err -> error err

